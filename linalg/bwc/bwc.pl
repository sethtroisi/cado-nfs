#!/usr/bin/env perl
use warnings;
use strict;

use POSIX qw/getcwd/;
use File::Basename;
use File::Temp qw/tempdir/;
use List::Util qw[max];
use Data::Dumper;
use Fcntl;

# This companion program serves as a helper for running bwc programs.

# MPI_BINDIR=/opt/mpich2-1.1a2-1.fc10.x86_64/usr/bin ./bwc.pl :complete matrix=c72 interval=256 splits=64 mn=64 ys=0..64 mpi=4x4

# {{{ usage
sub usage {
    print STDERR <<EOF;
Usage: ./bwc.pl <action> <parameters>

# 20140901: some revamping in progress. The instructions below may be
# partly inaccurate.

The action to be performed is either:
- the path to a program, in which case the command line eventually called
  is quite similar to ``<action> <parameters>'', with some substitutions
  performed. Depending on the program name, some options are discarded
  because they are not meaningful.
- one of the special actions defined by the script. These actions are
  prepended by a colon.
    :complete -- do a complete linear system solving. The solution ends
                 up in \$wdir/W.
    :wipeout  -- quite like rm -rf \$wdir, but focuses on bwc files. Does
                 not wipe out matrix file, nor cached in-memory files.
    :bench    -- does :wipeout, dispatch, prep, secure, and krylov with
                 an exceedingly large finish bound in order to perform
                 some timings.

Parameters are specified in the form <key>=<value>.
All parameters having a meaning for the bwc programs are accepted. Some
parameters control the wrapping script itself.

mpi thr              same meaning as for bwc programs.
mn m n prime         same meaning as for bwc programs.
nullspace            same meaning as for bwc programs.
interval             same meaning as for bwc programs.
ys splits            same meaning as for bwc programs (krylov and mksol)

matrix               input matrix file. Must be a complete path. Binary, 32-b LE
rhs                  rhs file for inhomogeneous systems (in binary, 32-b LE)
nrhs                 number of rhs columns
wdir                 working directory.
                     If unset, derived from 'matrix', 'mpi' 'thr'
mpiexec              Command to invoke for mpiexec
hosts                Comma-separated list of hosts. 'hosts' may be supplied
                     several times, so that lists concatenate. Thus one
                     may use the shell to supply hosts=node{1,3,5,7}.
hostfile             path to a host list file to be understood by mpiexec

Amongst parameters, some options with leading dashes may be set:

-d  Do a dry run ; show the commands that would be executed.
    (-d is relevant to the script only)
-v  Increase verbosity level ; understood by bwc, but meaningless at the moment
-h  Show this help.

Some environment variables have an impact on this script.

MPI_BINDIR   The path to the directory holding the mpiexec program (the
             mpiexec= parameter allows to specify the same thing. Note that
             MPI alone serves as an alias for MPI_BINDIR, but this may
             become legacy. Normally mpiexec is obtained by cmake
             substitution, at least this is the default.).
BWC_BINDIR   The directory holding the bwc binaries. Defaults to the
             directory holding the current script.

EOF

    if (scalar @_) {
        print STDERR "Error messages:\n", join("\n", @_), "\n";
    }

    exit 1;
}
# }}}

# $program denotes either a simple command, or something prepended by ':'
# which indicates something having a special meaning for the script.
my $main = shift @ARGV or usage;
print "$0 $main @ARGV\n";

# ----- cmake substituted variables -----
## mpiexec is substituted by cmake in case mpi has been used for the
## compilation. NOTE that this means that a priori, mpiexec _must_ be
## used for running all programs.
my $mpiexec='@MPIEXEC@';


my $bindir;
if (!defined($bindir=$ENV{'BWC_BINDIR'})) {
    $bindir=$0;
    $bindir =~ s{/[^/]*$}{};
    if ($bindir eq $0) {
        $bindir='.';
    }
}


##################################################
# Parse the command line, and fill the following:
#  - %$param hash, with special rules for the hosts argument.
#  - $show_only, @extra_args
#  - obey -h
my @extra_args;
my $show_only=0;
my $param={};
my $param_defaults={
    # some defaults.
    prime => 2,         # for factoring.
    thr => '1x1',
    mpi => '1x1',
    shuffled_product => 1,
    interval => 200,
};
# {{{
while (defined($_ = shift @ARGV)) {
    if ($_ eq '--') {
        push @extra_args, splice @ARGV;
        last;
    }
    # -d will never be found as a bw argument, it is only relevant here.
    if (/^(-d|--show|--dry-run)$/) { $show_only=1; next; }
    if (/^(-h|--help)$/) { usage; }
    my ($k,$v);
    if (/^--([^=]+)$/ && scalar @ARGV) { $k=$1; $v=shift(@ARGV); }
    elsif (/^(-.*)$/) { $k=$1; $v=undef; }
    elsif (/^([^=]+)=(.*)$/) { $k=$1; $v=$2; }
    if (!defined($k)) {
        usage "Garbage not understood on command line: $_";
    }
    if ($k eq 'bwc_bindir') { $bindir=$v; next; }
    if (!defined($param->{$k})) {
        $param->{$k}=$v;
        next;
    }
    if ($k eq 'hosts') {
        if (ref $param->{$k} eq '') {
            $param->{$k} = [$param->{$k}, $v];
        } elsif (ref $param->{$k} eq 'ARRAY') {
            $param->{$k} = [@{$param->{$k}}, $v];
        } else {
            die "\$param->{$k} has gone nuts.\n";
        }
    } else {
        usage "parameter $k may not be specified more than once";
    }
}
for my $k (keys %$param_defaults) {
    next if exists $param->{$k};
    $param->{$k} = $param_defaults->{$k};
}

# }}}

# {{{ global variables. Much resides in $param, but we do have globals
# too.
my @main_args;

my @mpi_split=(1,1);
my @thr_split=(1,1);
my $nh;
my $nv;

my $wdir;
my $matrix;
my $rhs;
my $nrhs;
my $balancing;
my $balancing_hash;

## The mpi_extra_args argument is used to pass information to the mpiexec 
## command. The idea is that mpiexec, the mpi driver program, may need
## additional info to properly setup communications between jobs. As an 
## example, if libvirt-bin is installed on a linux system, openmpi will
## need the following additions.
##    
## ./build/x86_64/linalg/bwc/bwc.pl :complete matrix=/net/tiramisu/localdisk/thome/mats/c90b wdir=/local/rsa768/tmp/c59 mn=64 mpi=2x1 thr=2x4 hosts=patate,tiramisu interval=10 mpi_extra_args='--mca btl_tcp_if_exclude lo,virbr0'
##    
## This tells mpi not to try routing traffic through either the lo or the 
## virbr0 interface. For the former, it's already openmpi's default
## behaviour (unless /etc/hosts is screwed up, which this file does tend to
## be sometimes). For virbr0, it's again some local mess, but each machine 
## having this class C network defined, openmpi can't really tell whether 
## they're connected or not -- and of course they're not.
my $mpi_extra_args;
 
my ($m, $n);
my $prime;
my $splitwidth;
my @splits=();

my $force_complete;

my $hostfile;
my @hosts=();
my $mpi;
my $mpi_ver;
my $needs_mpd;
# }}}
#
# Some more command line interpretation. Here we build @main_args, taking
# out the arguments which are only of interest to us.
# Here we read (from $param) the
# following things.
#  - $hostfile $mpi_extra_args $force_complete
#  - @hosts $m $n @splits
#  - @mpi_split @thr_split
#  - $nh $nv
#    --> dimensions of the split of the matrix (horizontal chunks,
#    vertical chunks).
#  - $prime
#    --> Linear system is over GF(p) for that prime. For factoring, which
#    is the default, we have prime==2
#  - $splitwidth
#    --> Number of vector which are put together in a file. This is only
#    inferred from $prime, and set to 64 for prime==2, 1 otherwise.
#
# {{{
sub set_mpithr_param { # {{{ utility
    my $v = shift @_;
    if ($v=~/^(\d+)$/) {
        my $nthreads = $1;
        my $s = int(sqrt($nthreads));
        for (; $s >= 1; $s--) {
            $nthreads % $s or return ($s, $nthreads / $s);
        }
        die "No possible thread split?";
    }
    $v=~/(\d+)x(\d+)$/ or usage "bad splitting value '$v'";
    return ($1, $2);
} # }}}

while (my ($k,$v) = each %$param) {
    # Filter out parameters which are _not_ relevant to the bwc programs.
    if ($k eq 'matrix') { $matrix=$v; next; }
    if ($k eq 'hostfile') { $hostfile=$v; next; }
    if ($k eq 'mpi_extra_args') { $mpi_extra_args=$v; next; }
    if ($k eq 'force_complete') { $force_complete=$v; next; }
    if ($k eq 'hosts') {
        $v=[$v] if (ref $v eq '');
        for (@$v) { push @hosts, split(',',$_); }
        next;
    }
    # The rest is passed to subprograms, unless explicitly discarded on a
    # per-program basis.
    if (defined($v)) {
        push @main_args, "$k=$v";
    } else {
        # This is for switches (like -v) which don't have a "value"
        # associated with them
        push @main_args, "$k";
    }

    # Yet it does make sense to read some parameters, as we are
    # interested by their value.

    if ($k eq 'wdir') { $wdir=$v; next; }
    if ($k eq 'mn') { $m=$n=$v; next; }
    if ($k eq 'n') { $n=$v; next; }
    if ($k eq 'm') { $m=$v; next; }
    if ($k eq 'prime') { $prime=$v; next; }
    if ($k eq 'rhs') { $rhs=$v; next; }
    if ($k eq 'nrhs') { $nrhs=$v; next; }
    if ($k eq 'splits') { @splits=split ',', $v; next; }
    if ($k eq 'mpi') { @mpi_split = set_mpithr_param $v; $param->{$k} = join "x", @mpi_split;}
    if ($k eq 'thr') { @thr_split = set_mpithr_param $v; $param->{$k} = join "x", @mpi_split;}

    # Ok, probably there's a fancy argument we don't care about.
}
$nh = $mpi_split[0] * $thr_split[0];
$nv = $mpi_split[1] * $thr_split[1];
$splitwidth = ($prime == 2) ? 64 : 1;
# }}}

# {{{ Some important argument checks
{
    my @miss = grep { !defined $param->{$_}; } (qw/m n prime matrix interval/);
    die "Missing argument(s): @miss" if @miss;
}

$param->{'matrix'} =~ s/~/$ENV{HOME}/g;
die "$param->{'matrix'}: $!" unless -f $param->{'matrix'};
if ($prime == 2 && (($m % 64 != 0) || ($n % 64 != 0))) {
    die "Currently for p=2 bwc supports only block sizes which are multiples of 64";
}
if (($m % $splitwidth != 0) || ($n % $splitwidth != 0)) {
    die "Currently bwc supports only block sizes which are multiples of $splitwidth";
}
if ($prime == 2 && (defined($rhs) || defined($nrhs))) {
    die "inhomogeneous systems currently supported only for p>2";
}
if (defined($nrhs) && !defined($rhs)) {
    die "nrhs may only be specified together with rhs";
}
if (defined($rhs) && !defined($nrhs)) {
    # Try to read the first line.
    my $header = eval { open F, $rhs or die "$rhs: $!"; <F>; };
    if ($header !~ /^(\d+) (\d+)$/) {
        die "$rhs does not seem to be ascii with header, we can't proceed.\n";
    }
    $param->{'nrhs'} = $nrhs = $2;
}
if (defined($nrhs) && $nrhs > $n) {
    die "nrhs > n is not supported";
}


# }}}

# {{{ wdir preliminary (pre-mpi) checks.
#  - Check that wdir is defined (or define it to some default);
#  - Check that it exists if it has to exist (:complete is entitled to
#  create $wdir

# Selection of the working directory -- there are three relevant options.
# The 'matrix' option gives the file name of the matrix to be used.
# The working directory used is either specified as a parameter on its
# own named 'wdir', or inferred from the basename of the matrix file,
# appended with some information on the splitting, e.g.  c160-4x4.


if (!defined($param->{'wdir'}) && defined($param->{'matrix'})) {
    # ok, we're doing this a bit ahead of time
    $wdir=$param->{'wdir'}="$param->{'matrix'}-${nh}x${nv}";
    push @main_args, "wdir=$wdir";
}

if ($main !~ /^:(?:complete|bench)$/ && !-d $param->{'wdir'}) {
    die "$param->{'wdir'}: no such directory";
}
# }}}

# Some quick shortcuts in case the user failed to specify some arguments
#  - splits=
#    --> splitting of vector blocks into sub-blocks to be treated as
#    independent sequences. This parameter is relevant globally, as it
#    has an impact on the bookkeeping around the prep and lingen steps.
#    For most common uses, we have splits of width 64 in the binary case,
#    or 1 in the GF(p) case. An SSE-2 kylov would need splits of width
#    128
#  - ys=
#    --> this is valid only for krylov or mksol, for *one* sequence run.
#    In general, this must correspond to one of the intervals defined in
#    "splits". For interleaving, this is a bit trickier, since two
#    contiguous intervals are to be treated.
# What the code below is just setting sensible defaults.
# {{{
##################################################
# Some more default values to set. Setting separately splits=, ys= and so
# on is quite awkward in the case of a single-site test.
if ((!defined($m) || !defined($n)) && $main !~ /^:wipeout$/) {
    usage "The parameters m and n must be set";
}
if (!defined($param->{'splits'})) {
    @splits=(0);
    my $x=0;
    while ($x<$n) {
        $x+=$splitwidth;
        push @splits, $x;
    }
    # This is probably ill-advised.
    # if ($param->{'interleaving'}) { @splits=(0,int($n/2),$n); }
    push @main_args, "splits=" . join(",", @splits);
}

if ($param->{'interleaving'}) {
    my @t = @splits;
    my $b = shift @t;
    while (scalar @t) {
        my $b0 = shift @t;
        my $b1 = shift @t;
        die "Interleaving can't work with splits ",
            join(",",@splits),
            "(intervals $b..$b0 and $b0..$b1 have different widths)"
            unless $b0-$b == $b1-$b0;
        $b = $b1;
    }
}

# After all we're probably better off decicing on what we do directly
# from with task_krylov.
#
# if (!defined($param->{'ys'})) { push @main_args, "ys=0..$n"; }
# }}}

# Done playing around with what we've been given.

##################################################
# System interaction.
# - set environment variables with my_setenv (which does a bit more)
# - run commands with dosystem
# - ssh to nodes with ssh_program ; this is MPI-backend-dependent.
# {{{
my $env_strings = "";
sub my_setenv
{
    return if (exists($ENV{$_[0]}) && $ENV{$_[0]} eq $_[1]);
    $env_strings .= "$_[0]=$_[1] ";
    $ENV{$_[0]}=$_[1];
}

sub dosystem
{
    my $nofatal=0;
    my $prg = shift @_;
    if ($prg eq '-nofatal') {
        $nofatal=1;
        $prg=shift @_;
    }
    my @args = @_;
    print STDERR '#' x 77, "\n";
    my $msg = "$env_strings$prg " . (join(' ', @args)) . "\n";

    if ($show_only) {
        print $msg;
        return 0;
    }

    print STDERR $msg;
    my $rc = system $prg, @args;
    return if $rc == 0;
    if ($rc == -1) {
        print STDERR "Cannot execute $prg\n";
    } elsif ($rc & 127) {
        my $sig = $rc & 127;
        my $coreinfo = ($rc & 128) ? 'with' : 'without';
        print STDERR "$prg: died with signal $sig, $coreinfo coredump\n";
    } else {
        my $ret = $rc >> 8;
        print STDERR "$prg: exited with status $ret\n";
    }
    exit 1;
}

sub ssh_program
{
    my $ssh;
    if ($ssh = $ENV{'SSH'}) {
        return $ssh;
    }
    $ssh='ssh';
    if (exists($ENV{'OAR_JOBID'})) {
        $ssh="/usr/bin/oarsh";
    }
    return $ssh
}
# }}}

if (!-d $wdir) {        # create $wdir on script node.
    if ($show_only) {
        print "mkdir $wdir\n";
    } else {
        mkdir $wdir;
    }
}

##################################################
# MPI-specific detection code. Beyond probing the environment variables
# to see which MPI middleware we're running, this code sets the following
# important argument lists:
#  - @mpi_precmd @mpi_precmd_single
#    --> Simply put, the former is prepended before all important
#    mpi-level programs (dispatch secure krylov mksol gather), while the
#    other is of course for leader-node-only programs.
#  - @mpi_precmd_nopthreads
#    --> exactly the same, but here we leave even thread management to
#    mpi. This applies to plingen_pz.
#  - $mpi_needed
#    --> to be used to check whether we want mpi. This is important
#    before calling dosystem.
# {{{

# If we've been built with mpi, then we _need_ mpi for running. Otherwise
# we run into shared libraries mess.
my $mpi_needed = $mpiexec ne '';

# @mpi_precmd_single is something silly; we want provision for the case
# where mpi is used for running non-mpi jobs. It's something which does
# not officially work, yet it always does. And some programs do turn out
# to be compiled with mpi, so we need the mpi libraries at runtime... So
# short of a more accurate solution, this is a hack.
my @mpi_precmd;
my @mpi_precmd_single;
my @mpi_precmd_nopthreads;


sub detect_mpi {
    if (defined($_=$ENV{'MPI_BINDIR'}) && -x "$_/mpiexec") {
        $mpi=$_;
    } elsif (defined($_=$ENV{'MPI'}) && -x "$_/mpiexec") { 
        $mpi=$_;
    } elsif (defined($_=$ENV{'MPI'}) && -x "$_/bin/mpiexec") { 
        $mpi="$_/bin";
    } elsif ($mpiexec && -x $mpiexec) {
        my($basename, $dirname) = fileparse($mpiexec);
        $mpi=$dirname;
    } else {
        my @path = split(':',$ENV{'PATH'});
        for my $d (@path) {
            if (-x "$d/mpiexec") {
                $mpi=$d;
                print STDERR "Auto-detected MPI_BINDIR=$mpi\n";
                last;
            }
        }
    }

    my $maybe_mvapich2=1;
    my $maybe_mpich2=1;
    my $maybe_openmpi=1;

    if (defined($mpi)) {
        SEVERAL_CHECKS: {
            # first check the alternatives system, which is fairly
            # commonplace.
            my $mpiexec = "$mpi/mpiexec";
            while (-l $mpiexec) {
                my $target=readlink($mpiexec);
                print STDERR "readlink($mpiexec)->$target\n";
                if ($target =~ m{^/}) {
                    $mpiexec = $target;
                } else {
                    $mpiexec =~ s{[^/]+$}{$target};
                }
                if ($mpiexec =~ /openmpi/) {
                    print STDERR "Auto-detecting openmpi based on alternatives\n";
                    $maybe_mvapich2=0;
                    $maybe_mpich2=0;
                    last;
                } elsif ($mpiexec =~ /mpich2/) {
                    print STDERR "Auto-detecting mpich2(old) based on alternatives\n";
                    $maybe_mvapich2=0;
                    $maybe_openmpi=0;
                    last;
                } elsif ($mpiexec =~ /hydra/) {
                    # Newer mvapich2 uses hydra as well...
                    print STDERR "Auto-detecting mpich or mvapich2 (hydra) based on alternatives\n";
                    $maybe_mvapich2='hydra';
                    $maybe_mpich2='hydra';
                    $maybe_openmpi=0;
                } elsif ($mpiexec =~ /mvapich2/) {
                    print STDERR "Auto-detecting mvapich2 based on alternatives\n";
                    $maybe_mpich2=0;
                    $maybe_openmpi=0;
                    last;
                }
            }
            CHECK_MVAPICH2: {
                if ($maybe_mvapich2 && -x "$mpi/mpiname") {
                    my $v = `$mpi/mpiname -n -v`;
                    chomp($v);
                    if ($v =~ /MVAPICH2\s+([\d\.]+)((?:\D\w*)?)/) {
                        # Presently all versions of mvapich2 up
                        # until 1.6rc3 included need mpd daemons.
                        # Released version 1.6 uses hydra.
                        $mpi_ver="mvapich2-$1$2";
                        if (($1 < 1.6) || ($1 == 1.6 && $2 =~ /^rc\d/)) {
                            $needs_mpd = 1;
                        } else {
                            $mpi_ver .= "+hydra" unless $mpi_ver =~ /hydra/;
                        }
                        last SEVERAL_CHECKS;
                    }
                }
            }
            CHECK_MPICH2_VERSION: {
                if ($maybe_mpich2 && -x "$mpi/mpich2version") {
                    my $v = `$mpi/mpich2version -v`;
                    chomp($v);
                    if ($v =~ /MPICH2 Version:\s*(\d.*)$/) {
                        $mpi_ver="mpich2-$1";
                        # Versions above 1.3 use hydra.
                        $needs_mpd=($mpi_ver =~ /^mpich2-(0|1\.[012])/);
                    } else {
                        $mpi_ver="mpich2-UNKNOWN";
                        $needs_mpd=1;
                    }
                    $v = `$mpi/mpich2version -c`;
                    chomp($v);
                    if ($v =~ /--with-pm=hydra/) {
                        $mpi_ver .= "+hydra";
                        $needs_mpd=0;
                    }
                    if ($maybe_mpich2 eq 'hydra') {
                        $mpi_ver .= "+hydra" unless $mpi_ver =~ /hydra/;
                        $needs_mpd=0;
                    }
                    last SEVERAL_CHECKS;
                }
            }
            CHECK_OMPI_VERSION: {
                if ($maybe_openmpi && -x "$mpi/ompi_info") {
                    my @v = `$mpi/ompi_info`;
                    my @vv = grep { /Open MPI:/; } @v;
                    last CHECK_OMPI_VERSION unless scalar @vv == 1;
                    $needs_mpd=0;
                    if ($vv[0] =~ /Open MPI:\s*(\d\S*)$/) {
                        $mpi_ver="openmpi-$1";
                        last SEVERAL_CHECKS;
                    } else {
                        $mpi_ver="openmpi-UNKNOWN";
                        last SEVERAL_CHECKS;
                    }
                }
            }
        }

        if (defined($mpi_ver)) {
            print STDERR "Using $mpi_ver, MPI_BINDIR=$mpi\n";
        } else {
            print STDERR "Using UNKNOWN mpi, MPI_BINDIR=$mpi\n";
            if (defined($needs_mpd=getenv("needs_mpd"))) {
                warn "Assuming needs_mpd=$needs_mpd as per env variable.\n";
            } else {
                $needs_mpd=1;
                warn "Assuming needs_mpd=$needs_mpd ; " .
                    "modify env variable \$needs_mpd to " .
                    "change fallback behaviour\n";
            }
        }
    }

    if (!defined($mpi)) {
        print STDERR <<EOMSG;
***ERROR*** No mpi library was detected. Arrange for mpiexec to be in
***ERROR*** your PATH, or set the MPI_BINDIR environment variable.
EOMSG
        exit 1;
    }
}
# 
# Starting daemons for mpich2 1.[012].x and mvapich2 ; we're assuming
# this works the same for older mpich2's, although this has never been
# checked.
sub check_mpd_daemons
{
    return unless $needs_mpd;

    my $ssh = ssh_program();

    my $rc = system "$mpi/mpdtrace > /dev/null 2>&1";
    if ($rc == 0) {
        print "mpi daemons seem to be ok\n";
        return;
    }

    if ($rc == -1) {
        die "Cannot execute $mpi/mpdtrace";
    } elsif ($rc & 127) {
        die "$mpi/mpdtrace died with signal ", ($rc&127);
    } else {
        print "No mpi daemons found by $mpi/mpdtrace, restarting\n";
    }

    open F, $hostfile;
    my %hosts;
    while (<F>) {
        $hosts{$_}=1;
    }
    close F;
    my $n = scalar keys %hosts;
    print "Running $n mpi daemons\n";
    dosystem "$mpi/mpdboot -n $n -r $ssh -f $hostfile -v";
}

sub get_mpi_hosts_torque {
    my @x = split /^/, eval {
		local $/=undef;
		open F, "$ENV{PBS_NODEFILE}";
		<F> };
    my $nth = $thr_split[0] * $thr_split[1];
    @hosts=();
    while (scalar @x) {
        my @z = splice @x, 0, $nth;
        my %h=();
        $h{$_}=1 for @z;
        die "\$PBS_NODEFILE not consistent mod $nth\n" unless scalar
			keys %h == 1;
        my $c = $z[0];
        chomp($c);
        push @hosts, $c;
    }
}

sub get_mpi_hosts_sge {
    print STDERR "Building hosts file from $ENV{PE_HOSTFILE}\n";
    my @x = split /^/, eval { local $/=undef; open F, "$ENV{PE_HOSTFILE}"; <F> };
    my $cores_on_node = {};
    for my $line (@x) {
        my ($node, $ncores, $toto, $tata) = split ' ', $line;
        print STDERR "$node: +$ncores cores\n";
        $cores_on_node->{$node}+=$ncores;
    }
    my $values_for_cores_on_node = {};
    local $_;
    $values_for_cores_on_node->{$_}=1 for (values %$cores_on_node);

    die "Not always the same number of cores obtained on the different nodes, as per \$PE_HOSTFILE" if keys %$values_for_cores_on_node != 1;

    my $ncores_obtained = (keys %$values_for_cores_on_node)[0];
    my $nnodes = scalar keys %$cores_on_node;
    print STDERR "Obtained $ncores_obtained cores on $nnodes nodes\n";

    my $nthr = $thr_split[0] * $thr_split[1];
    my $nmpi = $mpi_split[0] * $mpi_split[1];

    die "Not enough cores ($ncores_obtained) obtained: want $nthr\n" if $nthr > $ncores_obtained;

    @hosts=();
    push @hosts, $_ for keys %$cores_on_node;

    die "Not enough mpi nodes ($nnodes): want $nmpi\n" if $nmpi > $nnodes;
}

if ($mpi_needed) {
    # This is useful for debugging in case we see new MPI environments.
    print STDERR "Inherited environment:\n";
    print STDERR "$_=$ENV{$_}\n" for keys %$ENV;
    detect_mpi;

    push @mpi_precmd, "$mpi/mpiexec";

    # Need hosts.
    if (exists($ENV{'OAR_JOBID'}) && !defined($hostfile) && !scalar @hosts) {
	    print STDERR "OAR environment detected, setting hostfile.\n";
	    system "uniq $ENV{'OAR_NODEFILE'} > /tmp/HOSTS.$ENV{'OAR_JOBID'}";
	    $hostfile = "/tmp/HOSTS.$ENV{'OAR_JOBID'}";
    } elsif (exists($ENV{'PBS_JOBID'}) && !defined($hostfile) && !scalar @hosts ) {
	    print STDERR "Torque/OpenPBS environment detected, setting hostfile.\n";
	    get_mpi_hosts_torque;
    } elsif (exists($ENV{'PE_HOSTFILE'}) && exists($ENV{'NSLOTS'})) {
            print STDERR "Oracle/SGE environment detected, setting hostfile.\n";
	    get_mpi_hosts_sge;
    }

    if (scalar @hosts) {
        # Don't use an uppercase filename, it would be deleted by
        # wipeout.
        $hostfile = "$wdir/hosts";
        open F, ">$hostfile" or die "$hostfile: $!";
        for my $h (@hosts) { print F "$h\n"; }
        close F;
        print STDERR "Created $hostfile\n";
    }
    if (defined($hostfile)) {
        if ($needs_mpd) {
            # Assume daemons do the job.
        } elsif ($mpi_ver =~ /^\+hydra/) {
            my_setenv 'HYDRA_HOST_FILE', $hostfile;
            # I used to have various setups using --hostfile for openmpi,
            # --file in some other cases and so on. I think that
            # -machinefile is documented in the published standard, so
            # it's better to stick to it.
#        } elsif ($mpi_ver =~ /^openmpi/) {
#            push @mpi_precmd, "--hostfile", $hostfile;
#        } else {
#            push @mpi_precmd, "-file", $hostfile;
        } else {
            push @mpi_precmd, "-machinefile", $hostfile;
        }

    }
    if (!defined($hostfile)) {
        # At this point we're going to run processes on localhost.
        if ($mpi_ver =~ /^\+hydra/) {
            my_setenv 'HYDRA_USE_LOCALHOST', 1;
        }
        # Otherwise we'll assume that the simple setup will work fine.
    }
    check_mpd_daemons();
    if (!$needs_mpd) {
        # Then we must configure ssh
        # Note that openmpi changes the name of the option pretty
        # frequently.
        if ($mpi_ver =~ /^openmpi-1\.2/) {
            push @mpi_precmd, qw/--mca pls_rsh_agent/, ssh_program();
        } elsif ($mpi_ver =~ /^openmpi-1\.[34]/) {
            push @mpi_precmd, qw/--mca plm_rsh_agent/, ssh_program();
        } elsif ($mpi_ver =~ /^openmpi-1\.[56]/) {
            push @mpi_precmd, qw/--mca orte_rsh_agent/, ssh_program();
        } elsif ($mpi_ver =~ /^openmpi]/) {
            push @mpi_precmd, qw/--mca plm_rsh_agent/, ssh_program();
        } elsif ($mpi_ver =~ /^mpich2/ || $mpi_ver =~ /^mvapich2/) {
            # Not older mpich2's, which need a daemon.
            push @mpi_precmd, qw/-launcher ssh -launcher-exec/, ssh_program();
        }
    }
    if (defined($mpi_extra_args)) {
        push @mpi_precmd, split(' ', $mpi_extra_args);
    }

    @mpi_precmd_single = @mpi_precmd;
    @mpi_precmd_nopthreads = @mpi_precmd;
    push @mpi_precmd, '-n', $mpi_split[0] * $mpi_split[1];
    push @mpi_precmd_nopthreads, '-n', $nh * $nv;
    push @mpi_precmd_single, '-n', 1;

#} elsif (defined(my $mpi_path=$ENV{'MPI'})) {
#    my $ldlp;
#    if (defined($ldlp=$ENV{'LD_LIBRARY_PATH'})) {
#        $ldlp .= ':';
#    } else {
#        $ldlp = '';
#    }
#    $ldlp .= "$mpi_path/lib";
#    my_setenv 'LD_LIBRARY_PATH', $ldlp;
}

# }}}

if ($mpi_needed) {
    # mkdir must not be marked fatal, because if the command terminates
    # without having ever tried to join in an mpi collective like
    # mpi_init(), there's potential for the mpirun command to complain.
    dosystem('-nofatal', @mpi_precmd, split(' ', "mkdir -p $wdir"));
}

##################################################
### ok -- now @main_args is something relatively useful.
# print "main_args:\n", join("\n", @main_args), "\n";

push @main_args, splice @extra_args;
push @main_args, "matrix=$matrix";

# Now we have one function for each of the following macroscopic steps of
# the program
#
# dispatch
# prep
# krylov
# lingen
# mksol
# gather

# {{{ Some pretty-printing
my $current_task;

my $terminal_colors = {
    BLACK	=> "\e[01;30m",
    RED 	=> "\e[01;31m",
    GREEN	=> "\e[01;32m",
    YELLOW	=> "\e[01;33m",
    BLUE	=> "\e[01;34m",
    VIOLET	=> "\e[01;35m",
    black	=> "\e[00;30m",
    red         => "\e[00;31m",
    green	=> "\e[00;32m",
    yellow	=> "\e[00;33m",
    blue	=> "\e[00;34m",
    violet	=> "\e[00;35m",
};
$terminal_colors = {} if $ENV{'TERM'} !~ /^(xterm|screen|linux)/;

sub task_begin_message {
    my $blue = $terminal_colors->{'BLUE'} || '';
    my $normal = $terminal_colors->{'black'} || '';
    print "## Entering task: ${blue}$current_task${normal}\n";
}

sub task_check_message {
    my $status = shift;
    my $normal = $terminal_colors->{'black'} || '';
    my $color = {
        'ok' => $terminal_colors->{'green'} || '',
        'missing' => $terminal_colors->{'YELLOW'} || '',
        'error' => $terminal_colors->{'RED'} || '',
    };
    my @lines;
    push @lines, split(/^/m, $_) for @_;
    my $head = shift @lines;
    chomp($head);
    for (@lines) {
        chomp($_);
        $_ = "\t" . $_;
    }
    unshift @lines, $head;
    print "## $color->{$status}Check for $current_task$normal: $_\n" for @lines;
}
# }}}

# Inspection code for what is currently in $wdir

# Some of the inspection code correspond to queries which are done
# several times, so it's better to cache the results as long as they are
# expected to remain valid.
# {{{ cache hash, and corresponding queries
my $cache = {
    # This corresponds to data which isn't supposed to expire with
    # immediate effect.
    # Here we'll find things which are inferred from reading the data
    # files in $wdir, but would be cumbersome to find otherwise.
};
sub expire_cache_entry {
    my $k = shift;
    delete $cache->{$k};
}
sub store_cached_entry {
    my ($k, $v) = @_;
    if (defined($cache->{$k}) && $cache->{$k} != $v) {
        die "Fatal error: conflict for cache entry $k: previously had $cache->{$k}, now detected $v";
    }
    $cache->{$k} = $v;
}


sub get_cached_balancing_header {
    my $key = 'balancing_header';
    if (defined(my $z = $cache->{$key})) { return @$z; }
    sysopen(my $fh, $balancing, O_RDONLY) or die "$balancing: $!";
    sysread($fh, my $bhdr, 16);
    my @x = unpack("LLLL", $bhdr);
    $cache->{$key} = \@x;
    close($fh);
    return @x;
}
sub get_cached_leadernode_filelist {
    my $key = 'leadernode_filelist';
    my $opt = shift;
    my @x;
    if (defined(my $z = $cache->{$key})) {
        @x = @$z;
    } else {
        my $foo = join(' ', @mpi_precmd_single, "find $wdir -type f -a -printf '%s %p\\n'");
        for my $line (`$foo`) {
            $line =~ s/^\s*//;
            chomp($line);
            $line =~ s/^(\d+)\s+//;
            my $s = $1;
            push @x, [basename($line), $s];
        }
    }
    if ($opt && $opt eq 'HASH') { # we could also use wantarray
        my $h = {};
        $h->{$_->[0]} = $_->[1] for @x;
        return $h;
    } else {
        return @x;
    }
}
# }}}

# {{{ obtain_bfile -> check for balancing file.
sub obtain_bfile {
    # We're checking on the leader node only. Because the other nodes
    # don't really mind if they don't see the balancing file.
    # sense.
    my $pat;
    my $x = $matrix;
    $x =~ s{^(?:.*/)?([^/]+)$}{$1};
    $x =~ s/\.(?:bin|txt)$//;
    my $pre_pat = "$x\.${nh}x${nv}";
    if ($param->{'shuffled_product'}) {
        $pat = qr/\.([0-9a-f]{7}[13579bdf])\.bin\s*$/;
    } else {
        $pat = qr/\.([0-9a-f]{7}[02468ace])\.bin\s*$/;
    }
    if ($balancing) {
        $balancing =~ /$pat/ or die "$balancing does not match pattern";
        $balancing_hash = $1;
        return;
    }
    my @bfiles;
    my $hash;
    for my $fileinfo (get_cached_leadernode_filelist) {
        my ($file, $size) = @$fileinfo;
        $file =~ /^(.*)$pat/ && $1 eq $pre_pat or next;
        $hash = $2;
        push @bfiles, $file;
    }
    return unless scalar @bfiles;
    $balancing=$wdir . "/" . $bfiles[0];
    $balancing_hash = $hash;
    return @bfiles;
}
# }}}

# {{{ List the V files in $wdir -- output is given as a hash y-range =>
# iterations found. We also return the file size, which has to be
# constant.
sub list_files_generic {
    my $files = {};
    my $filesize;
    my ($n, $pattern) = @_;
    for my $fileinfo (get_cached_leadernode_filelist) {
        my ($file, $size) = @$fileinfo;
        $file =~ /^$pattern$/ or next;
        my $n0 = scalar @-;
        die "Found $n0 < $n matches. Bad pattern $pattern ??\n" if $n0 < $n;
        # See man perlvar
        my @matches = map { 0+substr $file, $-[$_], $+[$_] - $-[$_] } (1..$n0-1);
        my @kmatches = splice(@matches, 0, $n);
        push @{$files->{join("..", @kmatches)}}, \@matches;
        $filesize = $size if !defined($filesize);
        if ($filesize != $size) {
            task_check_message 'error', "Inconsistency detected for the sizes of the ${pattern} files. We have seen at least $filesize and $size (last seen: $_). Please fix.\n";
            die;
        }
    }
    return $files, $filesize;
}
sub list_vfiles {
    my ($f, $filesize) = list_files_generic(2, qr/V(\d+)-(\d+)\.(\d+)/);
    if ($filesize) {    # take the occasion to store it.
        my ($bnh, $bnv, $bnrows, $bncols) = get_cached_balancing_header;
        my $N = $bncols > $bnrows ? $bncols : $bnrows;
        store_cached_entry('nbytes_per_splitwidth', $filesize / $N);
    }
    my $flatten = sub { local $_; map { shift @$_ } @_; };
    $f->{$_} = [&$flatten(@{$f->{$_}})] for keys %$f;
    @{$f->{$_}} = sort { $a <=> $b } @{$f->{$_}} for keys %$f;
    return $f;
}
sub lexcmp {
    my $a = shift;
    my $b = shift;
    for my $k (0..$#$a) {
        my $z = $a->[$k] <=> $b->[$k];
        return $z if $z;
    }
    return 0;
}

sub list_afiles {
    my ($f, $filesize) = list_files_generic(2, qr/A(\d+)-(\d+)\.(\d+)-(\d+)/);
    if ($filesize) {    # take the occasion to store it.
        (my $k = (keys %$f)[0]) =~ /^(\d+)\.\.(\d+)$/;
        my $length = $m * ($2-$1) * ($f->{$k}->[0]->[1] - $f->{$k}->[0]->[0]);
        store_cached_entry('nbytes_per_splitwidth', $filesize / ($length / $splitwidth));
    }
    @{$f->{$_}} = sort { lexcmp($a, $b) } @{$f->{$_}} for keys %$f;
    return $f;
}

sub list_sfiles {
    my ($f, $filesize) = list_files_generic(4, qr/S\.sols(\d+)-(\d+).(\d+)-(\d+)\.(\d+)/);
    if ($filesize) {    # take the occasion to store it.
        my ($bnh, $bnv, $bnrows, $bncols) = get_cached_balancing_header;
        my $N = $bncols > $bnrows ? $bncols : $bnrows;
        store_cached_entry('nbytes_per_splitwidth', $filesize / $N);
    }
    my $flatten = sub { local $_; map { shift @$_ } @_; };
    $f->{$_} = [&$flatten(@{$f->{$_}})] for keys %$f;
    @{$f->{$_}} = sort { $a <=> $b } @{$f->{$_}} for keys %$f;
    return $f;
}
# }}}

# {{{ how many bytes per ($splitwidth) element ?
sub get_cached_nbytes_per_splitwidth {
    my $key = 'nbytes_per_splitwidth';
    if (defined(my $z = $cache->{$key})) { return $z; }
    # This function should practically never go beyond this point, as
    # presumably the list_* functions below, which all fill the cache as
    # a side effect, should have been run beforehand.
    list_vfiles;
    list_afiles;
    list_sfiles;
    if (defined(my $z = $cache->{$key})) { return $z; }
    die "Cannot find the number of bytes needed per finite field element";
}
# }}}


# {{{ inferring max iteration indices.
sub max_krylov_iteration {
    # read matrix dimension from the balancing header.
    my ($bnh, $bnv, $bnrows, $bncols) = get_cached_balancing_header;
    my $length = $bncols > $bnrows ? $bncols : $bnrows;
    $length = int(($length+$m-1)/$m) + int(($length+$n-1)/$n);
    $length += 2 * int(($m+$n-1)/$m);
    $length += 2 * int(($m+$n-1)/$n);
    $length += 10;
    return $length;
}
# This one is of course much harder to guess.
sub max_mksol_iteration {
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    my $nbytes_per_splitwidth = get_cached_nbytes_per_splitwidth;
    my @x;
    for my $file (keys %$leader_files) {
        $file =~ /^F\.sols(\d+)-(\d+)\.(\d+)-(\d+)$/ or next;
        my $size = $leader_files->{$file};
        return $size / (($2-$1)*($4-$3)/$splitwidth*$nbytes_per_splitwidth);
    }
    die "can't find any F file, cannot infer the mksol max iteration";
} # }}}

# {{{ task_common_run is just a handy proxy.
sub task_common_run {
    my $program = shift @_;
    expire_cache_entry 'leadernode_filelist';
    # Some arguments are relevant only to some contexts.
    @_ = grep !/^splits?=/, @_ unless $program =~ /split$/;
    @_ = grep !/^ys=/, @_ unless $program =~ /(?:krylov|mksol|dispatch)$/;
    @_ = grep !/^n?rhs=/, @_ unless $program =~ /(?:prep|gather|plingen.*|mksol)$/;
    @_ = grep !/shuffled_product/, @_ unless $program =~ /mf_bal/;

    $program="$bindir/$program";
    unshift @_, $program;

    if ($mpi_needed) {
        if ($program =~ /\/plingen[^\/]*$/) {
            unshift @_, @mpi_precmd_nopthreads;
        } elsif ($program =~ /\/(?:split|acollect|lingen|mf_bal|cleanup)$/) {
            unshift @_, @mpi_precmd_single;
        } elsif ($program =~ /\/(?:prep|secure|krylov|mksol|gather|dispatch)$/) {
            unshift @_, @mpi_precmd;
        } else {
            die "Don't know the parallel status of program $program ... ?";
        }
    }
    dosystem @_;
}
# }}}

# {{{ dispatch
sub task_dispatch {
    task_begin_message;
    # This does the matrix dispatch, or ensures that it's been already
    # done.
 
    # input files:
    #   well, nothing. The matrix, of course, but it is outside $wdir.
    #
    # output files, created if needed.
    #   - $balancing
    #   - cache files
    #
    # side-effect files, which have no impact on whether this step gets
    # re-run or not:
    #   - H1 Hx.0 ; these are sanity check vectors, used within dispatch (only)
    #   - cachelist.$balancing_hash.txt ; recomputed each time, just to
    #   see which cache files are where.
 
    my @bfiles = obtain_bfile;
    if (scalar @bfiles == 0) {
        print STDERR "No bfile found, we need a new one\n";
        my @mfbal=("mfile=$matrix", "out=$wdir/", $nh, $nv);
        unshift @mfbal, "--shuffled-product" if $param->{'shuffled_product'};
        unshift @mfbal, "--withcoeffs" if $prime ne '2';
        task_common_run "mf_bal", @mfbal;
        @bfiles = obtain_bfile; 
    }
    if (scalar @bfiles != 1) {
        print STDERR "Expected 1 bfile, found " . scalar @bfiles . ":\n";
        print "$_\n" foreach (@bfiles);
        die;
    }
    print "# Using balancing file $balancing\n";
    # XXX (TODO another way: inject $balancing within @main_args)
    # XXX Shall we accept @main_args as being possibly tinkered with by
    # the task_* functions ? This looks damn ugly...
    push @main_args, "balancing=$balancing";

    my @more_args = 'sequential_cache_build=1';

    # Note the ys=0..$splitwidth below. It's here only to match the width
    # which is used in the main krylov/mksol programs. The reason is that
    # the data width for vector coefficients gets passed to the cache
    # building code, which may take decisions depending on processor
    # cache size, which translate into a given number of vector
    # coefficients.
    # XXX For an SSE-2 binary, we would need something else.
    push @more_args, "ys=0..$splitwidth";

    # Do this always. This is a trivial operation, and we use it to see
    # who has what as far as cache files are concerned.
    task_common_run 'dispatch', (grep { !/^ys=/ } @main_args), @more_args, "export_cachelist=$wdir/cachelist.$balancing_hash.txt";

    # TODO: make sure that the matrix cache files are present everywhere.
    my $cachelist={};
    for my $line (eval { open F, "$wdir/cachelist.$balancing_hash.txt"; my @x=<F>; @x; }) {
        my ($token, $node, @files) = split(' ', $line);
        if ($token eq 'get-cache') {
            $cachelist->{$node}->{$_}=0 for @files;
        } elsif ($token eq 'has-cache') {
            $cachelist->{$node}->{$_}='found' for @files;
        }
    }
    my $missing=0;
    for my $node (keys %{$cachelist}) {
        my @bad = ();
        for my $file (keys %{$cachelist->{$node}}) {
            push @bad, $file unless $cachelist->{$node}->{$file};
        }
        if (@bad) {
            $missing += scalar @bad;
            print STDERR "## missing cache files for $node (" . (scalar @bad) . "): @bad\n";
        }
    }
    # Now we require that $missing is either 0 of $nh*$nv. Anthing
    # inbetween deserves inspection.
    if ($missing > 0 && $missing < $nh * $nv) {
        task_check_message 'error', "$missing missing cache files. We don't fix this automatically. Please investigate.";
        die;
    }

    if ($missing == $nh * $nv) {
        task_check_message 'missing', "no cache files have been found, need to run dispatch";
    } elsif ($missing == 0) {
        task_check_message 'ok', "all cache files have been found, good";
    }

    # Note that the sanity check vectors Hx.0 and H1 are just
    # side-effects of the dispatch program, but there is no need to keep
    # them around. They're never reused.

    return unless $missing;

    push @more_args, "sanity_check_vector=H1" if $prime == 2;


    task_common_run('dispatch', @main_args, @more_args);

} # }}}

# {{{ prep
sub task_prep {
    task_begin_message;
    # This prepares the starting vectors for block Wiedemann.
 
    # input files:
    #   - output files from previous tasks: dispatch
    #
    # output files, created if needed.
    #   - X
    #   - V${x0}-${x1}.0 ; the different starting vectors.
    #
    # side-effect files, which have no impact on whether this step gets
    # re-run or not:
    #   - Y.0 ; this is a concatenation of the starting vectors in the
    #   binary case. This must now be viewed as an intermediary file
    #   only, and will disappear at some point.

    my @s = @splits;
    my @starts;
    while (scalar @s) {
        my $x0 = shift @s;
        my $x1 = shift @s or last;
        unshift @s, $x1;
        push @starts, "V${x0}-${x1}.0";
    }
    my @missing;
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    for my $file ("X", @starts) {
        push @missing, $file unless exists $leader_files->{$file};
    }
    if (!@missing) {
        task_check_message 'ok', "all starting files are here, good";
        return;
    }
    task_check_message 'missing', "missing starting files: @missing";

    if ($prime == 2) {
        task_common_run('prep', @main_args);
        task_common_run('split',
            (grep { /^(?:m|n|wdir|prime|splits)=/ } @main_args),
            qw{--ifile Y.0 --ofile-fmt V%u-%u.0});
    } else {
        # The prime case is somewhat different. We'll generate random
        # data, that's it. We might prepend the RHS if it exists.
        task_common_run('prep', @main_args, "ys=0..$splitwidth");
    }

} # }}}

# {{{ secure -- this is now just a subtask of krylov or mksol
sub subtask_secure {
    return if $param->{'skip_online_checks'};
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    my @x = grep { !exists($leader_files->{$_}); } "C.$param->{'interval'}";
    if (@x) {
        task_check_message 'missing', "missing check vector @x\n";
        task_common_run('secure', @main_args);
    }
}
# }}}

# {{{ subtask_krylov_mksol_todo - the hard checkpoint recovery work
sub subtask_krylov_mksol_todo {
    # For each sequence, this finds the "most advanced" V file. Being
    # advanced takes two things. First this means an existing file, and
    # the code here checks for this. But we also request that some extra
    # function returns true, and this check is provided by the caller.
    my $length = shift;
    my $morecheck = shift || sub{1;};
    # We unconditionally do a check of the latest V checkpoint found in
    # $wdir, for all vector.
    my @all_ys = map { [ $_*$splitwidth, ($_+1)*$splitwidth ] } (0..$n/$splitwidth - 1);

    my $vfiles = list_vfiles;
    my $vstarts = {};

    for my $y (@all_ys) {
        my $yrange = $y->[0] . ".." . $y->[1];
        print "## V files for $yrange:";
        my @for_this_y =
            grep { $_ == 0 || &$morecheck($yrange, $_); } @{$vfiles->{$yrange}}
            or do { print " none\n"; next; };
        $vstarts->{$yrange} = $for_this_y[$#for_this_y];
        $vstarts->{$yrange} .= " (DONE!)" if $vstarts->{$yrange} >= $length;
        print " last $current_task checkpoint is $vstarts->{$yrange}\n";
    }

    my @ys;
    if (@ys = grep { /^ys/ } @main_args) {
        @ys = map { /^(?:ys=)?(\d+)\.\.(\d+)$/; $_=[$1, $2]; } @ys;
    } elsif ($param->{'interleaving'}) {
        my @t = @all_ys;
        # merge 2 by 2.
        @ys = ();
        while (scalar @t >= 2) {
            my ($a, $b0) = @{shift @t};
            my ($b1, $c) = @{shift @t};
            push @ys, [ $a, $c ];
        }
        # Maybe there's a last one.
        push @ys, @t;
    } else {
        @ys = @all_ys;
    }

    my @todo;
    my @impossible;
    for my $ab (@ys) {
        # most ys are double ranges, except maybe the last one.
        my ($a, $b) = @{$ab};
        my $yrange = $a . ".." . $b;
        my $start;
        if ($param->{'interleaving'} && !($b - $a == $splitwidth && $b == $n)) {
            # interleaving is quite special. At the moment we must make
            # sure that hte start points for both vectors agree, although
            # admittedly this restriction is quite artificial.
            next if $b - $a == $splitwidth && $b == $n;
            die unless $b-$a == 2 * $splitwidth;
            my @subs = ([$a, int(($a+$b)/2)], [int(($a+$b)/2), $b]);
            my @sub_starts;
            for (@subs) {
                my $r = $_->[0]."..".$_->[1];
                my $s = $vstarts->{$r} or do {
                    task_check_message 'error',
                    "No starting vector for range $r"
                    . " (interleaved sub-range from $yrange)\n";
                    die;
                };
                push @sub_starts, $s;
            } 
            if ($sub_starts[0] ne $sub_starts[1]) {
                task_check_message 'error',
                "Inconsistent start vectors for the two"
                . " interleaved sub-ranges from $yrange"
                . " (first at $sub_starts[0],"
                . " second at $sub_starts[1])\n";
                die;
            };
            $start = $sub_starts[0];
        } else {
            $start = $vstarts->{$yrange};
        }
        if (!defined($start)) {
            print STDERR "## Can't schedule any work for $yrange, as *NO* checkpoint is here\n";
            push @impossible, $yrange;
            next;
        }
        next if $start =~ /DONE/;
        push @todo, "ys=$yrange start=$start";
    }
    if (!@todo && @impossible) {
        task_check_message 'error', "Cannot schedule remaining work. No checkpoints for ranges:\n", @impossible;
        die;
    }
    return @todo;
}
# }}}

# {{{ krylov
sub task_krylov {
    task_begin_message;

    # input files:
    #   - output files from previous tasks: dispatch prep
    #     More accurately, we need V${x0}-${x1}.$i for some iteration i.
    #
    # if ys is found in the command line, then we focus on that sequence
    # specifically. Otherwise, we process all sequences one after
    # another.
    #
    # Other:
    #   - C.$interval ; check vector for online checks. Computation
    #   skipped if explicitly requested (skip_online_checks=1).
    #
    # side-effect files, which have no impact on whether this step gets
    # re-run or not:
    #   - A${x0}-${x1}.* ; dot product files which are fed to lingen.
    #   These are write-only as far as this task is concerned.

    subtask_secure;

    my $length = max_krylov_iteration;
    print "## krylov max iteration is $length\n";

    my @todo = subtask_krylov_mksol_todo $length;

    if (!@todo) {
        # Note that we haven't checked for the A files yet !
        task_check_message 'ok',
                    "All krylov tasks are done, good.\n";
        return;
    }

    task_check_message 'missing',
                "Pending tasks for krylov:\n" . join("\n", @todo);
    for my $t (@todo) {
        # take out ys from main_args, put the right one in place if
        # needed.
        my @args = grep { !/^ys/ } @main_args;
        push @args, split(' ', $t);
        task_common_run 'krylov', @args;
    }
} # }}}

# {{{ lingen
sub task_lingen {
    task_begin_message;
    
    # input files:
    #   - A.*-*.*-* ; all A files created by krylov.
    #     -> this task accumulates them in a single file.
    # output_file
    #   - F.sols*-*.*-* ; all polynomial entries of the generating
    #     matrix.

    my $length = max_krylov_iteration;
    my $afiles = list_afiles;
    my $astrings = {};
    my $ok = 1;
    my $nb_afiles = 0;
    my $rlength = int(($length + $param->{'interval'} - 1) / $param->{'interval'}) * $param->{'interval'};
    for my $k (keys %$afiles) {
        my @strings;
        my ($z0, $z1);
        for my $seg (@{$afiles->{$k}}) {
            $nb_afiles++;
            if (defined($z1)) {
                if ($seg->[0] == $z1) {
                    $z1 = $seg->[1];
                } else {
                    push @strings, [$z0, $z1];
                    $z0 = undef;
                    $z1 = undef;
                }
                next;
            }
            ($z0, $z1) = @$seg;
        }
        $ok = 0 unless scalar @strings == 0 && $z0 == 0 && $z1 >= $length;
        push @strings, [$z0, $z1] if defined $z1;
        $astrings->{$k} = \@strings;
    }
    if (!$ok) {
        my $msg='';
        for my $k (sort {$a cmp $b} keys %$afiles) {
            $msg .= "$k, computed range(s): " . join(" ", map { $_->[0]."..".$_->[1]; } @{$astrings->{$k}}) . "\n";
        }

        task_check_message 'error', "Incomplete set of A files found\n" . $msg;
        die;
    }
    my @all_astrings;
    task_check_message 'ok', "Complete set of A files found ($nb_afiles files), good.";

    if ($nb_afiles > 1) {
        task_check_message 'missing', "Running acollect to create A0-$n.0-$rlength";
        task_common_run("acollect", @main_args, "--remove-old");
    }

    # We expect lingen to split its output into pieces of predictable
    # size.
    my $leader_files = get_cached_leadernode_filelist 'HASH';
    my @missing;
    my @expected;
    for my $i (0..$n/$splitwidth-1) {
        my $i0 = $splitwidth * $i;
        my $i1 = $splitwidth + $i0;
        for my $j (0..$n/$splitwidth-1) {
            my $j0 = $splitwidth * $j;
            my $j1 = $splitwidth + $j0;
            my $f = "F.sols${i0}-${i1}.${j0}-${j1}";
            push @expected, $f;
            push @missing, $f unless exists $leader_files->{$f};
        }
    }
    if (@missing == 0) {
        task_check_message 'ok', "lingen result found, good.";
        return;
    } elsif (@missing < @expected) {
        task_check_message 'error', "Incomplete lingen output found. Missing files:\n" , @missing;
        die;
    }

    task_check_message 'missing', "lingen has not run yet. Running now.";
    # Now run lingen itself. Which binary we'll run is not totally
    # obvious though.
    if ($prime == 2) {
        task_common_run("lingen", @main_args, qw/--lingen-threshold 64/);
        task_common_run('split',
            (grep { /^(?:m|n|wdir|prime|splits)=/ } @main_args),
            split(' ', "--ifile F --ofile-fmt F.sols0-$n.%u-%u"));
    } else {
        # NOTE: It may be worthwhile to run specifically this step, but
        # with adapted mpi and thr parameters.
        my @args;
        my $lt = $param->{'lingen_threshold'} || 10;
        my $lmt = $param->{'lingen_mpi_threshold'} || 100;
        push @args, "lingen-threshold=$lt";
        push @args, "lingen-mpi-threshold=$lt";
        push @args, "afile=A0-$n.0-$rlength";
        push @args, grep { /^(?:m|n|wdir|prime|nrhs|mpi|thr)=/ } @main_args;
        task_common_run("plingen_pz", @args);
        # Some splitting work needed...
        @args=();
        push @args, "splits=" . join(",",@splits);
        push @args, grep { /^(?:m|n|wdir|prime)=/ } @main_args;
        task_common_run "split", @args,
                "ifile=A0-$n.0-$rlength.gen", "ofile-fmt=F.%u-%u";
        for my $j (0..$n/$splitwidth-1) {
            my $j0 = $splitwidth * $j; 
            my $j1 = $splitwidth + $j0;
            my $dseq = $j0 . "-" . $j1;
            task_common_run "split", @args, "ifile=F.$dseq",
            "ofile-fmt=F.sols%u-%u.$dseq";
        }
    }
}
# }}}

# {{{ mksol
sub task_mksol {
    task_begin_message;

    # input files:
    #   - output files from dispatch
    #   - V*-*.<some iteration>, e.g. maybe the starting one generated by
    #     prep, or any other.
    #   - F.sols*-*.*-* ; all polynomial entries of the generating
    #     matrix.
    #
    # Note that the computation will start at the first iteration number
    # for which a V file is here, and no S file yet (for any of the
    # solutions considered, which is a number limited by $nrhs if
    # specified).
    #
    # Other:
    #   - C.$interval ; check vector for online checks. Computation
    #   skipped if explicitly requested (skip_online_checks=1).
    #
    # side-effect files, which have no impact on whether this step gets
    # re-run or not:
    #   - S.sols*-*.*-*.* ; partial sum which are fed to gather.
    #   These are write-only as far as this task is concerned.

    subtask_secure;

    # We choose to require the S files for considering checkpoints as
    # valid. This is because in all likelihood, the V files from krylov
    # will still be there, so even though this is a hint as to what can
    # be started right now, we can't really use it...
    my $sfiles = list_sfiles;
    $sfiles->{$_} = eval { my %x=map {$_=>1;} @{$sfiles->{$_}};\%x;} for keys %$sfiles;

    my $length = max_mksol_iteration;

    my $n_or_rhs = $nrhs || $n;
    my @sols = map { $_*=$splitwidth; $_."..".($_+$splitwidth); } (0..$n_or_rhs/$splitwidth-1);
    print "## mksol considering the following solution blocks: @sols\n";
    print "## mksol max iteration is $length\n";

    my @todo = subtask_krylov_mksol_todo $length, sub {
        my ($range, $cp) = @_;
        ## return if $cp > $length;
        for (@sols) {
            return unless $sfiles->{$_ . ".." . $range}->{$cp};
        }
        return 1;
    };

    if (!@todo) {
        # Note that we haven't checked for the A files yet !
        task_check_message 'ok', "All $current_task tasks are done, good.\n";
        return;
    }

    task_check_message 'missing',
                "Pending tasks for $current_task:\n" . join("\n", @todo);
    for my $t (@todo) {
        # take out ys from main_args, put the right one in place if
        # needed.
        my @args = grep { !/^(ys|n?rhs)/ } @main_args;
        push @args, split(' ', $t);
        my $nsolvecs = $nrhs || $splitwidth;
        push @args, "nsolvecs=$nsolvecs";

        task_common_run 'mksol', @args;
    }
}
# }}}

sub task_gather {
    task_begin_message;

    # input files:
    #   - S.sols*-*.*-*.* ; partial sums computed by mksol

    my $n_or_rhs = $nrhs || $n;
    my @sols = map { $_*=$splitwidth; $_."..".($_+$splitwidth); } (0..$n_or_rhs/$splitwidth-1);
    print "## gather considering the following solution blocks: @sols\n";
    my $sfiles = list_sfiles;
    my @missing;
    my $sfile_count_matrix = {};
    for my $k (keys %$sfiles) {
        $k =~ /^(\d+)\.\.(\d+)\.\.(\d+)\.\.(\d+)$/ or die;
        $sfile_count_matrix->{($1/$splitwidth).",".($3/$splitwidth)} = scalar @{$sfiles->{$k}};
    }
    my $c = 4;
    {
        my $i = -1;
        for my $j (0..$n/$splitwidth-1) {
            $sfile_count_matrix->{"$i,$j"}=$splitwidth*$j."-".($splitwidth*($j+1));
            my $l = length($sfile_count_matrix->{"$i,$j"});
            $c = $l if $l > $c;
        }
    }
    for my $i (0..$n_or_rhs/$splitwidth-1) {
        $sfile_count_matrix->{"$i,-1"}=$splitwidth*$i."-".($splitwidth*($i+1));
        my $l = length($splitwidth*$i);
        $c = $l if $l > $c;
        for my $j (0..$n/$splitwidth-1) {
            my $n = $sfile_count_matrix->{"$i,$j"} || 'NONE';
            push @missing, "S.sols$i-".($i+$splitwidth).".$j".($j+$splitwidth)
                if $n eq 'NONE';
            my $l = length($sfile_count_matrix->{"$i,$j"} || 'NONE');
            $c = $l if $l > $c;
        }
    }
    my @count_matrix_rows;
    for my $i (-1..$n_or_rhs/$splitwidth-1) {
        my @row;
        for my $j (-1..$n/$splitwidth-1) {
            my $t = $sfile_count_matrix->{"$i,$j"};
            $t = '' if $i==-1 && $j==-1;
            $t = sprintf("%${c}s", $t);
            $t .= "|" if $j == -1;
            push @row, $t;
        }
        $row[0] = ($i>=0 ? "sols" : "    ") . $row[0];
        push @count_matrix_rows, "##    @row\n";
        if ($i == -1) {
            my $sep = "";
            $sep .= "-" x (length($row[0])-1);
            $sep .= "+";
            $sep .= "-" x (length("@row[1..$#row]")+1);
            push @count_matrix_rows, "##    $sep\n";
        }
    }
    print "## Number of S files found:\n" . join("", @count_matrix_rows);
    if (@missing) {
        task_check_message 'error', "Missing files for $current_task:", @missing;
        die;
    }

    my $leader_files = get_cached_leadernode_filelist 'HASH';

    my $rhs_companion;
    if ($param->{'rhs'}) {
        my $length = max_krylov_iteration;
        my $rlength = int(($length + $param->{'interval'} - 1) / $param->{'interval'}) * $param->{'interval'};
        $rhs_companion = "A0-$n.0-$rlength.gen.rhs";
        if (!exists($leader_files->{$rhs_companion})) {
            task_check_message 'error', "Missing file (for RHS): $rhs_companion\n";
            die;
        }
    }

    task_check_message 'ok', "All required files for gather seem to be present, good.";

    my @args = grep { !/^(ys|rhs)/ } @main_args;
    push @args, "nsolvecs=$splitwidth";
    if ($param->{'rhs'}) {
        push @args, "rhs=$rhs_companion";
    }
    task_common_run 'gather', @args;
}

my @tasks = (
    [ 'dispatch', \&task_dispatch],
    [ 'prep', \&task_prep],
    [ 'krylov', \&task_krylov],
    [ 'lingen', \&task_lingen],
    [ 'mksol', \&task_mksol],
    [ 'gather', \&task_gather],
);

for my $tc (@tasks) {
    if ($main eq $tc->[0] || $main eq ':complete') {
        $current_task = $tc->[0];
        &{$tc->[1]}(@main_args);
    }
}
