#!/usr/bin/perl
use warnings;
use strict;

use POSIX qw/getcwd/;
use File::Basename;
use File::Temp qw/tempdir/;
use List::Util qw[max];

# This companion program serves as a helper for running bwc programs.

# MPI_BINDIR=/opt/mpich2-1.1a2-1.fc10.x86_64/usr/bin ./bwc.pl :complete matrix=c72 interval=256 splits=64 mn=64 ys=0..64 mpi=4x4

sub usage {
    print STDERR <<EOF;
Usage: ./bwc.pl <action> <parameters>

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
mn m n               same meaning as for bwc programs.
nullspace            same meaning as for bwc programs.
interval             same meaning as for bwc programs.
ys splits            same meaning as for bwc programs (krylov and mksol)

matrix matpath       input matrix file. Either 'matrix' is a full
                     path, or a basename relative to 'matpath'
wdir                 working directory.
                     If unset, derived from 'matrix', 'mpi' 'thr'
mode                 Family of binaries to be used for e.g. krylov/mksol
mpiexec              Command to invoke for mpiexec
hosts                Comma-separated list of hosts. 'hosts' may be supplied
                     several times, so that lists concatenate. Thus one
                     may use the shell to supply hosts=node{1,3,5,7}.
hostfile             path to a host list file to be understood by mpiexec
                     NOTE: multi-hosts support is incomplete in this
                     script at the moment, and specific to mpich2+pm=hydra.

Amongst parameters, some options with leading dashes may be set:

-d  Do a dry run ; show the commands that would be executed.
    (-d is relevant to the script only)
-v  Increase verbosity level ; understood by bwc, but meaningless at the moment
-h  Show this help.

Some environment variables have an impact on this script.

HOME         The default value for matpath is \$HOME/Local/mats
MPI_BINDIR   The path to the directory holding the mpiexec program (the
             mpiexec= parameter allows to specify the same thing. Note that
             MPI alone serves as an alias for MPI_BINDIR, but this may
             become legacy).
BWC_BINDIR   The directory holding the bwc binaries. Defaults to the
             directory holding the current script.

EOF

    if (scalar @_) {
        print STDERR "Error messages:\n", join("\n", @_), "\n";
    }

    exit 1;
}

# $program denotes either a simple command, or something prepended by ':'
# which indicates something having a special meaning for the script.
my $main = shift @ARGV;

# First we're going to parse globally the list of arguments, to see if
# there's something to infer.

my $matpath="$ENV{HOME}/Local/mats";
my @main_args = ();
my @extra_args = ();
my @mpi_split=(1,1);
my @thr_split=(1,1);
my $matrix;
my $balancing;
my $balancing_hash;
my $hostfile;

## mpiexec is substituted by cmake in case mpi has been used for the
## compilation. NOTE that this means that a priori, mpiexec _must_ be
## used for running all programs.
my $mpiexec='@MPIEXEC@';

## The mpi_extra_args argument is used to pass information to the mpiexec 
## command. The idea is that mpiexec, the mpi driver program, may need
## additional info to properly setup communications between jobs. As an 
## example, on Fedora machines with openmpi:
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
 
my $wdir;
my $m;
my $n;
my @splits=();
my @hosts=();
my $mode='u64';
my $show_only=0;
my $nh;
my $nv;
my $mpi;
my $mpi_ver;
my $tmpdir;
my $interleaving;
my $force_complete;

print $0, " ", join(" ", @ARGV), "\n";

# {{{ MPI detection

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

    if (defined($mpi)) {
        CHECK_MPICH2_VERSION: {
            if (-x "$mpi/mpich2version") {
                my $v = `$mpi/mpich2version -v`;
                chomp($v);
                if ($v =~ /MPICH2 Version:\s*(\d.*)$/) {
                    $mpi_ver="mpich2-$1";
                } else {
                    $mpi_ver="mpich2-UNKNOWN";
                }
                $v = `$mpi/mpich2version -c`;
                chomp($v);
                if ($v =~ /--with-pm=hydra/) {
                    $mpi_ver .= "+hydra";
                }
            }
        }
        CHECK_OMPI_VERSION: {
            if (-x "$mpi/ompi_info") {
                my @v = `$mpi/ompi_info`;
                my @vv = grep { /Open MPI:/; } @v;
                last CHECK_OMPI_VERSION unless scalar @vv == 1;
                if ($vv[0] =~ /Open MPI:\s*(\d\S*)$/) {
                    $mpi_ver="openmpi-$1";
                } else {
                    $mpi_ver="openmpi-UNKNOWN";
                }
            }
        }

        if (defined($mpi_ver)) {
            print STDERR "Using $mpi_ver, MPI_BINDIR=$mpi\n";
        } else {
            print STDERR "Using UNKNOWN mpi, MPI_BINDIR=$mpi\n";
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
# }}}


my $bindir;
if (!defined($bindir=$ENV{'BWC_BINDIR'})) {
    $bindir=$0;
    $bindir =~ s{/[^/]*$}{};
    if ($bindir eq $0) {
        $bindir='.';
    }
}

my $param={};

while (defined($_ = shift @ARGV)) {
    # -d will never be found as a bw argument, it is 
    if ($_ eq '--') {
        push @extra_args, $_, splice @ARGV;
        last;
    }
    if (/^(-d|--show|--dry-run)$/) { $show_only=1; next; }
    if (/^(-h|--help)$/) { usage; }
    my ($k,$v);
    if (/^(-.*)$/) { $k=$1; $v=1; }
    if (/^([^=]+)=(.*)$/) { $k=$1; $v=$2; }
    if (!defined($k)) {
        usage "Garbage not undertood on command line: $_";
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

sub set_mpithr_param {
    my $v = shift @_;
    $v=~/(\d+)x(\d+)$/ or usage "bad splitting value '$v'";
    return ($1, $2);
}

# Some parameters are more important than the others because they
# participate to the default value for wdir.
if ($param->{'mpi'}) { @mpi_split = set_mpithr_param $param->{'mpi'}; }
if ($param->{'thr'}) { @thr_split = set_mpithr_param $param->{'thr'}; }
if ($param->{'wdir'}) { $wdir=$param->{'wdir'}; }
if ($param->{'matrix'}) { $matrix=$param->{'matrix'}; }
if ($param->{'matpath'}) { $matpath=$param->{'matpath'}; }
if ($param->{'tmpdir'}) { $tmpdir=$param->{'tmpdir'}; }
if ($param->{'interleaving'}) { $interleaving=$param->{'interleaving'}; }

$nh = $mpi_split[0] * $thr_split[0];
$nv = $mpi_split[1] * $thr_split[1];

# Selection of the working directory -- there are three relevant options.
# The 'matrix' option gives the file name of the matrix to be used. If
# 'matrix' gives only a basename, then this basename is interpreted
# relative to the 'matpath' parameter. The working directory used is
# either specified as a parameter on its own named 'wdir', or inferred
# from the basename of the matrix file, appended with some information on
# the splitting, e.g.  c160-4x4.

if (!defined($wdir) && defined($matrix)) {
    if ($matrix !~ m{/} && defined($matpath)) {
        $matrix = "$matpath/$matrix";
    }
    $matrix =~ s/~/$ENV{HOME}/g;
    $wdir="$matrix-${nh}x${nv}";
    # In most cases, a non-existing wdir is an error, but not for
    # :complete, which is entitled to create it.
    if ($main !~ /^:(?:complete|bench)$/ && !-d $wdir) {
        die "No directory $wdir found";
    }
    push @main_args, "wdir=$wdir";

    # Put back in our parameter list.
    $param->{'matrix'}=$matrix;
    $param->{'wdir'}=$wdir;
}

if (!defined($wdir) && !defined($matrix)) {
    die "At least one of matrix or wdir must be defined";
}

# Assuming that the parameters given on the command line are sufficient
# to determine an existing wdir, then this wdir is scanned for an
# existing bw.cfg file. This file is read. However, parameters read from
# the bw.cfg have _lower_ precedence than the ones on the command line.

if (defined($wdir) && -d $wdir && -f "$wdir/bw.cfg") {
    open F, "$wdir/bw.cfg";
    while (<F>) {
        chomp($_);
        next if /^#/;
        next if /^\s*$/;
        my ($k,$v);
        if (/^(-.*)$/) { $k=$1; $v=1; }
        if (/^([^=]+)=(.*)$/) { $k=$1; $v=$2; }
        if (!defined($k)) {
            usage "Garbage not undertood in $wdir/bw.cfg:\n$_";
        }
        if (!defined($param->{$k})) {
            $param->{$k}=$v;
        }
    }
    close F;
}

while (my ($k,$v) = each %$param) {
    # First the ones that are _not_ relevant to the bwc programs.
    if ($k eq 'matrix') { $matrix=$v; next; }
    if ($k eq 'matpath') { $matpath=$v; next; }
    if ($k eq 'hostfile') { $hostfile=$v; next; }
    if ($k eq 'mpi_extra_args') { $mpi_extra_args=$v; next; }
    if ($k eq 'force_complete') { $force_complete=$v; next; }
    if ($k eq 'mode') { $mode=$v; next; }
    if ($k eq 'hosts') {
        $v=[$v] if (ref $v eq '');
        for (@$v) { push @hosts, split(',',$_); }
        next;
    }
    # The rest is passed to subprograms, unless explicitly discarded on a
    # per-program basis.
    push @main_args, "$k=$v";

    # Yet it does make sense to read some parameters, are we are
    # interested by their value.

    if ($k eq 'wdir') { $wdir=$v; next; }
    if ($k eq 'mn') { $m=$n=$v; next; }
    if ($k eq 'n') { $n=$v; next; }
    if ($k eq 'm') { $m=$v; next; }
    if ($k eq 'splits') { @splits=split ',', $v; next; }
    if ($k eq 'mpi') { @mpi_split = set_mpithr_param $v; }
    if ($k eq 'thr') { @thr_split = set_mpithr_param $v; }

    # Ok, probably there's a fancy argument we don't care about.
}

# refresh. This may seem silly, but if wdir was non-default, perhaps
# we've had _no_ information on nh and nv so far !
$nh = $mpi_split[0] * $thr_split[0];
$nv = $mpi_split[1] * $thr_split[1];


##################################################
# Some more default values to set. Setting separately splits=, ys= and so
# on is quite awkward in the case of a single-site test.
if ((!defined($m) || !defined($n)) && $main !~ /^:wipeout$/) {
    usage "The parameters m and n must be set";
}
if (!defined($param->{'splits'})) {
    if ($param->{'interleaving'}) {
        @splits=(0,int($n/2),$n);
    } else {
        @splits=(0,$n);
    }
    push @main_args, "splits=" . join(",", @splits);
}
if (!defined($param->{'ys'})) {
    push @main_args, "ys=0..$n";
}


##################################################
# Done playing around with what we've been given.

##################################################
##################################################
# System interaction.

my $env_strings = "";
sub my_setenv
{
    return if (exists($ENV{$_[0]}) && $ENV{$_[0]} eq $_[1]);
    $env_strings .= "$_[0]=$_[1] ";
    $ENV{$_[0]}=$_[1];
}

sub dosystem
{
    my $prg = shift @_;
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

##################################################
# Do we need mpi ??

# Starting daemons for mpich2 1.0.x
sub check_mpd_daemons()
{
    return if !defined $mpi_ver;
    return if $mpi_ver !~ /^mpich2/;
    return if $mpi_ver =~ /^mpich2.*\+hydra/;

    my $rsh='ssh';
    if (exists($ENV{'OAR_JOBID'})) {
        $rsh="/usr/bin/oarsh";
    }

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
    dosystem "$mpi/mpdboot -n $n -r $rsh -f $hostfile -v";
}

# my $mpi_needed = $mpi_split[0] * $mpi_split[1] != 1;
# && $main !~ /(?:split|acollect|lingen)$/;

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

if (!-d $wdir) {
    if ($show_only) {
        print "mkdir $wdir\n";
    } else {
        mkdir $wdir;
    }
}

if ($mpi_needed) {
    detect_mpi;

    push @mpi_precmd, "$mpi/mpiexec";
    push @mpi_precmd_single, "$mpi/mpiexec";

    # Need hosts.
    if (exists($ENV{'OAR_JOBID'}) && !defined($hostfile) && !scalar @hosts) {
        print STDERR "OAR environment detected, setting hostfile.\n";
        system "uniq $ENV{'OAR_NODEFILE'} > /tmp/HOSTS.$ENV{'OAR_JOBID'}";
        $hostfile = "/tmp/HOSTS.$ENV{'OAR_JOBID'}";
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
        if ($mpi_ver =~ /^\+hydra/) {
            my_setenv 'HYDRA_HOST_FILE', $hostfile;
        } elsif ($mpi_ver =~ /^openmpi/) {
            push @mpi_precmd, "--hostfile", $hostfile;
        } elsif ($mpi_ver =~ /^mpich2/) {
            # Assume daemons do the job.
        } else {
            push @mpi_precmd, "-file", $hostfile;
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
    push @mpi_precmd, '-n', $mpi_split[0] * $mpi_split[1];
    push @mpi_precmd_single, '-n', 1;

    if (defined($mpi_extra_args)) {
        push @mpi_precmd, split(' ', $mpi_extra_args);
    }
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

##################################################
### ok -- now @main_args is something relatively useful.
# print "main_args:\n", join("\n", @main_args), "\n";

push @main_args, splice @extra_args;
push @main_args, "matrix=$matrix";

sub obtain_bfile {
    return if $balancing;
    opendir(my $dh, $wdir);
    my $x = $matrix;
    $x =~ s/\.(bin|txt)$//;
    $x =~ s/^.*\/([^\/]+)$/$1/;
    my $pat = qr/^$x.${nh}x${nv}.([0-9a-f]{8}).bin$/;
    my @bfiles = grep { /$pat/ } readdir $dh;
    if (scalar @bfiles != 1) {
        print STDERR "Expected 1 bfile, found ", scalar @bfiles, ":\n";
        print "$_\n" foreach (@bfiles);
        die;
    }
    closedir $dh;
    $balancing=$bfiles[0];
    $balancing =~ /$pat/;
    $balancing_hash = $1;
}

sub last_cp {
    my $pat;
    return undef if $force_complete;
    if ($_[0] eq 'krylov') {
        $pat = qr/^V/;
    } elsif ($_[0] eq 'mksol') {
        $pat = qr/^S/;
    } else {
        die "please enlighten me";
    }
    opendir DIR, $wdir or die "Cannot open directory `$wdir': $!\n";
    my @cp= grep /$pat/, readdir DIR;
    close DIR;
    return undef unless @cp;
    obtain_bfile();
    @cp = grep { /\.(\d+)\.$balancing_hash$/ } @cp;
    @cp = map { /\.(\d+)\.$balancing_hash$/; $1 } @cp;
    my $l = max( @cp );
    return $l;
}
sub drive {
    my $program = shift @_;

    if (! -x "$bindir/$program" && $program !~ /^:/) {
        die "No program $program found";
    }

    if ($program eq ':wipeout') {
        # Special case, really.
        my $pwd = getcwd;
        die "No directory $wdir" unless -d $wdir;
        chdir $wdir;
        die "Won't wipe cwd" if $pwd eq getcwd;
        print STDERR "Doing cleanup in $wdir\n";
        if ($show_only) {
            print "find $wdir -name '[A-Z]*' | xargs -r rm";
            print "(cd $wdir ; rm -f bw.cfg)";
        } else {
            system "find $wdir -name '[A-Z]*' | xargs -r rm";
            system "rm -f bw.cfg";
        }
        chdir $pwd;
        return;
    }

    if ($program eq ':bench') {
        &drive(":complete", @_, "end=1000000");
        return;
    }


    if ($program eq ':complete') {
        my $cp;
        &drive(":wipeout", @_) if $force_complete;
        unless (defined($cp = last_cp('mksol'))) {
            unless (defined($cp = last_cp('krylov'))) {
                &drive(":wipeout", @_);
                &drive("mf_bal", "mfile=$matrix", "out=$wdir/", $nh, $nv);
                obtain_bfile();
                push @_, "balancing=$balancing";
                &drive("${mode}_dispatch", @_, "sequential_cache_build=1");
                &drive("u64n_prep", @_);
                &drive("u64_secure", @_) unless $param->{'skip_online_checks'};
                &drive("./split", @_, "--split-y");
                &drive("${mode}_krylov", @_);
            } else {
                obtain_bfile();
                push @_, "balancing=$balancing";
                &drive("${mode}_krylov", @_, "start=$cp");
            }
            &drive("./acollect", @_, "--remove-old");
            &drive("./lingen", @_, "--lingen-threshold", 64);
            &drive("./split", @_, "--split-f");
            &drive("${mode}_mksol", @_);
        } else {
            obtain_bfile();
            push @_, "balancing=$balancing";
            &drive("${mode}_mksol", @_, "start=$cp");
        }
        &drive("u64n_gather", @_);
        # The kernel is within the space K, MK, M^2K, and so on.
        # Restricting to the last non-zero member is not necessarily
        # correct, if we happen to have pathological matrices (with
        # zero rows in the case of left nullspaces, for instance).
        # Pending an improvement to the gather code, a quick workaround
        # is to return the full thing for the characters step.
        opendir D, $wdir;
        for my $f (grep { /^(?:K\.\d+|W)\.$balancing_hash$/ } readdir D) {
            my $g = $f;
            $g =~ s/\.$balancing_hash$//;
            &drive("mf_untwistvec", "$wdir/$balancing", "$wdir/$f", "--out", "$wdir/$g");
        }
        closedir D;
        # &drive("mf_untwistvec", "$wdir/$balancing", "$wdir/W.twisted", "--out", "$wdir/W");
        return;
    }

    @_ = grep !/^tmpdir=/, @_;

    if ($program eq ':ysplit') { $program="./split"; push @_, "--split-y"; }
    if ($program eq ':fsplit') { $program="./split"; push @_, "--split-f"; }

    # Some arguments are relevant only to some contexts.
    unless ($program =~ /split$/) {
        @_ = grep !/^(?:splits?=|--split-[yf])/, @_;
    }
    unless ($program =~ /(?:krylov|mksol|dispatch)$/) {
        @_ = grep !/^ys=/, @_;
    }

    $program="$bindir/$program";

    if ($mpi_needed) {
        unshift @_, $program;
        if ($program =~ /(?:split|acollect|lingen|mf_bal)$/) {
            unshift @_, @mpi_precmd_single;
        } else {
            unshift @_, @mpi_precmd;
        }
        $program = shift @_;
    }
##     if ($mpi_split[0] * $mpi_split[1] != 1 && $program !~ /(?:split|acollect|lingen)$/) {
##         unshift @_, $program;
##         unshift @_, @mpi_precmd;
##         $program = shift @_;
##     } elsif ($program =~ /(?:split|acollect|lingen)$/ && $mpi_needed) {
##         if (defined(my $mpi_path=$ENV{'MPI'})) {
##             my $ldlp;
##             if (defined($ldlp=$ENV{'LD_LIBRARY_PATH'})) {
##                 $ldlp .= ':';
##             } else {
##                 $ldlp = '';
##             }
##             $ldlp .= "$mpi_path/lib";
##             unshift @_, 'env', "LD_LIBRARY_PATH=$ldlp";
##             $program = shift @_;
##         }
##     }

    dosystem($program, @_);
}

&drive($main, @main_args);

