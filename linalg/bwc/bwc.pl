#!/usr/bin/perl
use warnings;
use strict;

use POSIX qw/getcwd/;

# This companion program serves as a helper for running bwc programs.

# MPI=/opt/mpich2-1.1a2-1.fc10.x86_64/usr/bin ./bwc.pl :complete mat=c72 interval=256 splits=64 mn=64 ys=0..64 mpi=4x4

# $program denotes either a simple command, or something prepended by ':'
# which indicates something having a special meaning for the script.
my $main = shift @ARGV;

# First we're going to parse globally the list of arguments, to see if
# there's something to infer.

my $matpath="$ENV{HOME}/Local/mats";
my @main_args = ();
my @mpi_split=(1,1);
my @thr_split=(1,1);
my $matrix_dir;
my $hostfile;
my $mpiexec='mpiexec';
my $wdir;
my $mn;
my @splits=();
my @hosts=();
my $mode='u64';
my $show_only=0;

if (defined($ENV{'MPI'})) {
    $mpiexec="$ENV{'MPI'}/mpiexec";
}

my $bindir;
if (!defined($bindir=$ENV{'BWC_BINDIR'})) {
    $bindir=$0;
    $bindir =~ s{/[^/]*$}{};
    if ($bindir eq $0) {
        $bindir='.';
    }
}

# Some command-line arguments have a special meaning.
#
# auto=<path> -> auto-create wdir from the mpi and thr splittings.
# hostfile=<path> -> pass to mpiexec if relevant.

while (defined($_ = shift @ARGV)) {
    if (/^(?:auto|mat|matrix)=(\S*)$/) { $matrix_dir=$1; next; }

    if (/^hostfile=(\S*)$/) { $hostfile = $1; next; }
    if (/^matpath=(\S*)$/) { $matpath=$1; next; }
    if (/^hosts=([\w:\.-]+(?:,[\w:\.-]+)*)$/) {
        my @nh=split ',', $1;
        push @hosts, @nh; 
        next;
    }
    if (/^mode=(\w+)$/) { $mode=$1; next; }
    if (/^(-d|--show)$/) { $show_only=1; next; }

    # These will be passed anyway.
    push @main_args, $_;
    
    if (/^wdir=(\S*)$/) { $wdir=$1; next; }
    if (/^mn=(\S*)$/) { $mn=$1; next; }
    if (/^splits?=(\d+(?:,\d+)*)$/) { @splits=split ',', $1; next; }
    if (/^mpi=(\d+)x(\d+)$/) { @mpi_split = ($1, $2); next; }
    if (/^thr=(\d+)x(\d+)$/) { @thr_split = ($1, $2); next; }
}

my $nh = $mpi_split[0] * $thr_split[0];
my $nv = $mpi_split[1] * $thr_split[1];

if (!defined($wdir) && defined($matrix_dir)) {
    if ($matrix_dir !~ m{/} && defined($matpath)) {
        $matrix_dir = "$matpath/$matrix_dir";
    }
    $matrix_dir =~ s/~/$ENV{HOME}/g;
    $wdir="$matrix_dir-${nh}x${nv}";
    if ($main !~ /^:(?:balance|complete)$/ && !-d $wdir) {
        die "No directory $wdir found";
    }
    push @main_args, "wdir=$wdir";
}

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


### ok -- now @main_args is something relatively useful.

sub drive {
    my $program = shift @_;

    if (! -x "$bindir/$program" && $program !~ /^:/) {
        die "No program $program found";
    }

    if ($program eq ':wipeout') {
        # Special case, really.
        my $pwd = getcwd;
        if ($show_only && ! -d $wdir) {
            print "mkdir $wdir\n";
        }
        die "No directory $wdir" unless -d $wdir;
        chdir $wdir;
        die "Won't wipe cwd" if $pwd eq getcwd;
        print STDERR "Doing cleanup in $wdir\n";
        if ($show_only) {
            print "/bin/rm -f [A-Z]* bw.cfg\n";
        } else {
            system "/bin/rm -f [A-Z]* bw.cfg";
        }
        chdir $pwd;
        return;
    }

    if ($program eq ':complete') {
        if (! -d $wdir) {
            &drive(':balance', @_);
        }
        &drive(":wipeout", @_);
        &drive("u64n/prep", @_);
        &drive("u64/secure", @_);
        &drive("./split", @_, "--split-y");
        &drive("$mode/krylov", @_);
        &drive("./acollect", @_, "--remove-old");
        &drive("lingen/lingen", @_, "--lingen-threshold", 64);
        &drive("./split", @_, "--split-f");
        &drive("$mode/mksol", @_);
        &drive("u64n/gather", @_);
        return;
    }

    if ($program eq ':balance') {
        mkdir $wdir;
        my $cmd="$bindir/../balance";
        my @args= (
                "-in", $matrix_dir,
                "-out", "$wdir/mat",
                "--nslices", "${nh}x${nv}",
                "--conjugate-permutations",
                "--square");
        dosystem($cmd, @args);
        return;
    }

    if ($program eq ':ysplit') { $program="./split"; push @_, "--split-y"; }
    if ($program eq ':fsplit') { $program="./split"; push @_, "--split-f"; }

    # Some arguments are relevant only to some contexts.
    unless ($program =~ /split$/) {
        @_ = grep !/^(?:splits?=|--split-[yf])/, @_;
    }
    unless ($program =~ /(?:krylov|mksol)$/) {
        @_ = grep !/^ys=/, @_;
    }

    $program="$bindir/$program";

    if ($mpi_split[0] * $mpi_split[1] != 1 && $program !~ /(?:split|acollect|lingen)$/) {
        unshift @_, $program;
        # Need hosts.
        if (defined($hostfile)) {
            my_setenv 'HYDRA_HOST_FILE', $hostfile;
        } elsif (scalar @hosts) {
            open F, ">$wdir/HOSTS";
            for my $h (@hosts) { print F "$h\n"; }
            close F;
            my_setenv 'HYDRA_HOST_FILE', "$wdir/HOSTS";
        } else {
            my_setenv 'HYDRA_USE_LOCALHOST', 1;
        }
        unshift @_, '-n', $mpi_split[0] * $mpi_split[1];
        $program=$mpiexec;
    }

    dosystem($program, @_);
}


&drive($main, @main_args);

