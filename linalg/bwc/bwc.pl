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

if (defined($ENV{'MPI'})) {
    $mpiexec="$ENV{'MPI'}/mpiexec";
}

# Some command-line arguments have a special meaning.
#
# auto=<path> -> auto-create wdir from the mpi and thr splittings.
# hostfile=<path> -> pass to mpiexec if relevant.

while (defined($_ = shift @ARGV)) {
    if (/^(?:auto|mat)=(\S*)$/) { $matrix_dir=$1; next; }

    if (/^hostfile=(\S*)$/) { $hostfile = $1; next; }
    if (/^matpath=(\S*)$/) { $matpath=$1; next; }
    if (/^hosts=([\w:\.]+(?:,[\w:\.]+)*)$/) { @hosts=split ',', $1; next; }

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
    push @main_args, "wdir=$wdir", "matrix=mat";
}

my $env_strings = "";

sub my_setenv
{
    return if (exists($ENV{$_[0]}) && $ENV{$_[0]} eq $_[1]);
    $env_strings .= "$_[0]=$_[1] ";
    $ENV{$_[0]}=$_[1];
}


### ok -- now @main_args is something relatively useful.

sub drive {
    my $program = shift @_;

    if (! -x $program && $program !~ /^:/) {
        die "No program $program found";
    }

    if ($program eq ':wipeout') {
        # Special case, really.
        my $pwd = getcwd;
        die "No directory $wdir" unless -d $wdir;
        chdir $wdir;
        die "Won't wipe cwd" if $pwd eq getcwd;
        print STDERR "Doing cleanup in $wdir\n";
        system "/bin/rm -f [A-Z]* bw.cfg";
        chdir $pwd;
        return;
    }

    if ($program eq ':split') {
        # Special case, really.
        my $pwd = getcwd;
        chdir $wdir;
        $program="$pwd/split";
        print STDERR "Running $program ", join(' ', @splits), "\n";
        system $program, @splits;
        chdir $pwd;
        return;
    }

    if ($program eq ':complete') {
        &drive(':wipeout', @_);
        &drive('u64n/prep', @_);
        &drive('u64/secure', @_);
        &drive(':split', @_);
        &drive('u64/krylov', @_);
        return;
    }

    if ($program eq ':balance') {
        mkdir $wdir;
        my $cmd="../balance"
                . " -in $matrix_dir"
                . " -out $wdir/mat"
                . " --nslices ${nh}x${nv}"
                . " --conjugate-permutations"
                . " --square"
                ;
        print "$cmd\n";
        system $cmd;
        return;
    }

    # Some arguments are relevant only to some contexts.
    @_ = grep !/^splits?=/, @_;
    unless ($program =~ /krylov$/) {
        @_ = grep !/^ys=/, @_;
    }

    if ($mpi_split[0] * $mpi_split[1] != 1) {
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

    print STDERR "Running $env_strings$program ", join(' ', @_), "\n";
    system $program, @_;
}


&drive($main, @main_args);

