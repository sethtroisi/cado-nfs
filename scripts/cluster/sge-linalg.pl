#!/usr/bin/perl -w

# This script is specific to the Oracle SGE (Sun Grid Engine) scheduler.
# It serves as a helper to submit a bwc (``:complete'', for now) job.

use strict;
use warnings;

# At some point we'll substitute this with some better mechanism... This
# gives the path where the (cmake-substituted) bwc.pl is found.
my $bwc_bindir="$ENV{HOME}/cado-nfs/build/tom_pouce/linalg/bwc";
my $matrix;
my $wdir;

# The node numbers on which we are going to run. These get substituted to
# real names by nodenumber_to_nodename below. This is just a default
# value, and can be overriden by a comma-separated list of node numbers
# on the command line, as in:
# ./sge-linalg.pl [...] nodes=1,2,3,4
#
my @nodes=qw/1 2 3 4/;

# mpi and thr are potentially set automatically, but is is legit to set
# them on the command-line as well.
my ($mpi,$thr);

while (defined($_=shift @ARGV)) {
        if (/^matrix=(.*)$/) {
                $matrix=$1;
                die "matrix file $matrix is not readable" unless -r $matrix;
                next;
        }

        if (/^wdir=(.*)$/) {
                $wdir=$1;
                die "working directory $wdir already exists" if -e $wdir;
                next;
        }
        if (/^nodes=(\d+(?:,\d+)*)$/) {
            @nodes = split ',',  $1;
            next;
        }
        if (/^mpi=(\d+x\d+)$/) { $mpi = $1; next; }
        if (/^thr=(\d+x\d+)$/) { $thr = $1; next; }
        die "unexpected arg $_\n";
}

die "wdir= is mandatory" unless $wdir;
die "matrix= is mandatory" unless $matrix;
die "nodes= is mandatory" unless scalar @nodes;

mkdir $wdir;

if (defined($mpi)) {
    $mpi =~ /^(\d+)x(\d+)$/ or die;
    my $x = $1 * $2;
    my $nnodes = scalar @nodes;
    die "Cannot use $nnodes nodes with mpi=$mpi (which requires $x nodes)\n" if $x != $nnodes;
}

my @args = qw|
        -cwd
        -N bwc_linalg
|;

push @args, "-e", "$wdir/stderr";
push @args, "-o", "$wdir/stdout";


# Site-specific:
#
# Add the magic commands to the submission so that we reserve exclusively
# the nodes we want.
my $cores_per_node = 12;
sub nodenumber_to_nodename { return sprintf("node%03d.cm.cluster", $_); }
my $qlist = join(",", map { "batch.q\@" . nodenumber_to_nodename($_); } @nodes);
push @args, "-q", $qlist;
push @args, "-pe", "openmpi", ($cores_per_node * scalar @nodes);

# This is a simple-minded strategy for now. We just pick a split which is
# as square as it can be.
sub best_split_for_size {
    my $n = shift;
    my $x = 1;
    my @possible = ($x);
    while ($x*$x <= $n) {
        my $y = int($n/$x);
        push @possible, $x if $x*$y == $n;
        $x++;
    }
    $x = $possible[$#possible];
    my $y = int($n/$x);
    my @reverse = reverse @possible;
    if ($x == $y) {
        shift @reverse;
        @reverse = map { int($n/$_); } @reverse;
    }
    push @possible, @reverse;
    return "${x}x${y}", @possible;
}

if (!defined($mpi)) {
    my $n = scalar @nodes;
    my @all_possible;
    ($mpi, @all_possible) = best_split_for_size($n);
    print STDERR "Possible splits for $n nodes: ", join(", ", map { $_ .
        "x" . int($n/$_); } @all_possible), "\n";
    print STDERR "Picking mpi=$mpi\n";
}

if (!defined($thr)) {
    my $n = $cores_per_node;
    my @all_possible;
    ($thr, @all_possible) = best_split_for_size($n);
    print STDERR "Possible splits for $n cores: ", join(", ", map { $_ .
        "x" . int($n/$_); } @all_possible), "\n";
    print STDERR "Picking thr=$thr\n";
} else {
    my $n = $cores_per_node;
    $thr =~ /^(\d+)x(\d+)$/ or die;
    my $x = $1 * $2;
    die "Cannot use $n cores with thr=$thr (which requires at least $x cores)\n" if $x > $n;
}


push @args, "$bwc_bindir/bwc.pl", qw/:complete seed=1 nullspace=left
mm_impl=bucket interleaving=0 interval=100 mn=64 shuffled_product=1/;
push @args, "mpi=$mpi";
push @args, "thr=$thr";
push @args, "matrix=$matrix";
push @args, "wdir=$wdir";
# Should no longer be needed, in fact.
push @args, "bwc_bindir=$bwc_bindir";

unshift @args, "qsub";

print join(" ", @args), "\n";

# matrix=/scratch/morain/CADO/M307/M307.small.bin wdir=/scratch/morain/CADO/M307/M307.bwc
# qsub -pe openmpi 48 /home/tanc/morain/cado-nfs/build/tom_pouce/linalg/bwc/bwc.pl :complete seed=1 thr=1x1 mpi=2x2 matrix=/scratch/morain/CADO/M307/M307.small.bin nullspace=left mm_impl=bucket interleaving=0 interval=100 mn=64 wdir=/scratch/morain/CADO/M307/M307.bwc shuffled_product=1 bwc_bindir=/home/tanc/morain/cado-nfs/build/tom_pouce/linalg/bwc >> /scratch/morain/CADO/M307/M307.bwc.log 2>&1
# qsub -N linalg -cwd -e foo.err -o foo.out -sync y -q "batch.q@node001.cm.cluster,batch.q@node002.cm.cluster,batch.q@node003.cm.cluster,batch.q@node004.cm.cluster" -pe openmpi 48 /home/tanc/morain/cado-nfs/build/tom_pouce/linalg/bwc/bwc.pl :complete seed=1 thr=1x1 mpi=2x2 matrix=/scratch/morain/CADO/M307/M307.small.bin nullspace=left mm_impl=bucket interleaving=0 interval=100 mn=64 wdir=/scratch/morain/CADO/M307/M307.bwc shuffled_product=1 bwc_bindir=/home/tanc/morain/cado-nfs/build/tom_pouce/linalg/bwc >> /scratch/morain/CADO/M307/M307.bwc.log 2>&1

