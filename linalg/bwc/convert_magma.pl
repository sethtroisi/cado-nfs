#!/usr/bin/perl

use warnings;
use strict;

# Converts one bwc file as given on stdin to magma-parseable text on
# stdout

my $mode = shift @ARGV;

if ($mode eq 'matrix') {
    die unless defined($_=<>);
    die unless /^(\d+) (\d+)$/;
    my $nr=$1;
    my $nc=$2;
    my $i = 0;
    my $ne = 0;
    my $needcomma = 0;
    print "var:=SparseMatrix($nr,$nc,[\n";
    while (<>) {
        s/^(\d+)// or die;
        my $l = $1;
        if ($l == 0) { $i++; next; }
        print ",\n" if $needcomma;
        $needcomma=0;
        for(my $j = 0;$j < $l;$j++) {
            s/^\s+(\d+)// or die;
            print ", " if $needcomma;
            $ne++;
            my $i1 = $i + 1;
            my $j1 = $1 + 1;
            print "<$i1,$j1,1>";
            $needcomma=1;
        }
        die unless /^\s*$/;
        $i++;
    }
    print "]);\n";
    if ($i != $nr) {
        print STDERR "Matrix nrows was wrong (read $i, expected $nr)\n";
    }
    exit;
}

if ($mode eq 'permutation') {
    # Dump a list of unsigned ints
    my @p=();
    while(sysread(STDIN, my $x, 4)) {
        my $v = unpack("L",$x);
        push @p, $v+1;
    }
    print "var:=[",join(',',@p),"];\n";
    exit;
}

if ($mode eq 'x') {
    die unless defined($_=<>);
    die unless /^(\d+)$/;
    my $nx = 1;
    my @p=();
    while(<>) {
        my @xs = split ' ', $_;
        @xs = map { $_+1; } @xs;
        my $xstring=join(",",@xs);
        push @p,"[$xstring]";
    }
    print "var:=[",join(',',@p),"];\n";
    exit;
}

if ($mode eq 'vector') {
    # Dump a list of uint64_t's ; this version is robust on 32-bits.
    my @p=();
    while(sysread(STDIN, my $x, 8)) {
        my @v = unpack("L2",$x);
        my $sv = join(',',@v);
        my $v = "Seqint([$sv],2^32)";
        push @p, $v;
    }
    print "var:=[",join(',',@p),"];\n";
    exit;
}

