#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

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

if ($mode eq 'bmatrix') {
    my $nr=0;
    my $nc=0;
    my $curr=0;
    my @coeffs=();
    while(sysread(STDIN, my $x, 4)) {
        my $v=unpack("L",$x);
        if ($curr-- == 0) {
            $nr++;
            $curr=$v;
        } else {
            $nc = $v+1 if $v >= $nc;
            push @coeffs, [$nr, $v+1];
        }
    }
    print "var:=SparseMatrix($nr,$nc,[\n";
    print join(", ", map { "<$_->[0], $_->[1],1>" } @coeffs);
    print "]);\n";
    die unless $curr==0;
    exit;
}

if ($mode eq 'balancing') {
    sysread(STDIN, my $x, 32);
    my ($nh,$nv,$nr,$nc,$ncoeffs,$checksum,$flags) = unpack("L4QLL", $x);
    $checksum = sprintf("%04x", $checksum);
    my $txflags="";
    $txflags .= ", colperm" if $flags & 1;
    $txflags .= ", rowperm" if $flags & 2;
    $txflags .= ", padding" if $flags & 4;
    $txflags .= ", replicate" if $flags & 8;
    my $tr = $nr;
    my $tc = $nc;
    if ($flags & 4) {
        $tr = $tc = $nr > $nc ? $nr : $nc;
    }
    print "// $nr rows $nc cols, split ${nh}x${nv}, checksum $checksum$txflags\n";
    print "nr:=$tr; // originally $nr\n";
    print "nc:=$tc; // originally $nc\n";
    print "nh:=$nh;\n";
    print "nv:=$nv;\n";
    my @p=();
    while(sysread(STDIN, my $x, 4)) {
        push @p, 1+unpack("L",$x);
    }
    print "var:=[", join(", ", @p), "];\n";
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

