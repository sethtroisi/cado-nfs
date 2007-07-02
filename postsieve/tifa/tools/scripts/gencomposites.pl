#!/usr/bin/perl -w

#
# Copyright (C) 2006, 2007 INRIA (French National Institute for Research in
# Computer Science and Control)
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
#

#-------------------------------------------------------------------------------
#                     Generate a list of composite numbers
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : gencomposites.pl
# Author        : Jerome Milan
# Created on    : Fri Feb  20 2007
# Last modified : Fri Feb  20 2007 by JM
#
# Version : 0.1.0
# License : GNU Lesser General Public License (LGPL)
#           Copyright (C) 2006, 2007 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# This script generate a file containing a list of composite numbers of
# different bit length.
#-------------------------------------------------------------------------------

use strict;
use Getopt::Long;
use File::Basename;
use Tifa::NumberGenerator;

my $minsize   = 60;
my $maxsize   = 100;
my $persize   = 10;
my $incremt   = 10;
my $outfile   = "numbers.txt";
my $nfactors  = 2;
my $norandom  = 0;
my $nocomment = 0;
my $help;

GetOptions(
    "minsize=i"  => \$minsize,
    "maxsize=i"  => \$maxsize,
    "persize=i"  => \$persize,
    "out=s"      => \$outfile,
    "nfactors=i" => \$nfactors,
    "incremt=i"  => \$incremt,
    "norandom"   => \$norandom,
    "nocomment"  => \$nocomment,
    "help"       => \$help
);

if ($help) {
    print_help();
    exit(-1);
}

check_options_validity();

my $gen = new Tifa::NumberGenerator();
$gen->use_random() unless ($norandom);

#
# Create the output file and write composites integers, one per line...
#
open(OUT, ">>$outfile") or die("Cannot open $outfile!\n");
for (my $size = $minsize; $size <= $maxsize; $size += $incremt) {
    if (! $nocomment) {
        print OUT "#\n";
        print OUT "# Numbers of $size bits with $nfactors prime factors\n";
        print OUT "#\n";
    }
    for my $index (1 .. $persize) {
        my $number = $gen->generate_composite($size, $nfactors);
        print OUT "$number\n";
    }
}
close(OUT);

#------------------------------------------------------------------------------
#                            Subroutines
#------------------------------------------------------------------------------
sub check_options_validity {
    my $error = 0;
    if ($minsize <= 0) {
        print("ERROR: minsize must be greater than 0.\n");
        $error++;
    }
    if ($maxsize <= 0) {
        print("ERROR: maxsize must be greater than 0.\n");
        $error++;
    }
    if ($persize <= 0) {
        print("ERROR: persize must be greater than 0.\n");
        $error++;
    }
    if ($nfactors <= 0) {
        print("ERROR: nfactors must be greater than 0.\n");
        $error++;
    } elsif ($nfactors > $minsize/3) {
        print("ERROR: nfactors must be lesser than minsize/3.\n");
        $error++;
    }
    if ($incremt <= 0) {
        print("ERROR: incremt must be greater than 0!\n");
        $error++;
    }
    if ($maxsize < $minsize) {
        print("ERROR: maxsize must be greater than minsize.\n");
        $error++;
    }
    if (-e $outfile) {
        print("WARNING: output file $outfile already exists. Appending ");
        print("new numbers.\n");
    }
    if ($error) {
        exit(-1);
    }
}
#------------------------------------------------------------------------------
sub print_help {
    my $name  = basename($0);
    my $title = "$name - Generate composite numbers";
    my $line  = '-' x length($title);
    print << "EOF";

$title
$line

The $name script generates composite numbers for bit sizes in a given
range, and writes them in a text file. These composites are obtained by
multiplying (roughly) equal-sized (pseudo-)prime factors.

Usage:
------

    $name [options]

Available options:
------------------

  --minsize=i
      Minimal size (in bits) of the composites to generate.
      Default: $minsize

  --maxsize=i
      Maximal size (in bits) of the composites to generate.
      Default: $maxsize

  --incremt=i
      Bit size increment. The size of the generated composites will spans from
      'minsize' to 'maxsize' by increments of 'incremt'.
      Default: $incremt

  --persize=i
      Number of composites to generate for each given size.
      Default: $persize

  --nfactors=i
      Number of prime factors of the generated composite numbers.
      Default: $nfactors

  --norandom
      Do not generate prime factors randomly: sequentially use the smallest
      primes of the 'good size' (i.e. size(composite) / 'nfactors') making
      sure no prime is used twice.

  --nocomment
      Do not write comments (i.e. size of composite and number of factors) in
      the output file.

  --out=s
      Output file name. Composite numbers will be appended at the end of the
      output file if it already exists.
      Default: $outfile

  --help
      Print this help message.

EOF
}
