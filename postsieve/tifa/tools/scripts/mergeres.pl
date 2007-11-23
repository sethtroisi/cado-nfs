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
#          Merge data from result files into one output result file
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
#
# File          : mergeres.pl
# Author        : Jerome Milan
# Created on    : Mon Feb 19 2007
# Last modified : Mon Feb 19 2007 by JM
#
# Version : 0.1.0
# License : GNU Lesser General Public License (LGPL) v2.1 or later
#           Copyright (C) 2006, 2007 INRIA
#-------------------------------------------------------------------------------
# History :
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# This script can be used after executing a benchmark macros created by
# the benchmarker.pl script. It reads the two result files from the benchmarks
# to create a merged result file. This script is generic and does not make
# assumptions on the result files other than being in the proper format (as
# defined by the module Tifa::DataDescriptor). Note however than each entry
# in the input files will contribute to one entry in the output file (i.e.
# individual entry data are not merged together).
#-------------------------------------------------------------------------------

use strict;
use Fcntl;
use Getopt::Long;
use File::Basename;
use File::Copy;

use Tifa::DataDescriptor;

#
# Options/arguments exclusively passed via the command line.
#
my @input_files = ();
my $output_file;
my $help;

GetOptions(
    "in=s"  => \@input_files,
    "out=s" => \$output_file,
    "help"  => \$help
);

if (defined $help ) {
    print_help();
    exit();
}

check_options_validity();

my @entries = ();
my %descr_hash = ();
#
# Read files and extract descriptions and data entries
#
foreach my $infile (@input_files) {
    open(my $IN, "<$infile") or die("Cannot open $infile\n");

    my $descr_in = new Tifa::DataDescriptor();

    $descr_in->load_descriptions($IN);

    my @params_in = $descr_in->get_all_field_names();

    foreach my $param (@params_in) {
        $descr_hash{$param} = $descr_in->get_field_description($param);
    }
    while (my %entry_in = $descr_in->read_next_data_entry($IN)) {
       push(@entries, \%entry_in);
    }
    close($IN);
}
#
# Create ouput file and write data for all entries
#
open(my $OUT, ">$output_file") or die("Cannot open $output_file\n");

my $descr_out = new Tifa::DataDescriptor();

my @params = sort keys %descr_hash;

foreach my $param (@params) {
    $descr_out->add_field_description($param, $descr_hash{$param});
}
$descr_out->write_descriptions($OUT);

foreach my $entry (@entries) {
   my @data_entry = ();

   foreach my $param (@params) {
       #
       # Replace field with '-1' if entry doesn't contain such a field
       #
       if (defined $$entry{$param}) {
           push(@data_entry, $$entry{$param});
       } else {
           push(@data_entry, -1);
       }
   }
   $descr_out->write_data_entry($OUT, @data_entry);
}

close($OUT);

#-------------------------------------------------------------------------------
#   Subroutines
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
sub print_help {
    my $name = basename($0);

    my $title = "$name - Merge data from result files";
    my $line  = '-' x length($title);
    print << "EOF";

$title
$line

Usage:
------

EOF

    printf("%15s --in <input_file_1> --in <input_file_2>\n", $name);
    printf("%15s [--in <other_input_file>]* --out <output_file>\n", '');
    printf("%15s [--help]\n", '');

print << "EOF";

The $name script can be used after executing benchmark macros created by
the benchmarker.pl script. It reads the result files from the benchmarks
to create a merged result file. This script is generic and does not make
assumptions on the result files other than being in the proper format (as
defined by the module Tifa::DataDescriptor). Note however than each entry
in the input files will contribute to one entry in the output file (i.e.
individual entry data are not merged together).

General parameters/options:
---------------------------

  --in=s
      (mandatory)
      Names of the input result files as generated by the benchmarking script
      created by benchmarker.pl.

  --out=s
      (mandatory)
      Name of the output result file. This file will contain information
      from the original input files.

See also:
---------

  The documentation for the benchmarker.pl script and the Tifa::DataDescriptor
  module.

EOF

  exit();
}
#-------------------------------------------------------------------------------
sub check_options_validity {

    if ((scalar @input_files) == 0) {
        print("ERROR: No input file specified\n");
        exit();
    }
    if ((scalar @input_files) == 1) {
        print("ERROR: Only one input file specified\n");
        exit();
    }
    foreach my $file (@input_files) {
        if (! -e $file) {
            print("ERROR: Cannot find input file $file\n");
            exit();
        }
    }
    if (-e $output_file) {
        print("ERROR: Output file $output_file already exists\n");
        exit();
    }
}
#-------------------------------------------------------------------------------
