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
#      Benchmarks of the CFRAC/QS/SIQS/SQUFOF implementations of TIFA v0.2
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : benchmarker.pl
# Author        : Jerome Milan
# Created on    : Thu Apr  6 2006
# Last modified : Mon Feb 12 2007 by JM
#
# Version : 0.4.0
# License : GNU Lesser General Public License (LGPL) v2.1 or later
#           Copyright (C) 2006, 2007 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
#   0.4.0: Mon Feb 12 2007 by JM
#       - Completely rewritten to be (almost) completely generic. Now uses
#         the Tifa::ProgramRepository mechanism: adding new factoring programs
#         is now possible without editing this script. Also provides a generic
#         benchmarking framework as long as the program's arguments are given
#         in a file. If not, this script will make one from generated composite
#         numbers (the last relict from its non-generic, integer factorization-
#         oriented days).
#
#   0.3.3: Thu Nov  2 2006 by JM
#       - Added finer timings for different stages of the factoring algorithms
#
#   0.3.2: Thu Oct 19 2006 by JM
#       - Added filter option (select job to run according to an expression)
#
#   0.3.1: Tue Oct 4 2006 by JM
#       - Renamed as benchmarker.pl
#       - Added QS and SIQS support
#       - Some cleanup
#
#   0.3: Tue Jul 25 2006 by JM
#       - Modified number generation. Now use the external
#         Tifa::NumberGenerator module.
#       - Modified generated script to use the external module
#         Tifa::DataDescriptor.
#       - Now uses the Tifa::SimpleConfigReader to read the configuration file
#         WARNING: Configuration file format changed!
#       - WARNING: Data file format changed to use Tifa::DataDescriptor
#
#   0.2.1: Thu Jul 20 2006 by JM
#       - Modified slightly the way the numbers are generated (the previous code
#         was generating numbers 1 bit smaller than desired)
#
#   0.2: Thu May 18 2006 by JM
#       - Code cleanup
#       - No more nb_step_cfrac parameter
#       - Got rid of confusing default values (they were not helpful anyway)
#       - Added option to give the parameter values in the configuration
#         file as a list, eg: (128, 256, 512, 1024)
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# This script provides an almost generic framework for benchmarking programs.
# It takes a configuration file as input listing all parameters relevant to the
# program to bench, and another file listing the arguments the program will be
# benched with.
#
# If no argument file is provided, the benchmarker.pl script will make one by
# generatic composite numbers, a feature clearly inherited from its non-generic,
# integer factorization-oriented days.
#
# In any case, this script will create (but NOT execute!) the real benchmarking
# Perl script which will execute the program for a whole range of parameters as
# given in the configuration file.
#
# _WARNING_: This script needs the GMP Perl module GMP::Mpz. This module is
#            available with each distribution of GMP, but is not build nor
#            installed by default.
#-------------------------------------------------------------------------------

use strict;
use Fcntl;
use Getopt::Long;
use File::Basename;
use File::Copy;

use Tifa::NumberGenerator;
use Tifa::SimpleConfigReader;

eval "use GMP::Mpz qw(:all); 1" or die "ERROR: Module GMP::Mpz not found!";

use Tifa::Bencher;

#
# Fetch informations on available options and their types from the
# centralized Tifa::Repository
#
my @getopt_strings  = Tifa::ProgramRepository::get_getopt_strings();

#
# Options/arguments exclusively passed via the command line
#
my $config_file;   # name of the configuration file
my $args_file;     # name of the argument file
my $macro_file;    # name of generated macro file (the real bench macro)
my $output_dir;    # output directory for bench results and traces
my $help;

#
# Provide some default values
#
$config_file = "benchmarker.conf";
$macro_file  = "bench.pl";
$output_dir  = "bench_results";

my %param_values = ();

#
# Get the list of options from the command line and process them
#
GetOptions(
    \%param_values, @getopt_strings,
    "conf=s"   => \$config_file,
    "args=s"   => \$args_file,
    "macro=s"  => \$macro_file,
    "outdir=s" => \$output_dir,
    "filter=s" => \$param_values{"filter"},
    "mode=s"   => \$param_values{"mode"},
    "help"     => \$help
);

if (defined $help ) {
    print_help();
    exit();
}

#
# Performs some basic checks on the command line options passed
#
check_options_validity();

#
# Read the configuration file and extract parameter values
#
my $confreader = new Tifa::SimpleConfigReader();

$confreader->accept_unknown_param();
$confreader->read_config_file($config_file);

my %param_hash = $confreader->get_parameter_hash();

#
# Set the values of the parameters to the one just read, but give precedence
# to the values passed directly on the command line over the one defined in the
# configuration file
#
foreach my $name (keys %param_hash) {
    if (!defined $param_values{$name} && defined $param_hash{$name}) {
        $param_values{$name} = $param_hash{$name};
    }
}

#
# Delegate all the dirty work to the Tifa::Bencher module
#
my $bencher = new Tifa::Bencher;

my $argfile;
if (! defined $args_file) {
    #
    # Here is the ugly part that makes this script not so generic anymore: if
    # no argument file is provided, make a temporary one with composite numbers
    # as arguments: clearly a relict of the integer-factorization-only days
    # of benchmarker.pl
    #
    $argfile = sprintf(".tmp_%x%x.txt", int(rand(2**31)), int(rand(2**31)));
    generate_args_file($argfile);
} else {
    $argfile = $args_file;
}

$bencher->set_args_file($argfile);

$bencher->set_output_dir($output_dir);

$bencher->set_param_values(\%param_values);

#
# benchmarks can be constrained by the use of a conditionnal expression, what
# is called "filter" here
#
$bencher->set_filter($param_values{"filter"});

$bencher->create_bench_macro($macro_file);

#
# Delete the temporary argument file if it was generated by this script
#
unlink($argfile) if (! defined $args_file);

#-------------------------------------------------------------------------------
#   Subroutines
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
sub print_usage {
    my $name = basename($0);
    print("Usage:\n");
    print("------\n\n");

    printf("%15s [--conf   <config_file>]      [--macro <output_macro>]\n",
            $name);
    printf("%15s [--outdir <output_directory>] [--args  <argument_file>]\n","");
    printf("%15s [--args   <output_directory>] [--help]\n\n","");
}
#-------------------------------------------------------------------------------
sub print_help {
    my $name = basename($0);

    my $title = "$name - A (mostly) generic benchmark script";
    my $line  = '-' x length($title);
    print << "EOF";

$title
$line

EOF
    print_usage();

    print << "EOF";
The $name script provides an (almost!) generic framework for
benchmarking programs. It takes a configuration file as input listing the
program to bench and all it's relevant parameters, and another file listing the
arguments the program will be benched with. If no argument file is provided, the
benchmarker.pl script will make one by generatic composite numbers (a feature
clearly inherited from its non-generic, integer factorization-oriented days).

In any case, this script will create (but NOT execute!) the "real" benchmarking
Perl script which will (once manually launched) execute the program to benchmark
for the whole range of parameters given in the configuration file.

Note that this script can be used to bench any kind of programs as long as
they are wrapped in a Tifa::Program object and declared in the
Tifa::ProgramRepository module. Users willing to extend the list of available
programs should then have a closer look at these two modules.

Usage of this script can be hard to grasp, so the easiest way to get used to
it is probably to just read the provided benchmark.conf configuration file in
addition to this manual.

General parameters/options:
---------------------------

  --args=s
      (optional)
      Name of the argument file listing the arguments to be passed to the
      program to benchmark. If none are provided, composite numbers will be
      generated and used as argument (which is fine with TIFA's factorization
      programs). The number generation can be controlled by the
      bit_length_of_n and nprime_factors_in_n fields in the configuration
      file. See the provided benchmark.conf configuration file for more
      information.
      Default value: none

  --conf=s
      (optional, default value used if none provided)
      Name of the configuration file listing the values of the parameters
      relevant to the program to benchmark. This program will be benchmarked
      using each possible combination of the parameter values provided in this
      file. The configuration file syntax should match the one defined by the
      Tifa::SimpleConfigReader module. See the man page associated to that
      module for more information. Reading the provided benchmark.conf
      configuration file is also highly advised.
      Default value: $config_file

  --macro=s
      (optional, default value used if none provided)
      Name of the "real" benchmarking perl script that will be generated by
      $name. Note that this generated script has to be launched
      manually to proceed to the benchmarks.
      Default value: $macro_file

  --outdir=s
      (optional, default value used if none provided)
      Name of the directory where bench results and trace files will be saved.
      Default value: $output_dir

  --mode=s
      (optional, value in configuration file used if none provided)
      Program's mode to use.

  --filter=s
      (optional, value in configuration file used if none provided)
      Impose a condition on the set of parameters. A run will be performed only
      if this expression evaluates to non zero. This condition should be written
      as Perl 5.8 code.

  --help
      Prints this short manual.

See also:
---------

  The documentation for the Tifa::Program and Tifa::ProgramRepository modules.

EOF

  exit();
}
#-------------------------------------------------------------------------------
sub check_options_validity {

    if (! -e $config_file) {
        print("ERROR: Cannot find configuration file $config_file\n");
        exit();
    }
    if (! -r $config_file) {
        print("ERROR: Cannot read configuration file $config_file ");
        print("(bad permissions?)\n");
        exit();
    }
    if (-e $macro_file) {
        print("ERROR: Macro file $macro_file already exists\n");
        exit();
    }
    if (-e $output_dir) {
        if (! -d $output_dir) {
            print("ERROR: $output_dir is not a directory\n");
            exit();
        }
    }
    if ($args_file) {
        if (! -e $args_file) {
            print("ERROR: Cannot find argument file $args_file\n");
            exit();
        }
    }
}
#-------------------------------------------------------------------------------
sub generate_args_file {
    my $outfile = $argfile;

    my $generator = new Tifa::NumberGenerator();
    $generator->dont_use_random();

    open(OUT, ">$outfile")
        or croak( "ERROR: generate_args_file: Cannot create file $outfile");

    foreach my $nb_factors (@{$param_values{"nprime_factors_in_n"}}) {

        foreach my $length (@{$param_values{"bit_length_of_n"}}) {

            my $n = $generator->generate_composite($length, $nb_factors);

            print OUT "#\n";
            print OUT "# Bitsize: $length";
            print OUT " Number of prime factors: $nb_factors\n";
            print OUT "#\n";
            print OUT "$n\n";
        }
    }
    close(OUT);
}
#-------------------------------------------------------------------------------

