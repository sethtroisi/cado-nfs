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
#        Extracts timings from trace files and generates a result file
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : extractres.pl
# Author        : Jerome Milan
# Created on    : Fri Feb  9 2007
# Last modified : Fri Feb  9 2007 by JM
#
# Version : 0.1.0
# License : GNU Lesser General Public License (LGPL) v2.1 or later
#           Copyright (C) 2006, 2007 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# This script can be used after executing a benchmark macro created by
# the benchmarker.pl script. It reads the original result file and the trace
# files generated during the benchmarks to create a more complete result file
# by adding extra information found in the trace files. This script is not
# genering and is tailored to be used if and only if a TIFA program has been
# benchmarked since it expects to find specific information in the trace files.
#-------------------------------------------------------------------------------

use strict;
use Fcntl;
use Getopt::Long;
use File::Basename;
use File::Copy;

use Tifa::DataDescriptor;

eval "use GMP::Mpz qw(:all); 1" or die "ERROR: Module GMP::Mpz not found!";

#
# Options/arguments exclusively passed via the command line.
#
my $input_file;
my $output_file;
my $trace_dir;
my $nb_factors = -1;
my $help;

#
# Descriptions for the extra fields to be included in the output data file.
#
my %descr_hash_extra = (
    "nbits"       => "Number of bit of number to factor",
    "nfactors"    => "Number of prime factors composing number to factor",
    "n"           => "Number to factor",
    "inittime"    => "System time spent computing startup data",
    "cleantime"   => "System time spent cleaning up",
    "totaltime"   => "Total system time spent by factorization function",
    #
    # Extra timings for congruences of square methods
    #
    "relcoltime"  => "System time spent collecting relations",
    "relgentime"  => "System time spent generating relations",
    "relseltime"  => "System time spent selecting relations",
    "ptdivtime"   => "System time spent trial dividing residues",
    "factime"     => "System time spent factoring unfactored parts of residues",
    "lsrtime"     => "System time spent solving linear system",
    "deducetime"  => "System time spent deducing factors from collected data",
    #
    # Extra timings for SQUFOF
    #
    "fwdcyctime"  => "System time spent forward cycling to find a proper form",
    "invsqrttime" => "System time spent computing inverse square root of form",
    "revcyctime"  => "System time spent reverse cycling to find a factor",
    #
    # Extra timings for trial division
    #
    "tdivtime"    => "System time spent trial dividing the number to factor"
);

GetOptions(
    "in=s"        => \$input_file,
    "out=s"       => \$output_file,
    "tracedir=s"  => \$trace_dir,
    "nfactors=i"  => \$nb_factors,
    "help"        => \$help
);

if (defined $help ) {
    print_help();
    exit();
}

check_options_validity();

#
# Open the original result file and extract information via a DataDescriptor.
#
open(my $inhandle, "<$input_file") or die("Cannot open $input_file\n");

my $descr_in = new Tifa::DataDescriptor();

$descr_in->load_descriptions($inhandle);

my @params_in     = $descr_in->get_all_field_names();
my %descr_hash_in = ();

foreach my $param (@params_in) {
    $descr_hash_in{$param} = $descr_in->get_field_description($param);
}

#
# Create the new result file that will merge results from the original file
# and new information extracted from the trace files.
#
open(my $outhandle, ">$output_file") or die("Cannot open $output_file\n");

my $descr_out = new Tifa::DataDescriptor();

my %descr_hash_out = ();

@descr_hash_out{keys %descr_hash_extra} = values %descr_hash_extra;
@descr_hash_out{keys %descr_hash_in}    = values %descr_hash_in;

#
# Suppress the 'argmt' field as in our case, it is indeed a number to factor
# and it will be replaced by a 'n' field
#
delete $descr_hash_out{"argmt"};

foreach my $key (sort keys %descr_hash_out) {
    $descr_out->add_field_description($key, $descr_hash_out{$key});
}
$descr_out->write_descriptions($outhandle);

my @params_out = $descr_out->get_all_field_names();

my @data_entry;

my @trace_files = get_trace_filenames();

my $nfiles = scalar @trace_files;

my @entries = ();

while (my %entry_in = $descr_in->read_next_data_entry($inhandle)) {
   push(@entries, \%entry_in);
}

my $nentries = scalar @entries;

if ($nentries != $nfiles) {
    print("ERROR: found $nfiles trace files but there is $nentries entries!\n");
    exit();
}

for (my $ifile=0; $ifile <= $#trace_files; $ifile++) {
   @data_entry = ();
   #
   # Read next trace file and extract new information.
   #
   # _WARNING_: The trace files are sorted according to their number given
   #            in their filenames. It is _assumed_ that this order matches
   #            the order of the entries in the original result file (which
   #            will always be true unless the result file has been manually
   #            edited).
   #
   print("Processing file ", $ifile + 1, " of $nfiles...");
   my $entry_out_ref = make_new_entry($trace_files[$ifile], $entries[$ifile]);
   print("\n");

   foreach my $param (@params_out) {
       push(@data_entry, $$entry_out_ref{$param});
   }
   $descr_out->write_data_entry($outhandle, @data_entry);
}

close($inhandle);
close($outhandle);

#-------------------------------------------------------------------------------
#   Subroutines
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
sub print_usage {
    my $name = basename($0);
    print << "EOF";
Usage: $name [--conf <config_file>] [--macro <output_macro_name>]
             [--outdir <output_directory>] [--help]
EOF
  exit();
}
#-------------------------------------------------------------------------------
sub print_help {
    my $name = basename($0);

    my $title = "$name - Extracts timing information from trace files";
    my $line  = '-' x length($title);
    print << "EOF";

$title
$line

Usage:
------

EOF

    printf("%15s --in <input_file> --out <output_file>\n", $name);
    printf("%15s --tracedir <trace_directory> --nfactors <nb_factors>\n", '');
    printf("%15s [--help]\n", '');

print << "EOF";

The $name script is used after executing a benchmark macro created by
the benchmarker.pl script. It reads the original result file (in the
Tifa::DataDescriptor format) and the trace files generated during benchmarks to
create a more complete result file by adding extra information found in the
trace files. This script is not genering and is tailored to be used if and only
if a TIFA program has been benchmarked since it expects to find specific
information in the trace files.

General parameters/options:
---------------------------

  --in=s
      (mandatory)
      Name of the input result file as generated by the benchmarking script
      created by benchmarker.pl.

  --out=s
      (mandatory)
      Name of the output result file. This file will contain information
      from the original input file and extra information read in the
      trace files generated during the benchmarks.

  --tracedir=s
      (mandatory)
      Name of the directory containing the trace files.

  --nfactors=i
      (optional, meaningless value used if none provided)
      Number of prime factors composing the integers to factor that were used
      in the benchmarks. If not provided, -1 will be used to indicate that the
      real value is unknown.

See also:
---------

  The documentation for the benchmarker.pl script and the Tifa::DataDescriptor
  module.

EOF

  exit();
}
#-------------------------------------------------------------------------------
sub check_options_validity {

    if (! -e $input_file) {
        print("ERROR: Cannot find input file $input_file\n");
        exit();
    }
    if (-e $output_file) {
        print("ERROR: Output file $output_file already exists\n");
        exit();
    }
    if (! -e $trace_dir) {
        print("ERROR: Cannot find directory $trace_dir\n");
        exit();
    }
    if (! -d $trace_dir) {
        print("ERROR: $trace_dir is not a directory\n");
        exit();
    }
}
#-------------------------------------------------------------------------------
sub get_trace_filenames {
    my @files    = <$trace_dir/trace_*.txt>;
    my %returned = ();

    foreach my $file (@files) {
        if ($file =~ /trace_(\d+).txt/) {
            $returned{$file} = int($1);
        } else {
            printf("ERROR: unexpected file name found: $file\n");
            die();
        }
    }
    return (sort {$returned{$a} <=> $returned{$b}} keys %returned);
}
#-------------------------------------------------------------------------------
sub make_new_entry {
    #
    # This subroutine depends on the exact output of the benchmarked programs
    # since it uses regular expressions to extract the timings. It should be
    # changed/adapted whenever the output format of the benchmarked programs
    # changes...
    #
    # This subroutine could probably be isolated in a configuration file, which
    # could be then included right here...
    #
    my $infile        = shift(@_);
    my $old_entry_ref = shift(@_);

    my $start_up_time  = -1.0;
    my $collect_time   = -1.0;
    my $generated_time = -1.0;
    my $selected_time  = -1.0;
    my $ptdiv_time     = -1.0;
    my $bernfac_time   = -1.0;
    my $linalg_time    = -1.0;
    my $dedfac_time    = -1.0;
    my $fwdcycle_time  = -1.0;
    my $invsqrt_time   = -1.0;
    my $revcycle_time  = -1.0;
    my $tdiv_time      = -1.0;
    my $cleaning_time  = -1.0;
    my $total_time     = -1.0;

    my $start_up_str  = quotemeta("computing start up data...");
    my $collect_str   = quotemeta("collecting (X^2 = Y mod N) relations...");
    my $generated_str = quotemeta("residues generated in");
    my $selected_str  = quotemeta("residues selected in");
    my $ptdiv_str     = quotemeta("partial trial divisions...");
    my $bernfac_str   = quotemeta("full decomposition...");
    my $linalg_str    = quotemeta("resolving linear algebra system...");
    my $dedfac_str    = quotemeta("deducing factors...");
    my $cleaning_str  = quotemeta("cleaning...");
    my $fwdcycle_str  = quotemeta("forward cycling to find a proper form...");
    my $invsqrt_str   = quotemeta("computing inverse square root of form...");
    my $revcycle_str  = quotemeta("reverse cycling to find a factor...");
    my $tdiv_str      = 'tdiv:\ .*total';
    my $total_str     = '[^(tdiv)]:\ .*total';

    my $time_regexp      = '\d+\.\d+';
    my $start_up_regexp  = $start_up_str.'\D*('.$time_regexp.")";
    my $collect_regexp   = $collect_str.'\D*('.$time_regexp.")";
    my $generated_regexp = $generated_str.'\D*('.$time_regexp.")";
    my $selected_regexp  = $selected_str.'\D*('.$time_regexp.")";
    my $ptdiv_regexp     = $ptdiv_str.'\D*('.$time_regexp.")";
    my $bernfac_regexp   = $bernfac_str.'\D*('.$time_regexp.")";
    my $linalg_regexp    = $linalg_str.'\D*('.$time_regexp.")";
    my $dedfac_regexp    = $dedfac_str.'\D*('.$time_regexp.")";
    my $cleaning_regexp  = $cleaning_str.'\D*('.$time_regexp.")";
    my $fwdcycle_regexp  = $fwdcycle_str.'\D*('.$time_regexp.")";
    my $invsqrt_regexp   = $invsqrt_str.'\D*('.$time_regexp.")";
    my $revcycle_regexp  = $revcycle_str.'\D*('.$time_regexp.")";
    my $tdiv_regexp      = $tdiv_str.'\D*('.$time_regexp.")";
    my $total_regexp     = $total_str.'\D*('.$time_regexp.")";

    open(IN, "<$infile") or die("Cannot open $infile\n");

    while (my $line = <IN>) {
        #
        # Note that there may be several factoring traces in the output file.
        # Consequently, we increment the proper timing if it has already been
        # initialized...
        #
        if ($line =~ m/$start_up_regexp/) {
            if ($start_up_time > 0) {
                $start_up_time += $1;
            } else {
                $start_up_time = $1;
            }
            next;
        }
        if ($line =~ m/$generated_regexp/) {
            if ($generated_time > 0) {
                $generated_time += $1;
            } else {
                $generated_time = $1;
            }
            next;
        }
        if ($line =~ m/$selected_regexp/) {
            if ($selected_time > 0) {
                $selected_time = $1;
            } else {
                $selected_time = $1;
            }
            next;
        }
        if ($line =~ m/$collect_regexp/) {
            if ($collect_time > 0) {
                $collect_time += $1;
            } else {
                $collect_time = $1;
            }
            next;
        }
        if ($line =~ m/$ptdiv_regexp/) {
            if ($ptdiv_time > 0) {
                $ptdiv_time += $1;
            } else {
                $ptdiv_time = $1;
            }
            next;
        }
        if ($line =~ m/$bernfac_regexp/) {
            if ($bernfac_time > 0) {
                $bernfac_time += $1;
            } else {
                $bernfac_time = $1;
            }
            next;
        }
        if ($line =~ m/$linalg_regexp/) {
            if ($linalg_time > 0) {
                $linalg_time += $1;
            } else {
                $linalg_time = $1;
            }
            next;
        }
        if ($line =~ m/$dedfac_regexp/) {
            if ($dedfac_time > 0) {
                $dedfac_time += $1;
            } else {
                $dedfac_time = $1;
            }
            next;
        }
        if ($line =~ m/$cleaning_regexp/) {
            if ($cleaning_time > 0) {
                $cleaning_time += $1;
            } else {
                $cleaning_time = $1;
            }
            next;
        }
        if ($line =~ m/$fwdcycle_regexp/) {
            if ($fwdcycle_time > 0) {
                $fwdcycle_time += $1;
            } else {
                $fwdcycle_time = $1;
            }
            next;
        }
        if ($line =~ m/$invsqrt_regexp/) {
            if ($invsqrt_time > 0) {
                $invsqrt_time += $1;
            } else {
                $invsqrt_time = $1;
            }
            next;
        }
        if ($line =~ m/$revcycle_regexp/) {
            if ($revcycle_time > 0) {
                $revcycle_time += $1;
            } else {
                $revcycle_time = $1;
            }
            next;
        }
        if ($line =~ m/$tdiv_regexp/) {
            if ($tdiv_time > 0) {
                $tdiv_time += $1;
            } else {
                $tdiv_time = $1;
            }

            next;
        }
        if ($line =~ m/$total_regexp/) {
            #
            # We should NOT increment the total time here but just take the
            # last "total time" line in the output file to account for possible
            # recursive calls (the last "total time" line already takes into
            # account the time taken by these possible recursive calls).
            #
            $total_time = $1;
            next;
        }
    }
    close(IN);

    my %entry = ();

    $entry{"inittime"}    = $start_up_time;
    $entry{"relcoltime"}  = $collect_time;
    $entry{"relgentime"}  = $generated_time;
    $entry{"relseltime"}  = $selected_time;
    $entry{"ptdivtime"}   = $ptdiv_time;
    $entry{"factime"}     = $bernfac_time;
    $entry{"lsrtime"}     = $linalg_time;
    $entry{"deducetime"}  = $dedfac_time;
    $entry{"cleantime"}   = $cleaning_time;
    $entry{"fwdcyctime"}  = $fwdcycle_time;
    $entry{"invsqrttime"} = $invsqrt_time;
    $entry{"revcyctime"}  = $revcycle_time;
    $entry{"tdivtime"}    = $tdiv_time;
    $entry{"totaltime"}   = $total_time;

    $entry{"nfactors"} = $nb_factors;
    $entry{"n"}        = $$old_entry_ref{"argmt"};
    $entry{"nbits"}    = sizeinbase($$old_entry_ref{"argmt"}, 2);

    delete $$old_entry_ref{"argmt"};

    @entry{keys %$old_entry_ref} = values %$old_entry_ref;

    return \%entry;
}
#-------------------------------------------------------------------------------
