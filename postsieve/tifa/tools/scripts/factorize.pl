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
#          A simple wrapper for the TIFA library's factoring programs
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : factorize.pl
# Author        : Jerome Milan
# Created on    : Mon Mar 20 2006
# Last modified : Mon Jan 29 2007
#
# Version : 0.2.0
# License : GNU Lesser General Public License (LGPL)
#           Copyright (C) 2006, 2007 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
#  0.2.0: Mon Jan 29 2007 by JM
#      - Lots of refactoring to make addition of new algorithm easier and
#        cleaner. Now uses the Tifa::ProgramRepository mechanism.
#
#  0.1.4: Tue Oct  3 2006 by JM
#      - Renamed as factorize.pl
#      - Added qs_factors and siqs_factors support and updated help message
#
#  0.1.3: Thu Jul 27 2006 by JM
#      - Now uses the Tifa::SimpleConfigReader to read the configuration file
#        but this does not change the format of the configuration file (for
#        this script anyway...)
#
#  0.1.2: Thu Jul 13 2006 by JM
#      - Added option to read a list of numbers to factor from a file
#      - Added option to let the cfrac program choose the default parameters
#
#  0.1.1: Thu May 18 2006 by JM
#      - No more nb_step_cfrac parameter
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# The factorize.pl script is a very simple wrapper for the cfrac_factors,
# qs_factors, siqs_factors and squfof_factors programs. The script merely adds
# some syntaxic sugar for a more user friendly use. It also adds the possibility
# to use a file containing a list of numbers as input for the aforementionned
# programs.
#-------------------------------------------------------------------------------

use strict;

use Getopt::Long;
use File::Basename;
use Class::Struct;

use Tifa::SimpleConfigReader;
use Tifa::ProgramRepository;

#
# Fetch informations on available algorithms and programs...
#
my @getopt_strings = ();
my @param_names    = ();
my @programs       = ();

my %param_values = ();
my %param_descrs = ();
my %param_types  = ();

my %algo_to_program = ();

my $program;

@param_names     = Tifa::ProgramRepository::get_all_param_names();
@getopt_strings  = Tifa::ProgramRepository::get_getopt_strings();

%algo_to_program = Tifa::ProgramRepository::get_algo_to_program_hash();
%param_descrs    = Tifa::ProgramRepository::get_all_param_to_descrs_hash();
%param_types     = Tifa::ProgramRepository::get_all_param_to_types_hash();

#
# Options/arguments exclusively passed via the command line...
#
push(@param_names, "conf");
push(@param_names, "help");
push(@param_names, "numbers");
push(@param_names, "outdir");

push(@getopt_strings, "conf=s");
push(@getopt_strings, "help:s");
push(@getopt_strings, "numbers=s");
push(@getopt_strings, "outdir=s");
push(@getopt_strings, "cmd:s");

$param_descrs{"conf"}    = "Name of the configuration file";
$param_descrs{"help"}    = "Show complete help message";
$param_descrs{"numbers"} = "Name of input number file";
$param_descrs{"outdir"}  = "Name of ouput directory (to hold trace files)";
$param_descrs{"cmd"}     = "Show cmd";

$param_types{"conf"}    = "string";
$param_types{"help"}    = "switch";
$param_types{"numbers"} = "string";
$param_types{"outdir"}  = "string";
$param_types{"cmd"}     = "switch";

$param_values{"conf"}    = undef;
$param_values{"numbers"} = undef;
$param_values{"outdir"}  = undef;

GetOptions(\%param_values, @getopt_strings);

if (defined $param_values{"help"} ) {
    print_help($param_values{"help"});
    exit();
}

check_args_validity();
check_options_validity();

#
# Read configuration file if provided
#
if ($param_values{"conf"}) {
    my $confreader = new Tifa::SimpleConfigReader();

    $confreader->accept_unknown_param();
    $confreader->read_config_file($param_values{"conf"});

    my %param_hash = $confreader->get_parameter_hash();
    #
    # Set the values of the parameters to the one just read...
    #
    foreach my $name (keys %param_hash) {
        #
        # Give precedence to the options passed on the command line
        #
        if (!defined $param_values{$name} && defined $param_hash{$name}) {
            $param_values{$name} = $param_hash{$name};
        }
    }
}

#
# Infer algorithm used if the algo field is not a known algorithm
#
if (defined $algo_to_program{$param_values{"algo"}}) {
    $program = $algo_to_program{$param_values{"algo"}};
} else {
    if ((!defined $param_values{"algo"}) || ($param_values{"algo"} eq "")) {
        #
        # Infer the algorithm to use from the program names
        #
        foreach my $prog (values %algo_to_program) {
            if (basename($param_values{"exe"}) eq $prog->get_exe()) {
                $program = $prog;
                last;
            }
        }
        if (!defined $program) {
            print("ERROR: Could not infer which Tifa::Program to use!\n");
            exit();
        }
    } else {
        print("ERROR: ", $param_values{"algo"}, " is an unknown algorithm!\n");
        exit();
    }
}

#
# Set the mode of the Tifa::Program to use. Revert to the default mode if
# no mode or invalid mode provided.
#
$program->set_mode($program->get_default_mode());

if (defined $param_values{"mode"}) {
    #
    # This reverts to the defaults mode if the given mode name is unknown
    #
    $program->set_mode($param_values{"mode"});
}

#
# Set the name of the executable of the Tifa::Program to use. Revert to the
# default name if no name provided.
#
if (defined $param_values{"exe"}) {
    $program->set_exe($param_values{"exe"});
}

my @numbers = ();
if (defined $param_values{"numbers"}) {
    #
    # Reads the numbers from the file $number_file and stores then in the
    # array @numbers...
    #
    &read_numbers_from_file($param_values{"numbers"}, \@numbers);
    #
    # Create the output directory to hold the trace files...
    #
    if (! -e $param_values{"outdir"}) {
        mkdir($param_values{"outdir"}, 0700)
            or die("ERROR: Cannot create directory ", $param_values{"outdir"},
                   "!\n");
    }
    my $nb_numbers = @numbers;         # Number of numbers to factor
    my $size_field = length(@numbers); # Size of the counter field in the trace
                                       # file's name

    #
    # Now, calls the factoring program for each number to factor
    #
    for (my $i = 0; $i < $nb_numbers; $i++) {
        #
        # Save the ouput in a file named trace_XXXXX.txt...
        #
        my $trace_file = sprintf("trace_%0$size_field"."d", $i).".txt";
        my $postproc   = " > ".$param_values{"outdir"}."/$trace_file";

        if (defined $param_values{"cmd"} ) {
            print $program->make_cmd(\%param_values, $numbers[$i], $postproc);
            print "\n";
        } else {
            $program->execute(\%param_values, $numbers[$i], $postproc);
        }
    }
} else {
    #
    # There is only one number to factor: the one provided as argument
    # on the command line. The ouput is just printed on stdout. No trace
    # file is kept.
    #
    if (defined $param_values{"cmd"} ) {
        print $program->make_cmd(\%param_values, $ARGV[0])."\n";
    } else {
        $program->execute(\%param_values, $ARGV[0]);
    }
}
#
# End of script.
#

#-------------------------------------------------------------------------------
#   Subroutines
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
sub print_usage {
    my $help = shift;
    my $name = basename($0);

    printf "Usage:\n";
    printf "------\n\n";
    printf "%15s [parameters] [options] <number_to_factor>\n", $name;
    printf "\nor:\n\n";
    printf "%15s [parameters] [options] --numbers <file>", $name;
    printf " --outdir <directory>\n\n";

    print  "Type $name --help for more information.\n" if (defined $help);
}
#-------------------------------------------------------------------------------
sub print_help {
    my $help_algo = shift;
    my $name      = basename($0);

    if ($help_algo eq "config") {
        print_config_help();
        return;
    }

    foreach my $algo (keys %algo_to_program) {
        if ($help_algo eq $algo) {
            print_algo_help($help_algo);
            return;
        }
    }
    if ($help_algo ne "") {
        print("ERROR: print_help: unknown algorithm $help_algo\n");
        return;
    }

    my $title = "$name - A Perl wrapper for TIFA's factorization programs";
    my $line  = '-' x length($title);
    print << "EOF";

$title
$line

EOF

    print_usage;

    print << "EOF";
The $name script wraps the various TIFA factorization programs and offers
a user-friendly command line interface which can be easily extended to support
new factoring program, as far as the new programs are wrapped in a
Tifa::Program object and registered in the Tifa::ProgramRepository module.

This script can be used to factor a single composite integer given as an
argument on the command line or to factor a whole set of composite numbers
given in a text file.

General parameters/options:
---------------------------

  --algo=s?
      (optional, but --exe should then be defined)
      Name of the factorization algorithm to use. Should be one of the
      following value:
EOF

      foreach my $algo (sort {$a cmp $b} keys %algo_to_program) {
          printf "      %8s  (%s)\n", $algo,
                                      $algo_to_program{$algo}->get_descr();
      }
      my $inferinfo = "if empty, infers algorithm from program's default name";
      printf "      %8s  ($inferinfo)\n", "<empty>";

print << "EOF";

  --exe=s
        (optional, but --algo should then be defined)
        Name of the executable command line program. Note that the default value
        is only used if the algo parameter is set to a known algorithm name
        _and_ the exe parameter is not defined or left empty.
        Default: vary according to algorithm used.
                 Type $name --help <algo_name> for more information.

  --cmd
      Print the command line used to call the underlying C program and exit.

  --conf=s
      (optional, on command-line only)
      Name of the configuration file listing the values of the factorization
      algorithms' parameters written as: "<option_name> = <value>".
      Type $name --help config for information about the configuration
      file format.

  --numbers=s
      (optional, on command-line only)
      Do not read the number to factorize on the command line but read a
      set of numbers given in a text file. Numbers can be separated by any
      kind of non digit characters. Comments beginning by # and blank lines
      are allowed in the file. If this option is used the number passed as
      argument -if any- will be discarded.

  --outdir=s
      (mandatory if --numbers is used, on command-line only)
      Name of the directory to contain the various execution traces. Only
      used if the --numbers option is used. If the number to factor is
      directly given on the command line, no trace will be saved.

  --help=s?
      (optional, on command-line only)
      Without argument, prints this help message.
      With an algorithm name, prints the help message related to the use of
      this script with this particular algorithm.

Algorithm-specific parameters/options:
--------------------------------------

Type $name --help <algo_name> (where <algo_name> can take one of the
value mentionned in the above description of the --algo option) for more
information about algorithm specific parameters and options.

EOF
}
#-------------------------------------------------------------------------------
sub print_config_help {
    my $name = basename($0);
    my $title = "$name - A Perl wrapper for TIFA's factorization programs";
    my $line  = '-' x length($title);
    print << "EOF";

$title
$line

EOF

    print <<"EOF";
Most of $name\'s parameters and options can be given in a configuration
file whose syntax should conform to the one defined by the
Tifa::SimpleConfigReader module. (Type 'man Tifa::SimpleConfigReader' for
more information.)

These parameters and options should be written as: "<option_name> = <value>".
Note that in the configuration file, the option names should be written
without the introducing double dashes (--) that one would normally use on
a command line. For example, to set the name of the algorithm to 'squfof',
one can type:

    '--algo=squfof' on the command line
    'algo=squfof'   in a configuration file

Some options cannot be passed in a configuration file. See the general help
for more information.

Comments beginning by # and blank lines are allowed in the configuration file.
However, directives should not span multiple lines unless line breaks are
explicitely written as ' \\' (a space followed by a backslash).

See the provided 'factorize.conf' configuration file for a detailed example.

Type $name --help for more information.";
EOF
}
#-------------------------------------------------------------------------------
sub print_algo_help {
    my $algo     = shift;
    my $help_msg = $algo_to_program{$algo}->get_help();
    my $name     = basename($0);

    my $title = "$name - A Perl wrapper for TIFA's factorization programs";
    my $line  = '-' x length($title);
        print << "EOF";

$title
$line

EOF

    print_usage;

    my $descr = $algo_to_program{$algo}->get_descr();
    my $len   = (length $descr) % 80;

    print  "$descr\n";
    print  '-' x $len, "\n\n";

    print << "EOF"
This algorithm is selected by using the --algo=$algo option. Mandatory
parameters (if any) should be given on the command line (--option-name) or
in a configuration file (option-name=<value>).

For more information on the configuration file format, type
$name --help config.

For information on the non-mandatory options, type $name --help.

$help_msg
Type $name --help for more information.
EOF
}
#-------------------------------------------------------------------------------
sub check_args_validity {
    if (defined $param_values{"numbers"}) {
        return; # Numbers read from file: ignore any argument
    }
    if ($#ARGV+1 == 0) {
        print "ERROR: no argument provided\n\n";
        print_usage(1);
        exit();
    }
    if ($#ARGV+1 != 1) {
        print "ERROR: too many arguments\n\n";
        print_usage(1);
        exit();
    }
    #
    # Remove any leading or trailing whitespaces...
    #
    $ARGV[0] =~ s/^\s+//;
    $ARGV[0] =~ s/\s+$//;
    #
    # Check if $ARGV[0] is a positive integer...
    #
    if (! is_natural($ARGV[0])) {
      print "ERROR: $ARGV[0] is not a valid positive integer\n\n";
      print_usage(1);
      exit();
    }
}
#-------------------------------------------------------------------------------
sub check_options_validity {
    if (defined  $param_values{"numbers"}) {
        if (! -e $param_values{"numbers"}) {
            print "ERROR: Cannot find file ", $param_values{"numbers"}, "\n\n";
            print_usage(1);
            exit();
        }
        if (! defined $param_values{"outdir"}) {
            print "ERROR: No output directory specified\n\n";
            print_usage(1);
            exit();
        }
        if (-e $param_values{"outdir"}) {
            print "ERROR: Output directory ", $param_values{"outdir"};
            print " already exists\n\n";
            print_usage(1);
            exit();
        }
    }
}
#-------------------------------------------------------------------------------
sub read_numbers_from_file {
    #
    # Reads a file containing a list of number to factor and stores them
    # in an array. This function is quite lenient and will happily process
    # whatever you throw at it... The numbers just need to be separated
    # by non digit characters... Comments beginning by # are allowed...
    #
    my $filename  = shift(@_);
    my $arrayref  = shift(@_);

    open(IN, "<$filename") or die("ERROR: Cannot read file $filename\n");

    while (my $line = <IN>) {
        next if ($line =~ m/^(\s)*$/);     # Skip blank lines...
        next if ($line =~ m/^(\s)*\#.*$/); # Skip comments...

        chomp($line);      # Suppress trailing newline if any
        $line =~ s/#.*//;  # Suppress trailing comments
        $line =~ s/^\D*//; # Suppress leading non digit characters
        $line =~ s/\D*$//; # Suppress trailing non digit characters

        push(@$arrayref, split(/\D+/, $line));
    }
    close(IN);
}
#-------------------------------------------------------------------------------
sub is_natural {
    my $string = shift(@_);

    if ($string =~ /^[+]?\d+$/) {
        return 1;
    } else {
        return 0;
    }
}
#-------------------------------------------------------------------------------
