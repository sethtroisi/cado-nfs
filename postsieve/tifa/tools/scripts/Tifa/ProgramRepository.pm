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

package Tifa::ProgramRepository;

#-------------------------------------------------------------------------------
#                  A repository of TIFA programs.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : ProgramRepository.pm
# Author        : Jerome Milan
# Created on    : Early Feb 2007
# Last modified : Early Feb 2007
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
# This is a collection of available TIFA::Program's with a very simple
# interface to access them.
#-------------------------------------------------------------------------------

use 5.006002;
use strict;
use warnings;
use Class::Struct;
use File::Basename;
use Carp;

use Tifa::Program;

require Exporter;

our @ISA        = qw(Exporter);
our @EXPORT_OK  = qw();
our $VERSION    = '0.1.0';

my %algo_to_program = ();


my @param_names   = ();
my @param_descrs  = ();
my @param_types   = ();
my $cmdline       = '';
my $help          = '';

#-------------------------------------------------------------------------------
#                           CFRAC algorithm
#-------------------------------------------------------------------------------
my $cfrac_program = new Tifa::Program();

$algo_to_program{"cfrac"} = $cfrac_program;

$cfrac_program->set_algo("cfrac");
$cfrac_program->set_descr("The Continued FRACtion algorithm");
$cfrac_program->set_exe("cfrac_factors");

@param_names = (
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_nb",
    "nrelations",
    "lsr_method",
    "use_large_primes",
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes in factor base",
    "Number of primes to trial divide by the residues",
    "Number of relations",
    "Method of linear system resolution",
    "Whether to use the large prime variation or not",
    "Number of primes to trial divide the number to factor",
);
@param_types = (
    "int",
    "int",
    "int",
    "int",
    "switch",
    "int",
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_nb",
    "nrelations",
    "lsr_method",
    "use_large_primes",
    "nprimes_tdiv",
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$cfrac_program->add_mode("no_defaults",
                         "Default mode: specify all parameter values",
                         \@param_names, \@param_descrs, \@param_types,
                         $cmdline);

$cfrac_program->set_default_mode("no_defaults");
$cfrac_program->set_mode("no_defaults");

@param_names  = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
);
@param_types  = (
    "int",
);
$cmdline =  "exe nprimes_tdiv";
#
# use_defaults mode: Let CFRAC use the precomputed optimal parameter values
#                    depending on the size of the number to factor
#
$cfrac_program->add_mode("use_defaults",
                         "'Best' mode: let CFRAC use optimal values",
                         \@param_names, \@param_descrs, \@param_types,
                         $cmdline);
$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to cfrac _and_ the exe
      parameter is not defined or left empty.
      Default: cfrac_factors

  --nprimes_in_factor_base=i
      (mandatory in 'no_defaults' mode)
      Number of primes of the factor base.

  --nprimes_tdiv_smooth_nb=i
      (mandatory in 'no_defaults' mode)
      Number of primes to consider in the trial division of the
      smooth residues.

  --nrelations=i
      (mandatory in 'no_defaults' mode)
      Number of relations to find during the linear resolution
      phase of the CFRAC algorithm.

  --lsr_method=i
      (mandatory in 'no_defaults' mode)
      Linear system resolution method to use.
      Accepted values:
          0 : Smart gaussian elimination

  --use_large_primes
      (mandatory in 'no_defaults' mode)
      Use the large prime variation.

  --mode=s
      (optional, default value used if none provided)
      Sets the program mode to use.

      use_defaults: In this mode, the CFRAC program will determine the values
                    of the aforementioned parameters. These values are choosen
                    according to the size of the number to factor and are more
                    or less optimal if this size is between 60 and 200 bits.

      no_defaults:  In this mode, all parameters should be specified.

      Default value: no_defaults

  --nprimes_tdiv=i
      (mandatory in all modes)
      Number of primes to use to trial divide the integer to factor.

EOS

$cfrac_program->set_help($help);

#-------------------------------------------------------------------------------
#                       Fermat's algorithm (McKee's speedup)
#-------------------------------------------------------------------------------
my $fermat_program = new Tifa::Program();

$algo_to_program{"cfrac"} = $fermat_program;

$fermat_program->set_algo("fermat");
$fermat_program->set_descr("Fermat's algorithm (McKee's speedup)");
$fermat_program->set_exe("fermat_factors");

@param_names = (
    "nprimes_tdiv"
);
@param_descrs = (
    "Number of primes to trial divide the number to factor"
);
@param_types = (
    "int"
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_tdiv"
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$fermat_program->add_mode("no_defaults",
                          "Default mode: specify all parameter values",
                          \@param_names, \@param_descrs, \@param_types,
                          $cmdline);

$fermat_program->set_default_mode("no_defaults");
$fermat_program->set_mode("no_defaults");

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to fermat _and_ the exe
      parameter is not defined or left empty.
      Default: fermat_factors

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.

EOS

$fermat_program->set_help($help);

#-------------------------------------------------------------------------------
#                            QS algorithm
#-------------------------------------------------------------------------------
my $qs_program = new Tifa::Program();

$algo_to_program{"qs"} = $qs_program;

$qs_program->set_algo("qs");
$qs_program->set_descr("The Quadratic Sieve algorithm");
$qs_program->set_exe("qs_factors");

@param_names = (
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_nb",
    "nrelations",
    "lsr_method",
    "use_large_primes",
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes in factor base",
    "Number of primes to trial divide by the residues",
    "Number of relations",
    "Method of linear system resolution",
    "Whether to use the large prime variation or not",
    "Number of primes to trial divide the number to factor",
);
@param_types = (
    "int",
    "int",
    "int",
    "int",
    "switch",
    "int",
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_nb",
    "nrelations",
    "lsr_method",
    "use_large_primes",
    "nprimes_tdiv",
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$qs_program->add_mode("no_defaults",
                      "Default mode: specify all parameter values",
                      \@param_names, \@param_descrs, \@param_types,
                      $cmdline);

$qs_program->set_default_mode("no_defaults");
$qs_program->set_mode("no_defaults");

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to qs _and_ the exe
      parameter is not defined or left empty.
      Default: qs_factors

  --nprimes_in_factor_base=i
      (mandatory)
      Number of primes of the factor base.

  --nprimes_tdiv_smooth_nb=i
      (mandatory)
      Number of primes to consider in the trial division of the
      smooth residues.

  --nrelations=i
      (mandatory)
      Number of relations to find during the linear resolution
      phase of the QS algorithms.

  --lsr_method=i
      (mandatory)
      Linear system resolution method to use.
      Accepted values:
          0 : Smart gaussian elimination

  --use_large_primes
      (optional)
      Use the large prime variation

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.
EOS

$qs_program->set_help($help);
#-------------------------------------------------------------------------------
#                           SIQS algorithm
#-------------------------------------------------------------------------------
my $siqs_program = new Tifa::Program();

$algo_to_program{"siqs"} = $siqs_program;

$siqs_program->set_algo("siqs");
$siqs_program->set_descr("The Self-Initializing Quadratic Sieve algorithm");
$siqs_program->set_exe("siqs_factors");

@param_names = (
    "sieve_half_width",
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_nb",
    "nrelations",
    "lsr_method",
    "use_large_primes",
    "nprimes_tdiv",
);
@param_descrs = (
    "Half width of the sieving interval",
    "Number of primes in factor base",
    "Number of primes to trial divide by the residues",
    "Number of relations",
    "Method of linear system resolution",
    "Whether to use the large prime variation or not",
    "Number of primes to trial divide the number to factor",
);
@param_types = (
    "int",
    "int",
    "int",
    "int",
    "int",
    "switch",
    "int",
);
$cmdline = join(
    ' ',
    "exe",
    "sieve_half_width",
    "nprimes_in_factor_base",
    "nprimes_tdiv_smooth_nb",
    "nrelations",
    "lsr_method",
    "use_large_primes",
    "nprimes_tdiv",
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$siqs_program->add_mode("no_defaults",
                        "Default mode: specify all parameter values",
                        \@param_names, \@param_descrs, \@param_types,
                        $cmdline);

$siqs_program->set_default_mode("no_defaults");
$siqs_program->set_mode("no_defaults");

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to siqs _and_ the exe
      parameter is not defined or left empty.
      Default: siqs_factors

  --sieve_half_width=i
      (mandatory)
      Half width of the sieving interval. The complete sieving interval is
      thus given by [-sieve_half_width, +sieve_half_width].

  --nprimes_in_factor_base=i
      (mandatory)
      Number of primes of the factor base.

  --nprimes_tdiv_smooth_nb=i
      (mandatory)
      Number of primes to consider in the trial division of the
      smooth residues.

  --nrelations=i
      (mandatory)
      Number of relations to find during the linear resolution
      phase of the SIQS algorithms.

  --lsr_method=i
      (mandatory)
      Linear system resolution method to use.
      Accepted values:
          0 : Smart gaussian elimination

  --use_large_primes
      (optional)
      Use the large prime variation

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.
EOS
$siqs_program->set_help($help);
#-------------------------------------------------------------------------------
#                          SQUFOF algorithm
#-------------------------------------------------------------------------------
my $squfof_program = new Tifa::Program();

$algo_to_program{"squfof"} = $squfof_program;

$squfof_program->set_algo("squfof");
$squfof_program->set_descr("The SQUare FOrm Factorization algorithm");
$squfof_program->set_exe("squfof_factors");

@param_names = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes to trial divide the number to factor",
);
@param_types = (
    "int",
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_tdiv",
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$squfof_program->add_mode("no_defaults",
                          "Default mode: specify all parameter values",
                          \@param_names, \@param_descrs, \@param_types,
                          $cmdline);

$squfof_program->set_default_mode("no_defaults");
$squfof_program->set_mode("no_defaults");

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to squfof _and_ the exe
      parameter is not defined or left empty.
      Default: squfof_factors

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.

EOS
$squfof_program->set_help($help);
#-------------------------------------------------------------------------------
#                          Trial division algorithm
#-------------------------------------------------------------------------------
my $tdiv_program = new Tifa::Program();

$algo_to_program{"tdiv"} = $tdiv_program;

$tdiv_program->set_algo("tdiv");
$tdiv_program->set_descr("The naive trial division algorithm");
$tdiv_program->set_exe("tdiv_factors");

@param_names = (
    "nprimes_tdiv",
);
@param_descrs = (
    "Number of primes used to trial divide the number to factor",
);
@param_types = (
    "int",
);

$cmdline = join(
    ' ',
    "exe",
    "nprimes_tdiv",
);
#
# Default mode: We should specify all of the parameter values on
#               the command line
#
$tdiv_program->add_mode("no_defaults",
                        "Default mode: specify all parameter values",
                        \@param_names, \@param_descrs, \@param_types,
                        $cmdline);

$tdiv_program->set_default_mode("no_defaults");
$tdiv_program->set_mode("no_defaults");

$help = <<EOS;
Parameter(s):

  --exe=s
      (optional, see general help)
      Name of the executable command line program. Note that the default value
      is only used if the algo parameter is set to tdiv _and_ the exe
      parameter is not defined or left empty.
      Default: tdiv_factors

  --nprimes_tdiv=i
      (mandatory)
      Number of primes to use to trial divide the integer to factor.
EOS
$tdiv_program->set_help($help);
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Subroutines
#-------------------------------------------------------------------------------
sub get_algo_to_program_hash {
    return %algo_to_program;
}
#-------------------------------------------------------------------------------
sub get_algo_list {
    return keys %algo_to_program;
}
#-------------------------------------------------------------------------------
sub get_program_list {
    return values %algo_to_program;
}
#-------------------------------------------------------------------------------
sub get_all_param_names {
    my @all_params = ();
    foreach my $algo (keys %algo_to_program) {
        push(@all_params, $algo_to_program{$algo}->get_all_param_names());
    }
    my %hash = map {$_, 1} @all_params; # get rid of duplicate values

    return keys %hash;
}
#-------------------------------------------------------------------------------
sub get_all_param_to_descrs_hash {
    my %hash = ();
    foreach my $algo (keys %algo_to_program) {
        my %subhash = $algo_to_program{$algo}->get_all_param_to_descrs_hash();
        %hash = (%hash, %subhash);
    }
    return %hash;
}
#-------------------------------------------------------------------------------
sub get_all_param_to_types_hash {
    my %hash = ();
    foreach my $algo (keys %algo_to_program) {
        my %subhash = $algo_to_program{$algo}->get_all_param_to_types_hash();
        %hash = (%hash, %subhash);
    }
    return %hash;
}
#-------------------------------------------------------------------------------
sub get_getopt_strings {
    my $self = shift;

    my %hash = get_all_param_to_types_hash();

    my @getopt_strings = ();
    my $tmp = "";

    foreach my $param (keys %hash) {
        $tmp = $param;
        if ($hash{$param} ne "switch") {
            if ($hash{$param} eq "string") {
                $tmp .= '=s';
                push(@getopt_strings, $tmp);
            }
            if ($hash{$param} eq "int") {
                $tmp .= '=i';
                push(@getopt_strings, $tmp);
            }
            if ($hash{$param} eq "extint") {
                $tmp .= '=o';
                push(@getopt_strings, $tmp);
            }
            if ($hash{$param} eq "float") {
                $tmp .= '=f';
                push(@getopt_strings, $tmp);
            }
        } else {
            push(@getopt_strings, $tmp);
        }
    }
    push(@getopt_strings, "algo=s");

    return @getopt_strings;
}
#-------------------------------------------------------------------------------
sub get_program_from_name {
    my $name     = basename(shift);

    foreach my $prog (values %algo_to_program) {
        if ( $name eq basename($prog->get_exe()) ) {
            return $prog;
        }
    }
    carp("ERROR: ProgramRepository: $name is not a known program");
    return undef;
}
#-------------------------------------------------------------------------------

1;

__END__

#-------------------------------------------------------------------------------
# Stub documentation for this module.
#-------------------------------------------------------------------------------

=head1 NAME

Tifa::ProgramRepository - A repository of available Tifa::Program's

=head1 SYNOPSIS

 use Tifa::ProgramRepository;

 %algo_to_program = Tifa::ProgramRepository::get_algo_to_program_hash();

=head1 REQUIRE

Perl 5.006002, Carp, Class::Struct, Exporter and Tifa::Program.

=head1 SUMMARY

The Tifa::ProgramRepository module acts as a repository of all available
Tifa::Program objects.

=head1 DESCRIPTION

The Tifa::ProgramRepository module acts as a repository of all available
(factoring) Tifa::Program's. The listed Tifa::Program objects can then be used
in other scripts or modules.

=head2 Available methods

    get_algo_to_program_hash()
    get_algo_list()
    get_program_list()
    get_program_from_name($progname)
    get_all_param_names()
    get_all_param_to_descrs_hash()
    get_all_param_to_types_hash()

=head2 Methods description

    get_algo_to_program_hash()
        Returns a hashtable mapping algorithm names to the Tifa::Program
        objects implementing it.

    get_algo_list()
        Returns an array listing the names of the implemented algorithm.

    get_program_list()
        Returns an array listing the Tifa::Program objects available.

    get_program_from_name($progname)
        Returns a reference to the Tifa::Program object whose default
        program name is given by $progname (or the base name of $progname
        if $progname is a path).

    get_all_param_names()
        Returns an array listing all the parameter names from all
        programs.

    get_all_param_to_descrs_hash()
        Returns a hashtable mapping all the parameter names (from all
        programs and regardless of their current modes), to their
        descriptions.

    get_all_param_to_types_hash()
        Returns a hashtable mapping all the parameter types (from all
        programs and regardless of their current modes), to their
        descriptions.

=head1 EXAMPLE

This following code snippet gives an example of how the Tifa::ProgramRepository
module can be used.

 %algo_to_program = Tifa::ProgramRepository::get_algo_to_program_hash();

 $program = $algo_to_program{"squfof"};

 %param_vals = {
     "exe"       => "./squfof_factors",
     "nb_tdiv"   => 128,
 };

 $program->execute(\%param_vals, 816379, "> trace.txt");

=head1 EXPORT

No functions are exported from this package by default.

=head1 SEE ALSO

The Tifa::Program module's man page.

=head1 AUTHOR

Jerome Milan, E<lt>milanj@lix.polytechnique.frE<gt>

=head1 VERSION

0.1.0

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2006, 2007 INRIA (French National Institute for Research in
Computer Science and Control)

This module is part of the TIFA library (Tools for Integer FActorization).

The TIFA library is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

The TIFA library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with this library; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA

=cut
