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

package Tifa::Program;

#-------------------------------------------------------------------------------
#                  An abstraction of a TIFA program.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : Program.pm
# Author        : Jerome Milan
# Created on    : Wed Jan 31 2007
# Last modified : Wed Jan 31 2007
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
# This is an abstraction of what a program is from a user's point of view.
# The idea is to set-up a generic framework for program launchers and
# benchmarks to falicitate the addition of new programs to be launched and/or
# benchmarked.
#-------------------------------------------------------------------------------

use 5.006002;
use strict;
use warnings;
use Class::Struct;
use Carp;

require Exporter;

our @ISA        = qw(Exporter);
our @EXPORT_OK  = qw();
our $VERSION    = '0.1.0';

#-------------------------------------------------------------------------------
struct ProgramMode => {
    descr         => '$',  # brief description of mode
    param_names   => '@',  # list of program's parameter names
    param_descrs  => '@',  # list of program's parameter descriptions
    param_types   => '@',  # list of program's parameter types
    cmdline       => '$',  # command line of program in mode
};
#-------------------------------------------------------------------------------
sub new {
    my $type = shift(@_);
    my $self = {};

    $self->{algo}    = '';
    $self->{descr}   = '';
    $self->{exe}     = '';
    $self->{help}    = '';

    $self->{cur_mode}         = '';
    $self->{cur_param_names}  = ();
    $self->{cur_param_descrs} = ();
    $self->{cur_param_types}  = ();
    $self->{cur_cmdline}      = '';

    $self->{default_mode}  = '';
    $self->{modes}         = ();

    #
    # Array listing all of the parameter names
    #
    $self->{all_param_names} = ();
    #
    # Hashtable mapping the parameter names to their descriptions
    #
    $self->{all_param_to_descrs_hash} = ();
    #
    # Hashtable mapping the parameter names to their types
    #
    $self->{all_param_to_types_hash} = ();
    #
    # Hashtable mapping the parameter names to nothing. Used to quickly
    # check if a parameter name if already known
    #
    $self->{all_param_names_hash} = ();

    bless $self, $type;

    return $self;
}
#-------------------------------------------------------------------------------
sub set_help {
    my $self = shift;
    my $help = shift;
    $self->{help} = $help;
}
#-------------------------------------------------------------------------------
sub get_help {
    my $self = shift;
    return $self->{help};
}
#-------------------------------------------------------------------------------
sub set_algo {
    my $self = shift;
    my $algo = shift;
    $self->{algo} = $algo;
}
#-------------------------------------------------------------------------------
sub get_algo {
    my $self = shift;
    return $self->{algo};
}
#-------------------------------------------------------------------------------
sub set_descr {
    my $self  = shift;
    my $descr = shift;
    $self->{descr} = $descr;
}
#-------------------------------------------------------------------------------
sub get_descr {
    my $self = shift;
    return $self->{descr};
}
#-------------------------------------------------------------------------------
sub set_exe {
    my $self = shift;
    my $exe = shift;
    $self->{exe} = $exe;
}
#-------------------------------------------------------------------------------
sub get_exe {
    my $self = shift;
    return $self->{exe};
}
#-------------------------------------------------------------------------------
sub set_default_mode {
    my $self = shift;
    my $mode = shift;
    if (exists ${$self->{modes}}{$mode}) {
        $self->{default_mode} = $mode;
    } else {
        croak("ERROR: Program::set_default_mode: unknown mode $mode\n");
    }
}
#-------------------------------------------------------------------------------
sub get_default_mode {
    my $self = shift;
    return $self->{default_mode};
}
#-------------------------------------------------------------------------------
sub add_mode {
    my $self    = shift;
    my $name    = shift;
    my $descr   = shift;

    my $param_names_ref  = shift;
    my $param_descrs_ref = shift;
    my $param_types_ref  = shift;
    my $cmdline          = shift;

    my @param_names  = @$param_names_ref;
    my @param_descrs = @$param_descrs_ref;
    my @param_types  = @$param_types_ref;

    if (exists $self->{modes}{$name}) {
        croak("ERROR: Program::add_mode: mode $name already exists!");
    };

    #
    # Creates and adds a new mode
    #
    $self->{modes}{$name} = ProgramMode->new();

    $self->{modes}{$name}->descr($descr);
    $self->{modes}{$name}->param_names([@param_names]);
    $self->{modes}{$name}->param_descrs([@param_descrs]);
    $self->{modes}{$name}->param_types([@param_types]);
    $self->{modes}{$name}->cmdline($cmdline);

    #
    # Update program's data members
    #
    foreach my $i (0 .. $#param_names) {
        my $param = $param_names[$i];
        #
        # Provide an easy access to all of the parameters needed, regardless
        # of the current mode.
        #
        if ( exists ${$self->{all_param_names_hash}}{$param} ) {
            ${$self->{all_params_names_hash}}{$param}++;
        } else {
            ${$self->{all_params_names_hash}}{$param} = 1;
            push(@{$self->{all_param_names}}, $param);
        }
        ${$self->{all_param_to_descrs_hash}}{$param} = $param_descrs[$i];
        ${$self->{all_param_to_types_hash}}{$param}  = $param_types[$i];
    }
}
#-------------------------------------------------------------------------------
sub get_mode {
    my $self = shift;
    return $self->{cur_mode};
}
#-------------------------------------------------------------------------------
sub set_mode {
    my $self = shift;
    my $key  = shift;

    if (exists ${$self->{modes}}{$key}) {
        $self->{cur_mode} = $key;
    } else {
        #
        # Reverts to the defaults mode if the given mode name is unknown
        #
        $self->{cur_mode} = $self->{default_mode};
    }
    $self->{cur_param_names}
        = ${$self->{modes}}{$self->{cur_mode}}->param_names;
    $self->{cur_param_descrs}
        = ${$self->{modes}}{$self->{cur_mode}}->param_descrs;
    $self->{cur_param_types}
        = ${$self->{modes}}{$self->{cur_mode}}->param_types;
    $self->{cur_cmdline}
        = ${$self->{modes}}{$self->{cur_mode}}->cmdline;
}
#-------------------------------------------------------------------------------
sub get_all_param_names {
    my $self = shift;
    return (@{$self->{all_param_names}});
}
#-------------------------------------------------------------------------------
sub get_all_param_to_descrs_hash {
    my $self = shift;
    return %{$self->{all_param_to_descrs_hash}};
}
#-------------------------------------------------------------------------------
sub get_all_param_to_types_hash {
    my $self = shift;
    return %{$self->{all_param_to_types_hash}};
}
#-------------------------------------------------------------------------------
sub get_param_names {
    my $self = shift;
    return (@{$self->{modes}{$self->{cur_mode}}->param_names});
}
#-------------------------------------------------------------------------------
sub get_param_descrs {
    my $self = shift;
    return (@{$self->{modes}{$self->{cur_mode}}->param_descrs});
}
#-------------------------------------------------------------------------------
sub get_param_types {
    my $self = shift;
    return (@{$self->{modes}{$self->{cur_mode}}->param_types});
}
#-------------------------------------------------------------------------------
sub get_cmdline {
    my $self = shift;
    return ($self->{modes}{$self->{cur_mode}}->cmdline);
}
#-------------------------------------------------------------------------------
sub get_mode_descr {
    my $self = shift;
    return (@{$self->{modes}{$self->{cur_mode}}->descr});
}
#-------------------------------------------------------------------------------
sub get_all_modes {
    my $self = shift;
    return (keys %{$self->{modes}});
}
#-------------------------------------------------------------------------------
sub execute {
    my $self     = shift;
    my $hashref  = shift;
    my $args     = shift;
    my $postproc = shift;

    my $cmd = $self->make_cmd($hashref, $args, $postproc);

    return system($cmd);
}
#-------------------------------------------------------------------------------
sub make_cmd {
    my $self     = shift;
    my $hashref  = shift;
    my $args     = shift;
    my $postproc = shift;

    my $cmd = $self->{cur_cmdline};
    my $val = 0;

    foreach my $param (@{$self->{cur_param_names}}) {
        if ((! exists $$hashref{$param}) || ($$hashref{$param} eq "")) {
            print("ERROR: Program::make_cmd: $param is not defined\n");
            return undef;
        }
        $val = $$hashref{$param};
        $cmd =~ s/\b$param\b/$val/g;
    }
    $val =  $self->{exe};
    $cmd =~ s/\bexe\b/$val/g;

    $cmd .= " $args"     if (defined $args);
    $cmd .= " $postproc" if (defined $postproc);

    return $cmd;
}
#-------------------------------------------------------------------------------
sub get_getopt_strings {
    my $self    = shift;

    my $hashref = $self->{all_param_to_types_hash};

    my @getopt_strings = ();
    my $tmp = "";

    foreach my $param (keys %$hashref) {
        $tmp = $param;
        if ($$hashref{$param} ne "switch") {
            if ($$hashref{$param} eq "string") {
                $tmp .= '=s';
                push(@getopt_strings, $tmp);
            }
            if ($$hashref{$param} eq "int") {
                $tmp .= '=i';
                push(@getopt_strings, $tmp);
            }
            if ($$hashref{$param} eq "extint") {
                $tmp .= '=o';
                push(@getopt_strings, $tmp);
            }
            if ($$hashref{$param} eq "float") {
                $tmp .= '=f';
                push(@getopt_strings, $tmp);
            }
        } else {
            push(@getopt_strings, $tmp);
        }
    }
    #
    # In addition to the parameters per se, also add in this list options for
    # choosing the mode to use...
    #
    my @mode_list = $self->get_all_modes();
    foreach my $mode (@mode_list) {
        push(@getopt_strings, "$mode=s");
    }

    return @getopt_strings;
}
#-------------------------------------------------------------------------------

1;

__END__


#-------------------------------------------------------------------------------
# Stub documentation for this module.
#-------------------------------------------------------------------------------

=head1 NAME

Tifa::Program - An abstraction of a TIFA command line program

=head1 SYNOPSIS

  use Tifa::Program;
  $program = new Tifa::Program();

=head1 REQUIRE

Perl 5.006002, Carp, Class::Struct and Exporter.

=head1 SUMMARY

The Tifa::Program module is a weird module whose only raison d'etre is to
provide an abstraction of what a launchable and benchmarkable command line
program is. Most, if not all, of TIFA users do not need to known about this
module.

=head1 DESCRIPTION

The Tifa::Program module provide an abstraction of what a launchable and
benchmarkable command line program is. In less cryptic words, its goal is to
provide other TIFA's scripts and modules with a unified way to describe and
interact with the implemented factorization programs. The motivation for such a
module is to be able to add easily new (factorization) programs that can be
launched or benchmarked by TIFA's scripts while containing code changes to a few
places only. As a TIFA user you probably do not need to known anything about
this module unless you plan to use TIFA's benchmarking and plotting framework
for programs other than TIFA's included factorization programs.

=head2 Data members

A Program object is defined by the following attributes:

    $algo  : Name of the algorithm implemented.

    $descr : Short description of the algorithm.

    $exe   : Name of the command line program that is wrapped.

    $help  : Help message for the command line program that is wrapped.

    %modes : An hashtable mapping mode names to 'modes'. A mode
             represents a particular form of invocation of the
             program given by $exe, and is defined by:

                 - A description $descr
                 - A list of relevant parameter names
                 - A command line template

             For example, TIFA's CFRAC program can be called with
             different numbers of parameter according to whether or
             not we should let the program choose optimal parameter
             values. In this case, we could have a mode where all
             parameters have to be specified and another one where
             only a few of them are needed. (See example in the next
             section.)

    $cur_mode : Current mode name.

    @cur_param_names : List of parameter names relevant to the current
                       mode.

    @cur_param_descrs: List of parameter descriptions relevant to the
                       current mode. The descriptions are given in the
                       same parameter order than in @cur_param_names.

    @cur_param_types : List of parameter types relevant to the current
                       mode. The types are given in the same parameter
                       order than in @cur_param_names.

    $cur_cmdline: Command line template relevant to the current mode.

    $default_mode: Default mode name.

    @all_param_names : Array listing all of the parameter names

    %all_param_to_descrs_hash: Hashtable mapping the parameter names
                               to their descriptions;

    %all_param_to_types_hash : Hashtable mapping the parameter names
                               to their types

=head2 Available methods

This module provides the following methods:

    new()
    set_algo($algo)
    get_algo()
    set_descr($descr)
    get_descr()
    set_exe($exe)
    get_exe()
    set_help($help)
    get_help()
    set_default_mode($mode)
    get_default_mode()
    add_mode($key, $name, $descr,
             \@parnames, \@pardescs, \@partypes, $cmdline)
    set_mode($mode)
    get_mode()
    get_all_param_names()
    get_all_param_to_descrs_hash()
    get_all_param_to_types_hash()
    get_param_names()
    get_param_descrs()
    get_param_types()
    get_cmdline()
    get_mode_descr()
    get_all_modes()
    get_getopt_strings()
    execute(\%hash, $arg, $postcmd)
    make_cmd(\%hash, $arg, $postcmd)

=head2 Methods description

    new()
        Basic constructor allocating a Tifa::Program object.

    set_algo($algo)
        Sets the name of the algorithm implemented by the Program object
        to $algo.

    get_algo()
        Gets the name of the algorithm implemented by the Program object.

    set_descr($descr)
        Sets the description of the algorithm implemented to $descr.

    get_descr()
        Gets the description of the algorithm implemented.

    set_exe($exe)
        Sets the name of the command line program to $exe.

    get_exe()
        Gets the name of the command line program.

    set_help($help)
        Sets the help message to $help.

    get_help()
        Gets the help message as a string.

    set_default_mode($mode)
        Sets the name of the default mode of the Program object to $mode.
        Exits if mode $mode does not exist.

    get_default_mode()
        Gets the name of the default mode of the Program object.

    add_mode($key, $name, $descr,
             \@parnames, \@pardescs, \@partypes, $cmdline)
        Adds a mode to the Program object. The mode is identified by:

            - a name $name
            - a description string $descr,
            - a list of relevant parameter names given by a reference to
              an array @parnames
            - a list of relevant parameter descriptions given by a
              reference to an array @pardescs
            - a list of relevant parameter types given by a reference to
              an array @partypes
            - a template of the command line $cmdline.

        The parameter order used in the arrays @parnames, @pardescs,
        and @partypes are assumed to be the same, ie: $parnames[$i],
        $pardescs[$i] and $partypes[$i] refer to the same parameter.

        The type of a parameter is given by a string with one of the
        following value:

            - 'string'  if the parameter is a string
            - 'int'     if the parameter is a positive integer
            - 'extint'  if the parameter is an extended integer
            - 'float'   if the parameter is a float
            - 'switch'  if the parameter is just an option switch

    set_mode($mode)
        Sets the current mode to the one name $mode. Reverts to the
        default mode if no mode named $mode exists.

    get_mode()
        Gets the name of the current mode.

    get_all_param_names()
        Returns an array containing all the potentially relevant
        parameter names, regardless of the current mode.

    get_all_param_to_descrs_hash()
        Returns a hashtable mapping all the parameter names (regardless
        of the current mode), to their descriptions.

    get_all_param_to_types_hash()
        Returns a hashtable mapping all the parameter names (regardless
        of the current mode), to their types.

    get_param_names()
        Returns an array containing all the relevant parameter names
        in the current mode.

    get_param_descrs()
        Returns an array containing all the relevant parameter
        descriptions in the current mode.

    get_param_types()
        Returns an array containing all the relevant parameter types
        in the current mode.

    get_cmdline()
        Returns the command line template in the current mode.

    get_mode_descr($mode)
        Gets the description of the mode named $mode.

    get_all_modes()
        Returns an array containing all the mode names.

    get_getopt_strings()
        Returns an array of strings that can be passed as parameters to
        the GetOptions function of the Getopt::Long module. These strings
        are of the form <param_name>(=(s|i|o|f))?. For exemple: "help",
        "infile=s", etc. All of the parameters are included here,
        regardless of the current program's mode. Also options for
        selecting the program's mode are included. These are of the form
        <mode_name>=s.

    execute(\%param_values, $args, $postcmd)
        Executes the command line program $self->{exe} in the current
        mode with the parameter values being given by the hashtable
        %param_values (whose reference is passed as a parameter) for the
        arguments given in $args. Note that the %param_values hashtable
        should map the parameter names to their actual values.
        Optionnally, a third argument $postcmd can be passed. It will be
        appended to the command string and can be used to implement
        pipes, i/o redirections, etc.

        Returns 'undef' if some parameters were missing in the hashtable,
        in which case the program is not executed.

        Returns the return-value of the call to the underlying Perl
        'system' function otherwise.

    make_cmd(\%param_values, $args, $postcmd)
        Identical to the execute function except that the command
        obtained is not executed, but is returned as a string.

=head1 EXAMPLE

This following code snippet gives an example of how the Tifa::Program module
can be used.

     #
     # Let's wrap the CFRAC program: cfrac_factors
     #
 $cfrac_program = new Tifa::Program();

 $program->set_algo("cfrac");
 $program->set_descr("This is the Continued FRACtion algorithm");
 $program->set_exe("cfrac_factors");

 @params = (
     "exe",
     "nprimes_in_factor_base",
     "nprimes_tdiv_smooth_nb",
     "nrelations",
     "lsr_method",
     "use_large_primes",
     "nprimes_tdiv",
 );
     #
     # The command line template can be used down the road to produce
     # the exact command to launch the factoring program
     #
 $cmdline =  "exe nprimes_in_factor_base "
            ."nprimes_tdiv_smooth_nb nrelations lsr_method "
            ."use_large_primes nprimes_tdiv ";
     #
     # Default mode: We should specify all of the parameter values on
     # the command line
     #
 $program->add_mode("default",
                    "Default mode: specify all parameter values",
                    @params, $cmdline);

 $program->set_default_mode("default");

 @params  = ("exe", "nprimes_tdiv");
 $cmdline =  "exe nprimes_tdiv";
     #
     # Best mode: Let CFRAC use the precomputed optimal parameter values
     # depending on the size of the number to factor
     #
 $program->add_mode("best", "Best mode: let CFRAC choose optimal values",
                    @params, $cmdline);

 $program->set_mode("best");
     #
     # For example, all of the factoring programs are wrapped in such
     # a way in the Tifa::ProgramRepository module, so that the
     # other TIFA scripts and/or modules can use them without having to
     # know any specifics about the programs. For example, this makes
     # possible to reuse the benchmark framework with an user-defined
     # Tifa::Program without delving into the scripts' internals.
     #

=head1 EXPORT

No functions are exported from this package by default.

=head1 SEE ALSO

The Tifa::ProgramRepository module.

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
