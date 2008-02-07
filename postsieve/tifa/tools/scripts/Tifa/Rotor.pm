#
# Copyright (C) 2006, 2007, 2008 INRIA (French National Institute for Research 
# in Computer Science and Control)
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

package Tifa::Rotor;

#-------------------------------------------------------------------------------
#                A counting rotor for looping over combination.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : Rotor.pm
# Author        : Jerome Milan
# Created on    : Fri Feb  2 2007
# Last modified : Fri Feb  2 2007
#
# Version : 0.1.0
# License : GNU Lesser General Public License (LGPL) v2.1 or later
#           Copyright (C) 2006, 2007, 2008 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# An odometer-like counter that can be used to generate all possible
# combinations from sets of values. In TIFA, it is used in the benchmark.pl and
# plotmaker.pl scripts to loop through all of the possible combinations of the
# parameter values to either 1) generate an exhaustive benchmark for a wide
# set of inputs or 2) generate plots for each of these combinations.
#-------------------------------------------------------------------------------

use 5.006002;
use strict;
use warnings;
use Carp;

require Exporter;

our @ISA       = qw(Exporter);
our @EXPORT_OK = qw();
our $VERSION   = '0.1.0';

#-------------------------------------------------------------------------------
sub new {
    my $type     = shift(@_);
    my @wnstates = @_;
    my $self = {};

    $self->{nwheels}       = scalar(@wnstates);
    $self->{wheel_nstates} = ();
    $self->{state}         = ();

    bless $self, $type;

    push(@{$self->{wheel_nstates}}, @wnstates);

    $self->reset();

    return $self;
}
#-------------------------------------------------------------------------------
sub reset {
    my $self = shift;

    foreach my $istate (0 .. ($self->{nwheels} - 1)) {
        $self->{state}[$istate] = 0;
    }
}
#-------------------------------------------------------------------------------
sub get_state {
    my $self = shift;

    return @{$self->{state}};
}
#-------------------------------------------------------------------------------
sub set_state {
    my $self  = shift;
    my @new_state = @_;

    if (scalar(@new_state) != scalar(@{$self->{state}})) {
        croak("Rotor::set_state: incompatible input state length");
    }
    foreach my $iw (0 .. ($self->{nwheels} - 1)) {
        if ($new_state[$iw] < 0) {
            croak("Rotor::set_state: state values must be positive");
        }
        if ($new_state[$iw] >= $self->{wheel_nstates}[$iw]) {
            croak("Rotor::set_state: incompatible input state value");
        }
        $self->{state}[$iw] = $new_state[$iw];
    }
}
#-------------------------------------------------------------------------------
sub increment {
    my $self = shift;

    $self->__increment_wheel__($self->{nwheels} - 1);

    return;
}
#-------------------------------------------------------------------------------
sub __increment_wheel__ {
    my $self  = shift;
    my $wheel = shift;
    my $next  = 0;

    $self->{state}[$wheel] = $self->{state}[$wheel] + 1;

    if ($self->{state}[$wheel] == $self->{wheel_nstates}[$wheel]) {

        if ($wheel == 0) {
            $next = $self->{nwheels} - 1;
        } else {
            $next = $wheel - 1;
        }
        $self->__increment_wheel__($next);
        $self->{state}[$wheel] = 0;
    }
    return;
}
#-------------------------------------------------------------------------------
sub is_zero {
    my $self  = shift;

    foreach my $val (@{$self->{state}}) {
        if ($val != 0) {
            return 0;
        }
    }
    return 1;
}
#-------------------------------------------------------------------------------

1;

__END__


#-------------------------------------------------------------------------------
# Stub documentation for this module.
#-------------------------------------------------------------------------------

=head1 NAME

Tifa::Rotor - A basic odometer-like rotor/counter

=head1 SYNOPSIS

  use Tifa::Rotor;
  $rotor = new Tifa::Rotor(@array);

=head1 REQUIRE

Perl 5.006002, Carp, and Exporter.

=head1 SUMMARY

The Tifa::Rotor module is yet another weird TIFA module implementing an
odometer-like counter that can be used to generate all possible n-tuples from
n sets of values.

=head1 DESCRIPTION

The Tifa::Rotor module implements an odometer-like counter that can be used to
generate, step by step, all possible n-tuples (a, b, ..., n) from n sets: A={0,
1, ..., a_max}, B={0, 1, ..., b_max}, ..., N={0, 1, ..., n_max}.

In TIFA, it is used in the benchmark.pl and plotmaker.pl scripts to loop through
all of the possible combinations of the parameter values to either: 1) generate
an exhaustive benchmark for a wide set of inputs or 2) generate plots for each
of these combinations.

=head2 Available methods

This module provides the following methods:

    new(@array)
    reset()
    get_state()
    set_state(@array)
    increment()
    is_zero()

=head2 Methods description

    new(@array)
        Basic constructor allocating a Tifa::Rotor object from an array
        of integer @array. The length of @array gives the number of
        "wheels" of the rotor, or in other words, the size of the
        n-tuples generated. The values of @array sets the upper limit
        on each value of a given n-tuple. For example, if
        @array = (5, 2, 7) then the 3-tuples generated will be of
        the form:
            (a, b, c) with 0 <= a < 5; 0 <= b < 2 and 0 <= c < 7.

    reset()
        Resets the state of the rotor to the (0, ..., 0) n-tuple.

    get_state()
        Returns the current state of the rotor as an array of integer.

    set_state(@array)
        Sets the current state of the rotor to @array. @array must have
        the required length and its values must be consistent with the
        ranges provided to the constructor or an error will be raised.

    increment()
        Increments the rotor, changing its state and thus providing a
        new n-tuple.

    is_zero()
        Returns 1 if the rotor's current state is (0, ..., 0).
        Returns 0 otherwise.

=head1 EXAMPLE

This following code snippet gives an example of how (and why) the Tifa::Rotor
module can be used... No, it is not _that_ useless!

 #
 # In this example, we'd like to perform some operations on every
 # possible combination of the velocity, angle and height values.
 #
 %values = (
     "velocity"=>[1, 3, 5],    # velocity can take one of these 3 values
     "angle"   =>[20, 45, 70], # angle can take one of these 3 values
     "height"  =>[10, 20]      # height can take one of these 2 values
 );
 @params = keys %values;

 @sizes = ();
 foreach $p (@params) {
     push(@sizes, scalar @{$values{$p}} );
 }

 $rotor = new Tifa::Rotor(@sizes);

 do {
     @ntuple = $rotor->get_state(); # initial state is (0, 0, 0)
     %combo  = ();

     foreach $i (0 .. scalar(@ntuple)-1) {
         #
         # Fetch the parameter values using the indices given by the
         # rotor's state.
         #
         $parname = $params[$i];
         $combo{$parname} = $values{$parname}[$ntuple[$i]];
     }
     #
     # Now %combo is an hash holding a possible combination of
     # the parameter values. It could be represented as:
     #
     # %combo = ("velocity" => <X>, "angle" => <Y>, "height" => <Z>);
     #
     # Do something with it before stepping to the next combination...
     #
     $rotor->increment(); # step to the next rotor state
 } while (! $rotor->is_zero);

In this particularly simple example, there is no real use for the Tifa::Rotor
module as three "for" loops would do the job more straitforwardly. However in
more general cases the %values could be extracted from a configuration file,
passed as a parameter, etc., so using a Rotor object makes for a much more
reuseable code.

=head1 EXPORT

No functions are exported from this package by default.

=head1 SEE ALSO

=head1 AUTHOR

Jerome Milan, E<lt>milanj@lix.polytechnique.frE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2006, 2007, 2008 INRIA (French National Institute for Research
in Computer Science and Control)

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
