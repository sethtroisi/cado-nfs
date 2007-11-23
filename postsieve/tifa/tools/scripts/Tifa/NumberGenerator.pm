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

package Tifa::NumberGenerator;

#-------------------------------------------------------------------------------
#                  Generate prime or composite integers.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : NumberGenerator.pm
# Author        : Jerome Milan
# Created on    : circa late July / early August 2006
# Last modified : circa August 2006
#
# Version : 0.1
# Licence : GNU Lesser General Public License (LGPL) v2.1 or later
#           Copyright (C) 2006, 2007 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# This Perl module is part of the TIFA (Tools for Integer FActorization)
# library.  It uses the GNU Multi-Precision Perl modules GMP::Mpz and
# GMP::Rand to generate (possibly random) prime or composite non-square
# numbers.
#-------------------------------------------------------------------------------

use 5.006002;
use strict;
use warnings;
use Carp;
#
# Use the GMP::Mpz and GMP::Rand modules. It these modules are not available,
# die (well, croak actually) with an error message...
#
eval "use GMP::Mpz qw(:all); 1"
    or croak "ERROR: Module GMP::Mpz not found!\n";
eval "use GMP::Rand qw(:all); 1"
    or croak "ERROR: Module GMP::Rand not found!\n";

require Exporter;

our @ISA        = qw(Exporter);
our @EXPORT_OK  = qw();
our $VERSION    = '0.01';

#-------------------------------------------------------------------------------
sub new {
    my $type = shift(@_);
    my $self = {};

    $self->{rstate} = undef;

    bless $self, $type;
    return $self;
}
#-------------------------------------------------------------------------------
sub use_random {
    my $self    = shift(@_);
    $self->{rstate} = randstate();
    $self->{rstate}->seed(time());
}
#-------------------------------------------------------------------------------
sub dont_use_random {
    my $self = shift(@_);
    $self->{rstate} = undef;
}
#-------------------------------------------------------------------------------
sub generate_prime {
    #
    # Generate a prime of $size bits. The prime is either choosen randomly,
    # or in a deterministic fashion if $self->{rstate} is undefined. In the
    # latter case the prime returned is the smallest prime greater that
    # 2^($size-1).
    #
    my $self   = shift(@_);
    my $size   = shift(@_);
    my $rstate = $self->{rstate};

    my $n = mpz(2)**($size-1);

    if (defined $rstate) {
        $n += mpz_urandomm($rstate, mpz(2)**($size-1));
    }
    my $prime = nextprime($n);

    while (sizeinbase($prime, 2) != $size) {
        $prime = $self->generate_prime($size);
    }
    return $prime;
}
#-------------------------------------------------------------------------------
sub generate_composite {
    #
    # Generate a composite integer of $size bits from the multiplication of
    # $nb_factors primes. The primes are either choosen randomly, or in a
    # deterministic fashion if $self->{rstate} is undefined.
    #
    my $self       = shift(@_);
    my $size       = shift(@_);
    my $nb_factors = shift(@_);
    my $rstate     = $self->{rstate};

    if ($nb_factors == 1) {
        return $self->generate_prime($size);
    }

    #
    # With the way the prime numbers are generated here, we may "lose"
    # a bit for some multiplication steps. I.e. if size(n1) = x bits
    # and size(n2) = y bits, it is possible to have size(n1 * n2) = x + y - 1
    # bits... This could potentially lead to a larger last factor (to
    # compensate this loss and obtain the required composite's size) but we'd
    # like to keep each factors' size as similar as possible. The next trick
    # (adding a bit to the composite's size for each multiplication needed),
    # although not rigorous, helps to avoid too much discrepancies between
    # the factors' size...
    #
    my $length = $size + $nb_factors - 1;
    my $faclength = int($length/$nb_factors);
    my $lastfaclength = $length - (($nb_factors - 1) * $faclength);
    my $result;

    if ($rstate) {
        #
        # The primes are choosen at random, provided that they have the
        # required bit lengths.
        #
        my $composite = mpz(1);
        for (1..($nb_factors - 1)) {
            $composite *= $self->generate_prime($faclength);
        }
        $lastfaclength = $size - sizeinbase($composite, 2);
        $result = $composite * $self->generate_prime($lastfaclength);
        #
        # Make sure that the composite has the required size and is not
        # a square.
        #
        while ((sizeinbase($result, 2) != $size) || perfect_square_p($result)) {
            $result = $self->generate_composite($size, $nb_factors);
        }
    } else {
        #
        # Here, the prime selection is not random: the next few primes
        # following a power of 2 are used in a deterministic fashion. This
        # facilitates the comparison of benchmark results for different
        # parameter values.
        #
        my $composite = mpz(1);
        my $factor    = $self->generate_prime($faclength);
        for (1..($nb_factors - 1)) {
            $composite *= $factor;
            $factor = nextprime($factor);
        }
        if ($lastfaclength == $faclength) {
            #
            # This check is necessary to make sure we don't include in the
            # product a prime number already used. This could lead to the
            # composite being a square and cause problems with some algorithms
            # such as CFRAC. By making the composite a product of _different_
            # primes, we are sure that we will never run into this problem.
            #
            $result = $composite * $factor;
        } else {
            $factor = $self->generate_prime($lastfaclength);
            #
            # We are basically sure that $factor was not used in the product,
            # unless the number of factors is of the order of the number of
            # primes between 2^($faclength-1) and 2^($lastfaclength-1), a
            # situation which, for all intends and purposes, will never happen.
            #
            $result = $composite * $factor;
        }
    }
    return $result;
}
#-------------------------------------------------------------------------------

1;
__END__

#-------------------------------------------------------------------------------
# Stub documentation for the module.
#-------------------------------------------------------------------------------

=head1 NAME

Tifa::NumberGenerator - Generate prime or composite integers.

=head1 SYNOPSIS

  use Tifa::NumberGenerator;
  use GMP::Mpz;
  $generator = new Tifa::NumberGenerator();

=head1 REQUIRE

Perl 5.006002, Carp, Exporter, GMP::Mpz and GMP::Rand.

=head1 DESCRIPTION

This Perl module is part of the TIFA (Tools for Integer FActorization) library.
It uses the GNU Multi-Precision Perl modules GMP::Mpz and  GMP::Rand to generate
(possibly random) prime or composite non-square numbers.

=head2 Available methods

This module provides the following methods:

    Tifa::NumberGenerator()
    use_random()
    dont_use_random()
    generate_prime($size)
    generate_composite($size, $nb_factors)

=head2 Methods description

    Tifa::NumberGenerator()
        Basic constructor allocating a Tifa::NumberGenerator object.

    use_random()
        Invoke this method to generate numbers randomly.

    dont_use_random()
        Invoke this method to generate numbers deterministically.

    generate_prime($size)
        Returns a prime of $size bits as a GMP::Mpz object.

    generate_composite($size, $nb_factors)
        Returns a non-square composite integer of $size bits (obtained
        by multiplying $nb_factors prime integers of roughly equal sizes)
        as a GMP::Mpz object.

=head2 A note on deterministic generation

(Pseudo-) random generation can be toggled on and off by using the use_random()
and dont_use_random() method.

When a prime is generated deterministically (i.e. when not using a random
generator), the smallest prime of the required bit length is returned.

When a composite integer is generated deterministically, the smallest prime
numbers of the required bit length are multiplied together to yield the
final composite. In this case, each prime number is used only once in the
multiplication.

=head2 A note on primality

The Tifa::NumberGenerator() uses the GMP:Mpz module to check for the primality
of a given number. These tests are actually Miller-Rabin tests of compositeness
rather than primality tests strictly speaking. However, "prime" numbers
generated by the Tifa::NumberGenerator() have indeed a very good chance of
being truly prime, the probability of returning a composite being extremely
small.

=head1 EXAMPLE

This following code snippet gives an example of how to use the
Tifa::NumberGenerator() module.

    my $generator = new Tifa::NumberGenerator();
    #
    # Generate numbers randomly with $generator->use_random() or in
    # a deterministic fashion with $generator->dont_use_random().
    #
    $generator->use_random();
    #
    # Generate a prime of 120 bits.
    #
    my $prime = $generator->generate_prime(120);
    #
    # Generate a composite of 120 bits from the multiplication
    # of 4 prime numbers of (roughly) equal lengths.
    #
    my $composite = $generator->generate_composite(120, 4);

=head1 EXPORT

No functions are exported from this package by default.

=head1 SEE ALSO

The GNU Multi-Precision Perl modules GMP::Mpz and GMP::Rand.

=head1 AUTHOR

Jerome Milan, E<lt>milanj@lix.polytechnique.frE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2006 INRIA (French National Institute for Research in Computer
Science and Control)

This module is part of the TIFA (Tools for Integer FActorization) library.

The licence under which TIFA will be distributed is still to be defined.
In the interim the TIFA library should not be distributed outside of the TANC
group of the French National Institute for Research in Computer Science and
Control.

=cut
