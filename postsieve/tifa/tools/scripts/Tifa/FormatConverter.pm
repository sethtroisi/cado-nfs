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

package Tifa::FormatConverter;

#-------------------------------------------------------------------------------
#                A wrapper for format conversion utilities.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : FormatConverter.pm
# Author        : Jerome Milan
# Created on    : circa late July / early August 2006
# Last modified : circa August 2006
#
# Version : 0.1
# License : GNU Lesser General Public License (LGPL) v2.1 or later
#           Copyright (C) 2006, 2007, 2008 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# The Tifa::FormatConverter is a mere wrapper for several conversion
# utilities popular in the unix/unix-like world, namely convert (from
# the very powerful ImageMagick suite), ps2pdf, ps2eps, pdf2ps, epstopdf,
# dvips and dvipdf.
#-------------------------------------------------------------------------------

use 5.006002;
use strict;
use warnings;
use Carp;

require Exporter;

our @ISA        = qw(Exporter);
our @EXPORT_OK  = qw();
our $VERSION    = '0.01';

#-------------------------------------------------------------------------------
# Class variables
#-------------------------------------------------------------------------------

#
# Possible convertions for the convert utility from the ImageMagick suite
#
my %in_out_convert = (
    "pdf"  => ["png", "jpg", "jpeg","gif",  "tiff"],
    "eps"  => ["png", "jpg", "jpeg","gif",  "tiff"],
    "ps"   => ["png", "jpg", "jpeg","gif",  "tiff"],
    "png"  => ["pdf", "eps", "ps",  "jpeg", "gif",  "tiff"],
    "jpeg" => ["pdf", "eps", "ps",  "png",  "gif",  "tiff"],
    "gif"  => ["pdf", "eps", "ps",  "jpg",  "jpeg", "png", "tiff"],
    "tiff" => ["pdf", "eps", "ps",  "jpg",  "jpeg", "gif", "png"]
);
#
# Possible convertions for the epstopdf utility
#
my %in_out_epstopdf = (
    "eps"  => ["pdf"]
);
#
# Possible convertions for the ps2pdf utility
#
my %in_out_ps2pdf = (
    "ps"   => ["pdf"]
);
#
# Possible convertions for the ps2eps utility
#
my %in_out_ps2eps = (
    "ps"   => ["eps"]
);
#
# Possible convertions for the pdf2ps utility
#
my %in_out_pdf2ps = (
    "pdf"   => ["ps"]
);
#
# Possible convertions for the dvips utility
#
my %in_out_dvips = (
    "dvi"   => ["ps"]
);
#
# Possible convertions for the dvipdf utility
#
my %in_out_dvipdf = (
    "dvi"   => ["pdf"]
);
#
# Possible convertions for all programs
#
my %program_in_out = (
    "convert"  => \%in_out_convert,
    "epstopdf" => \%in_out_epstopdf,
    "ps2pdf"   => \%in_out_ps2pdf,
    "ps2eps"   => \%in_out_ps2eps,
    "pdf2ps"   => \%in_out_pdf2ps,
    "dvips"    => \%in_out_dvips,
    "dvipdf"   => \%in_out_dvipdf
);
#
# Place holders in command usage
#
my $inholder  = '<infile>';
my $outholder = '<outfile>';
#
# Possible convertions for all programs
#
my %program_usage = (
    "convert"  => "convert $inholder $outholder",
    "epstopdf" => "epstopdf $inholder --outfile=$outholder",
    "ps2pdf"   => "ps2pdf $inholder $outholder",
    "ps2eps"   => "ps2eps $inholder",
    "pdf2ps"   => "pdf2ps $inholder $outholder",
    "dvips"    => "dvips $inholder -o $outholder",
    "dvipdf"   => "dvips $inholder -o $outholder"
);
#
# Name of convertion programs
#
my @programs = (keys %program_in_out);

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
sub new {
    my $type = shift(@_);
    my $self = {};

    bless $self, $type;

    return $self;
}
#-------------------------------------------------------------------------------
sub convert {
    my $self = shift(@_);
    my $in   = shift(@_);
    my $out  = shift(@_);

    my $error = "Tifa::FormatConverter:convert:";

    (my $in_base, my $in_ext)   = ($in  =~ /(.*)\.([^\.]*)$/);
    (my $out_base, my $out_ext) = ($out =~ /(.*)\.([^\.]*)$/);

    my $program = $self->get_program($in_ext, $out_ext);

    if (!$program) {
        croak("$error Unsupported convertion\n");
    }
    #
    # Some programs (i.e. ps2eps) do not take the name of the output file
    # as a parameter (a real shame! :-) but infer it. So we have to make an
    # exception for them and explicitely rename the generated file to the
    # desired output name. Since I don't want to mantain a list of
    # exceptions, I just treat all programs on an equal footing.
    # Lame, I know...
    #
    my $tmp_out = "$in_base.$out_ext";
    my $command = $program_usage{$program};
    $command =~ s/$inholder/$in/;
    $command =~ s/$outholder/$tmp_out/;

    system("$command");

    if ($tmp_out ne $out) {
        rename($tmp_out, $out) or croak("$error Cannot rename file $tmp_out\n");
    }

}
#-------------------------------------------------------------------------------
sub get_program {
    #
    # Deduce the good program name from the extensions of the input and
    # output file.
    #
    my $self    = shift(@_);
    my $in_ext  = shift(@_);
    my $out_ext = shift(@_);

    foreach my $progname (keys %program_in_out) {
        my $hashref = $program_in_out{$progname};
        foreach my $in_format (keys %$hashref) {
            next if ($in_format ne $in_ext);
            my $arrayref = ${$hashref}{$in_format};
            foreach my $out_format (@$arrayref) {
                if ($out_format eq $out_ext) {
                    return $progname;
                }
            }
        }
    }
    return "";
}
#-------------------------------------------------------------------------------


1;

__END__

#-------------------------------------------------------------------------------
# Stub documentation for this module.
#-------------------------------------------------------------------------------

=head1 NAME

Tifa::FormatConverter - A wrapper for format conversion utilities

=head1 SYNOPSIS

 use Tifa::FormatConverter;

 my $converter = new Tifa::FormatConverter();
 $converter->convert($from, $to);

=head1 REQUIRE

Perl 5.006002, Carp, Exporter, Tifa::DataDescriptor and the following programs:
convert (from the ImageMagick suite), ps2pdf, ps2eps, pdf2ps, epstopdf, dvips
and dvipdf.

=head1 DESCRIPTION

The Tifa::FormatConverter is a mere wrapper for several conversion utilities
popular in the unix/unix-like world, namely convert (from the very powerful
ImageMagick suite), ps2pdf, ps2eps, pdf2ps, epstopdf, dvips and dvipdf.

File formats are directly inferred from the extension of the filenames. The
following extensions/conversions are supported:

 - Pixmap to pixmap/(pseudo-)vectorial (uses convert):

     "png"  to "pdf", "eps", "ps", "jpg", "jpeg", "gif" and "tiff"

     "gif"  to "pdf", "eps", "ps", "jpg","jpeg", "png" and "tiff"

     "tiff" to "pdf", "eps", "ps", "jpg","jpeg", "gif" and "png"

     "jpeg" to "pdf", "eps", "ps", "png", "gif" and "tiff"

     "jpg"  to "pdf", "eps", "ps", "png", "gif" and "tiff"

 - Vectorial to pixmap (uses convert):

     "ps"  to "png", "jpg", "jpeg", "gif" and "tiff"

     "eps" to "png", "jpg", "jpeg", "gif" and "tiff"

     "pdf" to "png", "jpg", "jpeg", "gif" and "tiff"

 - Vectorial to vectorial (uses either ps2eps, ps2pdf, epstopdf,
   pdf2ps, dvips or dvipdf):

     "ps"  to "eps" and "pdf"

     "eps" to "pdf"

     "pdf" to "ps"

     "dvi" to "ps" and "pdf"

The Tifa::FormatConverter module does not check whether or not these programs
are installed on your system. Consequently, one can still use this module if
some of the programs are not available, provided that the desired format
conversion does not rely on them.

=head2 Available methods

This module provides the following methods:

    new()
    convert($from, $to)

=head2 Methods description

    new()
        Basic constructor allocating a Tifa::FormatConverter object.

    convert($from, $to)
        Converts the file $from to a file $to, provided that the format
        conversion is supported. Formats are inferred from the filenames'
        extensions.

=head1 EXAMPLE

Using the Tifa::FormatConverter is completely straightforward. For example:

    my $converter = new Tifa::FormatConverter();

    my $eps_file   = "file.eps";
    my $pdf_file   = "file.pdf";
    my $png_file   = "file.png";
    my $jpg_file   = "file.jpg";

    $converter->convert($eps_file, $pdf_file);
    $converter->convert($pdf_file, $png_file);
    $converter->convert($png_file, $jpg_file);

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

