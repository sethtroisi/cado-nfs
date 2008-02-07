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

package Tifa::GnuPlotter;

#-------------------------------------------------------------------------------
#                Generate 2D plots from 3D data with gnuplot.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : GnuPlotter.pm
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
# The Tifa::GnuPlotter module facilitates the generation of 2D plots from
# 3D data by using the gnuplot program as its underlying engine. However,
# it is not intended to be a generic low-level gnuplot wrapper like
# Chart::Graph for example since it is way too limited for that.
#-------------------------------------------------------------------------------

use 5.006002;
use strict;
use warnings;
use Carp;
use File::Basename;

use Tifa::DataDescriptor;
use Tifa::FormatConverter;

require Exporter;

our @ISA        = qw(Exporter);
our @EXPORT_OK  = qw();
our $VERSION    = '0.01';

#-------------------------------------------------------------------------------
# Class variables
#-------------------------------------------------------------------------------

#
# Gnuplot macro extension.
#
our $macro_extension = "gp";
#
# Allowed output file format for the GnuPlotter. There's actually many more
# but the most useful are these ones.
#
our %allowed_output_format = (
    "pdf" => "generate Portable Document Format files",
    "eps" => "generate Encapsulated PostScript files",
    "png" => "generate Portable Network Graphics files"
);
#
# Pivot format for the output file. A temporary file for some format is
# created to bypass some technical restrictions (at the price of portability).
#
# For example, the PDF support via "set terminal pdf" is relatively new and
# not usable out of the box for most gnuplot versions (even the 4.0). So we
# keep the postscript terminal and explicitly convert the generated EPS file
# to the desired PDF format.
#
my %pivot_output_format = (
    "pdf" => "eps",
    "eps" => "",
    "png" => ""
);
#
# Extensions for the allowed file format.
#
my %extensions = (
    "pdf" => "pdf",
    "eps" => "eps",
    "png" => "png"
);
#
# The "set terminal" commands for each allowed output file format.
# Note that the PDF format will need an explicit conversion from EPS.
#
my %set_terminal_cmds = (
    "pdf" => "set terminal postscript color eps",
    "eps" => "set terminal postscript color eps",
    "png" => "set terminal png medium"
);
#
# Number of predefined line styles
#
my $nline_styles = 20;
#
# Gnuplot code predefining some line styles
#
my $style_defs = << 'EOF';
set style line 1 lt 1 lw 1 pt 1 ps 2
set style line 2 lt 2 lw 1 pt 2 ps 2
set style line 3 lt 3 lw 1 pt 3 ps 2
set style line 4 lt 4 lw 1 pt 4 ps 2
set style line 5 lt 5 lw 1 pt 5 ps 2
set style line 6 lt 6 lw 1 pt 6 ps 2
set style line 7 lt 7 lw 1 pt 7 ps 2
set style line 8 lt 8 lw 1 pt 8 ps 2
set style line 9 lt 9 lw 1 pt 9 ps 2
set style line 10 lt 10 lw 1 pt 10 ps 2
set style line 11 lt 11 lw 1 pt 11 ps 2
set style line 12 lt 12 lw 1 pt 12 ps 2
set style line 13 lt 13 lw 1 pt 13 ps 2
set style line 14 lt 14 lw 1 pt 14 ps 2
set style line 15 lt 15 lw 1 pt 15 ps 2
set style line 16 lt 16 lw 1 pt 16 ps 2
set style line 17 lt 17 lw 1 pt 17 ps 2
set style line 18 lt 18 lw 1 pt 18 ps 2
set style line 19 lt 19 lw 1 pt 19 ps 2
set style line 20 lt 20 lw 1 pt 20 ps 2
EOF

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
sub new {
    my $type = shift(@_);
    my $self = {};

    $self->{gnuplot_program}   = 'gnuplot';
    $self->{plot_format}       = 'pdf';
    $self->{pivot_plot_format} = 'eps';

    $self->{data_separator} = ':';

    $self->{xvar_name}   = '';
    $self->{yvar_name}   = '';
    $self->{zvar_name}   = '';

    $self->{xvar_descr}  = '';
    $self->{yvar_descr}  = '';
    $self->{zvar_descr}  = '';

    $self->{xvar_values} = '';
    $self->{yvar_values} = '';
    $self->{zvar_values} = '';

    $self->{xaxis_title} = '';
    $self->{yaxis_title} = '';
    $self->{graph_title} = '';

    $self->{main_label}   = '';
    $self->{second_label} = '';

    $self->{data_filename}   = '';
    $self->{macro_filename}  = '';
    $self->{output_filename} = '';
    $self->{pivot_filename}  = '';

    $self->{plot_style} = 'points';

    bless $self, $type;

    return $self;
}
#-------------------------------------------------------------------------------
sub set_gnuplot_program {
    my $self = shift(@_);

    $self->{gnuplot_program} = shift(@_);
}
#-------------------------------------------------------------------------------
sub set_output_format {
    my $self   = shift(@_);
    my $format = shift(@_);

    if (exists $allowed_output_format{$format}) {
        $self->{plot_format}       = $format;
        $self->{pivot_plot_format} = $pivot_output_format{$format};
    }
}
#-------------------------------------------------------------------------------
sub set_data_filename {
    my $self = shift(@_);
    my $name = shift(@_);

    $name =~ s/\s*$//;
    $name =~ s/^\s*//;

    $self->{data_filename} = $name;
}
#-------------------------------------------------------------------------------
sub set_macro_filename {
    my $self = shift(@_);
    my $name = shift(@_);

    $name =~ s/\s*$//;
    $name =~ s/^\s*//;

    if ($name !~ m/\.$macro_extension$/) {
        $name .= ".$macro_extension";
    }

    $self->{macro_filename} = $name;
}
#-------------------------------------------------------------------------------
sub set_output_filename {
    my $self = shift(@_);
    my $name = shift(@_);

    my $ext  = $extensions{$self->{plot_format}};

    $name =~ s/\s*$//;
    $name =~ s/^\s*//;

    if ($name !~ m/\.$ext$/) {
        $name .= ".$ext";
    }
    $self->{output_filename} = $name;

    (my $base, my $dirname, my $suffix) = fileparse($name);

    if ($self->{pivot_plot_format}) {
        #
        # The final output plot will be obtained from conversion
        # of a temporary file in the pivot format.
        #
        $self->{pivot_filename}  = "$dirname".".tmp_$base.";
        $self->{pivot_filename} .= "$extensions{$self->{pivot_plot_format}}";
    } else {
        $self->{pivot_filename} = $self->{output_filename};
    }
}
#-------------------------------------------------------------------------------
sub set_xvar_description {
    my $self        = shift(@_);
    my $name        = shift(@_);
    my $description = shift(@_);

    $self->{xvar_name}  = $name;
    $self->{xvar_descr} = $description;
}
#-------------------------------------------------------------------------------
sub set_yvar_description {
    my $self        = shift(@_);
    my $name        = shift(@_);
    my $description = shift(@_);

    $self->{yvar_name}  = $name;
    $self->{yvar_descr} = $description;
}
#-------------------------------------------------------------------------------
sub set_zvar_description {
    my $self        = shift(@_);
    my $name        = shift(@_);
    my $description = shift(@_);

    $self->{zvar_name}  = $name;
    $self->{zvar_descr} = $description;
}
#-------------------------------------------------------------------------------
sub set_xvar_values {
    my $self      = shift(@_);
    my $arrayref  = shift(@_);

    $self->{xvar_values} = $arrayref;
}
#-------------------------------------------------------------------------------
sub set_yvar_values {
    my $self      = shift(@_);
    my $arrayref  = shift(@_);

    $self->{yvar_values} = $arrayref;
}
#-------------------------------------------------------------------------------
sub set_zvar_values {
    my $self      = shift(@_);
    my $arrayref  = shift(@_);

    $self->{zvar_values} = $arrayref;
}
#-------------------------------------------------------------------------------
sub set_xaxis_title {
    my $self  = shift(@_);
    my $title = shift(@_);

    $title =~ s/^(\s)*//;
    $title =~ s/(\s)*$//;

    $self->{xaxis_title} = $title;
}
#-------------------------------------------------------------------------------
sub set_yaxis_title {
    my $self  = shift(@_);
    my $title = shift(@_);

    $title =~ s/^(\s)*//;
    $title =~ s/(\s)*$//;

    $self->{yaxis_title} = $title;
}
#-------------------------------------------------------------------------------
sub set_graph_title {
    my $self  = shift(@_);
    my $title = shift(@_);

    $title =~ s/^(\s)*//;
    $title =~ s/(\s)*$//;

    $self->{graph_title} = $title;
}
#-------------------------------------------------------------------------------
sub set_main_label {
    my $self  = shift(@_);
    my $label = shift(@_);

    $label =~ s/^(\s)*//;
    $label =~ s/(\s)*$//;

    $self->{main_label} = ($label);
}
#-------------------------------------------------------------------------------
sub set_second_label {
    my $self  = shift(@_);
    my $label = shift(@_);

    $label =~ s/^(\s)*//;
    $label =~ s/(\s)*$//;

    $self->{second_label} = $label;
}
#-------------------------------------------------------------------------------
sub write_data_file {
    my $self      = shift(@_);

    my $error     = "ERROR: Tifa::GnuPLotter:write_data_file:";
    my $localtime  = localtime();

    open(OUT, ">$self->{data_filename}")
        or croak("$error Cannot open $self->{data_filename}.\n");
    print OUT << "EOF";
#------------------------------------------------------------------------------
# File   : $self->{data_filename}
# Author : Automatically generated by Tifa::GnuPlotter
# Date   : $localtime
#------------------------------------------------------------------------------
EOF
    my $handle = *OUT;
    my $descr  = new Tifa::DataDescriptor();
    $descr->set_comment_style("Perl");
    $descr->set_field_separator($self->{data_separator});
    $descr->add_field_description($self->{xvar_name}, $self->{xvar_descr});
    $descr->add_field_description($self->{yvar_name}, $self->{yvar_descr});
    $descr->add_field_description($self->{zvar_name}, $self->{zvar_descr});
    $descr->write_descriptions($handle);

    my $main_label = $self->{main_label};
    $main_label =~ s/(\s)*\n(\s)*/\n# /g;
    $main_label =~ s/(\s)*\\n(\s)*/\n# /g;
    $main_label =~ s/^(\s)*/# /;

    my $second_label = $self->{second_label};
    $second_label =~ s/(\s)*\n(\s)*/\n# /g;
    $second_label =~ s/(\s)*\\n(\s)*/\n# /g;
    $second_label =~ s/^(\s)*/# /;

    print OUT << "EOF";
#------------------------------------------------------------------------------
$main_label
#------------------------------------------------------------------------------
$second_label
#------------------------------------------------------------------------------

EOF

    my $xdataref = \@{$self->{xvar_values}};
    my $ydataref = \@{$self->{yvar_values}};
    my $zdataref = \@{$self->{zvar_values}};

    for (my $i = 0; $i <= $#{@$xdataref}; $i++) {
        print OUT "${$xdataref}[$i]$self->{data_separator}";
        print OUT "${$ydataref}[$i]$self->{data_separator}";
        print OUT "${$zdataref}[$i]\n"
    }
    close(OUT);
}
#-------------------------------------------------------------------------------
sub generate_macro {
    my $self       = shift(@_);

    my $error     = "ERROR: Tifa::GnuPLotter:generate_macro:";
    my $localtime  = localtime();

    my $set_term_cmd = "set terminal postscript color eps";

    open(MAC, ">$self->{macro_filename}")
        or croak("$error Cannot create macro $self->{macro_filename}!\n");

	print MAC << "EOF";
#------------------------------------------------------------------------------
# File   : $self->{macro_filename}
# Author : Automatically generated by Tifa::GnuPlotter
# Date   : $localtime
#------------------------------------------------------------------------------
EOF

	print MAC << "EOF";
$set_terminal_cmds{$self->{plot_format}}
set output \"$self->{pivot_filename}\"
set datafile separator \"$self->{data_separator}\"

set autoscale
set multiplot
set grid
set key left Left reverse box

$style_defs
set tmargin 5

set label 1 "$self->{graph_title}" \\
    center \\
    at \\
    screen 0.5, \\
    screen 0.97
show label 1

set arrow nohead \\
lw 1 \\
lt -1 \\
from screen 0.0, screen 0.94 \\
to   screen 1.0, screen 0.94

set xlabel \"$self->{xaxis_title}\"
set ylabel \"$self->{yaxis_title}\"

EOF
    if ($self->{main_label}) {
        print MAC << "EOF";
set label 2 \"$self->{main_label}\" \\
    font "fixed" \\
    center \\
    at \\
    screen 0.5, \\
    screen 0.915
show label 2

EOF
    }
    if ($self->{second_label}) {
        print MAC << "EOF";
set label 3 \"$self->{second_label}\" \\
    font "fixed" \\
    left \\
    at \\
    screen 0.11, \\
    screen 0.6
show label 3

EOF
    }
    my %zhash = ();
    for (@{$self->{zvar_values}}) { $zhash{$_} = 1 }
    my @zvals  = sort {$a <=> $b} (keys %zhash);
    my $nzvals = scalar(@zvals);
    my $linestyle = 0;

    print MAC << "EOF";
plot \\
EOF
    foreach my $zvalue (0..($nzvals-2)) {
        $linestyle = 1 + ($zvalue % ($nline_styles - 1));
        print MAC << "EOF";
     \"$self->{data_filename}\" \\
     using 1:((\$3==$zvals[$zvalue]) ? \$2 : 1/0) \\
     ls $linestyle \\
     title \" $self->{zvar_name} = $zvals[$zvalue]\" \\
     with $self->{plot_style} , \\
\\
EOF
    }
    $linestyle = 1 + (($nzvals - 1) % ($nline_styles - 1));
    print MAC << "EOF";
     \"$self->{data_filename}\" \\
     using 1:((\$3==$zvals[$nzvals-1]) ? \$2 : 1/0) \\
     ls $linestyle \\
     title \" $self->{zvar_name} = $zvals[$nzvals-1]\" \\
     with $self->{plot_style}

EOF

    close(MAC);

}
#-------------------------------------------------------------------------------
sub execute_macro {
    my $self     = shift(@_);

    my $command  = "$self->{gnuplot_program} $self->{macro_filename}";
    $command    .= " >/dev/null";

    system("$command");

    if ($self->{pivot_plot_format}) {
        #
        # Convert the generated plot file to the desired output format and
        # remove the temporary file in pivot format.
        #
        my $converter = new Tifa::FormatConverter();
        $converter->convert($self->{pivot_filename}, $self->{output_filename});
        unlink($self->{pivot_filename});
    }
}
#-------------------------------------------------------------------------------

1;

__END__

#-------------------------------------------------------------------------------
# Stub documentation for this module.
#-------------------------------------------------------------------------------

=head1 NAME

Tifa::GnuPlotter - Generate 2D plots from 3D data with gnuplot

=head1 SYNOPSIS

  use Tifa::GnuPlotter;
  $reader = new Tifa::GnuPlotter();

=head1 REQUIRE

Perl 5.006002, Carp, Exporter, Tifa::DataDescriptor, Tifa::FormatConverter
and the gnuplot program.

=head1 SUMMARY

The Tifa::GnuPlotter module facilitates the generation of 2D plots from 3D
data by using the gnuplot program as its underlying engine. However,
it is not intended to be a generic low-level gnuplot wrapper like
Chart::Graph for example since it is way too limited for that.

=head1 DESCRIPTION

The Tifa::GnuPlotter module was written during the benchmarking phase of the
factoring tools of the TIFA library to facilitate the automatization of
timing plots. Although it was only ment to fulfills that very special
requirement, it can be helpful in other, similar data analysis problems.

Tifa::GnuPlotter is specifically designed to generate 2D plots from
3-dimensional data where the extra third dimension is supposed to be a
varying parameter: the output will be several graphs (actually one for
each value of the third parameter) drawn together on one plot.

Generated plots can be in Portable Document Format (PDF), Encapsulated
PostScript format (EPS) or Portable Network Graphics format (PNG).

=head2 Available methods

This module provides the following methods:

    new()
    set_gnuplot_program($path_to_gnuplot)
    set_output_format($output_format)
    set_data_filename($filename)
    set_macro_filename($filename)
    set_output_filename($filename)
    set_xvar_description($name, $description)
    set_yvar_description($name, $description)
    set_zvar_description($name, $description)
    set_xvar_values($values_as_array_ref)
    set_yvar_values($values_as_array_ref)
    set_zvar_values($values_as_array_ref)
    set_xaxis_title($title)
    set_yaxis_title($title)
    set_graph_title($title)
    set_main_label($label)
    set_second_label($label)
    write_data_file()
    generate_macro()
    execute_macro()

=head2 Methods description

    new()
        Basic constructor allocating a Tifa::GnuPlotter object.
        Plot output format is by default set to "pdf".

    set_gnuplot_program($path_to_gnuplot)
        Sets the name of the gnuplot executable, including its full path.

    set_output_format($output_format)
        Sets the file format of the generated plots.
        Must be one of "pdf", "eps" or "png".

    set_data_filename($filename)
        Sets the filename of the gnuplot data file to be generated.

    set_macro_filename($filename)
        Sets the filename of the gnuplot macro file to be generated.

    set_output_filename($filename)
        Sets the filename of the gnuplot plot file to be generated.
        If the provided $filename does not end with the proprer file
        extension, the good file extension will be appended according
        to the set output format.

    set_xvar_description($name, $description)
        Sets the name of the x variable together with a short
        description to be used in the generated data file.

    set_yvar_description($name, $description)
        Sets the name of the y variable together with a short
        description to be used in the generated data file.

    set_zvar_description($name, $description)
        Sets the name of the z variable together with a short
        description to be used in the generated data file.

    set_xvar_values(@x_values)
        Sets the values of the x variables to the ones in the array
        referenced by $values_as_array_ref.

    set_yvar_values(@y_values)
        Sets the values of the x variables to the ones in the array
        referenced by $values_as_array_ref.

    set_zvar_values(@z_values)
        Sets the values of the x variables to the ones in the array
        referenced by $values_as_array_ref.

    set_xaxis_title($title)
        Sets the title of the x axis.

    set_yaxis_title($title)
        Sets the title of the y axis.

    set_graph_title($title)
        Sets the title of the generated graph.

    set_main_label($label)
        Sets the main label that will appear under the graph title.

    set_second_label($label)
        Sets the second label that will appear on the plot.

    write_data_file()
        Writes the gnuplot data file.

    generate_macro()
        Generates the gnuplot macro file.

    execute_macro()
        Executes the generated gnuplot macro file.

=head1 EXAMPLE

The following code snippet shows how to create a 2D plot:

 #
 # The 3-dimensional toy-data...
 #
 my @x_array = (1.0, 2.1, 3.1, 4.0, 1.1, 2.1, 3.0,  3.9);
 my @y_array = (2.1, 3.9, 5.9, 7.8, 3.2, 6.4, 9.4, 12.2);
 my @z_array = (2,   2,     2,   2,   3,   3,   3,   3);

 $plotter = new Tifa::GnuPlotter();
 #
 # Set the relevant parameters...
 #
 $plotter->set_xvar_description("time", "Time Elapsed (s)");
 $plotter->set_yvar_description("distance", "Travelled Distance (m)");
 $plotter->set_zvar_description("acceleration",
                                "Theoretical Acceleration");
 #
 # Of course, one should in principle check that all the data
 # arrays have the same length.
 #
 $plotter->set_xvar_values(\@x_array);
 $plotter->set_yvar_values(\@y_array);
 $plotter->set_zvar_values(\@z_array);

 $plotter->set_xaxis_title("Time (s)");
 $plotter->set_yaxis_title("Distance (m)");
 $plotter->set_graph_title("Travelled distance as a function of time");
 $plotter->set_main_label('Projectile: book\nWhen: after the final');
 $plotter->set_second_label('Measured accelerations match theory!');

 $plotter->set_data_filename("flying_book.dat");
 $plotter->set_macro_filename("flying_book.gp");
 $plotter->set_output_filename("flying_book.pdf");

 $plotter->write_data_file();
 $plotter->generate_macro();
 #
 # Executing the flying_book.gp macro will produce a PDF file
 # displaying two graphs drawn on the same plot, one for each
 # value of the "acceleration" parameter.
 #
 $plotter->execute_macro();

=head1 EXPORT

No functions are exported from this package by default.

=head1 SEE ALSO

The documentation for the Tifa::DataDescriptor and Tifa::FormatConverter
Perl modules.

The gnuplot documentation at:
    http://www.gnuplot.info/documentation.html

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




