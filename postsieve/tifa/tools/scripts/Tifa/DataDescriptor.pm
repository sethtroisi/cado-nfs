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

package Tifa::DataDescriptor;

#-------------------------------------------------------------------------------
#               Write and read data description/entries in a file.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : DataDescriptor.pm
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
# This Perl module is part of the TIFA (Tools for Integer FActorization)
# library. Its goal is to provide a consistent interface to document and
# write a text data file and to retrieve these documentation and data entries
# from a data file in order to process it.
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
sub new {
    my $type = shift(@_);
    my $self = {};

    $self->{descr_open_tag}     = "<DataDescription>";
    $self->{descr_close_tag}    = "</DataDescription>";

    $self->{format_open_tag}   = "<format>";
    $self->{format_close_tag}   = "</format>";

    $self->{field_open_tag}     = "<field>";
    $self->{field_close_tag}    = "</field>";

    $self->{comment_open_tag}   = undef;
    $self->{comment_close_tag}  = undef;

    $self->{data_separator_tag} = ":";

    $self->{field_names}        = [];
    $self->{field_descriptions} = ();

    $self->{entry_open_tag}     = "<data>";
    $self->{entry_close_tag}    = "</data>";

    bless $self, $type;

    $self->set_comment_style("Perl");

    return $self;
}
#-------------------------------------------------------------------------------
sub set_comment_style {
    my $self  = shift(@_);
    my $style = shift(@_);

    if ($style eq "Perl") {
        $self->{comment_open_tag}  = "#";
        $self->{comment_close_tag} = "";
    } elsif ($style eq "C") {
        $self->{comment_open_tag}  = "/*";
        $self->{comment_close_tag} = "*/";
    } elsif ($style eq "C++") {
        $self->{comment_open_tag}  = "//";
        $self->{comment_close_tag} = "";
    } else {
        croak("ERROR: DataDescriptor::set_comment_style:"
              ." Unknown style $style\n");
    }
}
#-------------------------------------------------------------------------------
sub set_field_separator() {
    my $self = shift(@_);
    my $sep  = shift(@_);

    $self->{data_separator_tag} = $sep;
}
#-------------------------------------------------------------------------------
sub get_field_separator() {
    my $self = shift(@_);

    return $self->{data_separator_tag};
}
#-------------------------------------------------------------------------------
sub add_field_description {
    my $self        = shift(@_);
    my $field_name  = shift(@_);
    my $field_descr = shift(@_);

    push(@{$self->{field_names}}, $field_name);
    $self->{field_descriptions}{$field_name} = $field_descr;
}
#-------------------------------------------------------------------------------
sub get_field_description {
    my $self       = shift(@_);
    my $field_name = shift(@_);

    return $self->{field_descriptions}{$field_name};
}
#-------------------------------------------------------------------------------
sub get_all_field_names {
    my $self = shift(@_);

    return (@{$self->{field_names}});
}
#-------------------------------------------------------------------------------
sub write_descriptions {
    my $self   = shift(@_);
    local *OUT = shift(@_);

    print OUT ( "$self->{comment_open_tag} ");
    print OUT ( "$self->{descr_open_tag}");
    print OUT ( " $self->{comment_close_tag}\n");

    print OUT ( "$self->{comment_open_tag} ");
    print OUT ( "$self->{comment_close_tag} \n");

    my $descline = join($self->{data_separator_tag}, @{$self->{field_names}});

    print OUT ( "$self->{comment_open_tag}    ");
    print OUT ( "$self->{format_open_tag} $descline");
    print OUT ( " $self->{format_close_tag} ");
    print OUT ( "$self->{comment_close_tag}\n");

    print OUT ( "$self->{comment_open_tag}");
    print OUT ( " $self->{comment_close_tag} \n");

    for my $field (@{$self->{field_names}}) {
        print OUT ( "$self->{comment_open_tag}    ");
        print OUT ( "$self->{field_open_tag} ");
        print OUT ( "$field $self->{data_separator_tag}");
        print OUT ( " $self->{field_descriptions}{$field} ");
        print OUT ( "$self->{field_close_tag}");
        print OUT ( " $self->{comment_close_tag}");
        print OUT ( "\n");
    }
    print OUT ( "$self->{comment_open_tag} ");
    print OUT ( "$self->{comment_close_tag} \n");

    print OUT ( "$self->{comment_open_tag} ");
    print OUT ( "$self->{descr_close_tag}");
    print OUT ( " $self->{comment_close_tag}\n");
}
#-------------------------------------------------------------------------------
sub load_descriptions {
    my $self     = shift(@_);
    local *IN    = shift(@_);

    $self->{field_names}        = [];
    $self->{field_descriptions} = ();

    my $o_com  = quotemeta($self->{comment_open_tag});
    my $c_com  = quotemeta($self->{comment_close_tag});
    my $o_desc = quotemeta($self->{descr_open_tag});
    my $c_desc = quotemeta($self->{descr_close_tag});
    my $o_fmt  = quotemeta($self->{format_open_tag});
    my $c_fmt  = quotemeta($self->{format_close_tag});
    my $o_fld  = quotemeta($self->{field_open_tag});
    my $c_fld  = quotemeta($self->{field_close_tag});
    my $d_sep  = quotemeta($self->{data_separator_tag});

    my $format_line = "";

    while (my $line = <IN>) {
        last if ($line =~ m/^\s*$o_com\s*$o_desc\s*$c_com/);
    }
    while (my $line = <IN>) {

        last if ($line =~ m/^\s*$o_com\s*$c_desc\s*$c_com/);
        next if ($line =~ m/^\s*$o_com\s*$c_com?$/);

        $line =~ s/\s*$//g;
        $line =~ s/^\s*//g;

        if ($line =~ m/^\s*$o_com\s*$o_fmt(.*)$c_fmt\s*$c_com/) {
            $format_line = $1;
            $format_line =~ s/\s*$//g;
            $format_line =~ s/^\s*//g;
            @{$self->{field_names}} = split(/$d_sep/, $format_line);
        }
        if ($line =~ m/^$o_com\s*$o_fld(.*)$c_fld\s*$c_com/) {

            my $desc_line = $1;
            $desc_line =~ s/\s*$//g;
            $desc_line =~ s/^\s*//g;

            if ($desc_line =~ m/^(.+)$d_sep(.+)/) {

                my @splitted  = split(/$d_sep/, $desc_line);

                my $fieldname = shift(@splitted);
                my $fielddesc = join("$d_sep", @splitted);

                $fieldname =~ s/\s*$//g;
                $fieldname =~ s/^\s*//g;

                $fielddesc =~ s/\s*$//g;
                $fielddesc =~ s/^\s*//g;

                $self->{field_descriptions}{$fieldname} = $fielddesc;
            }
        }
    }
}
#-------------------------------------------------------------------------------
sub write_data_entry() {
    my $self   = shift(@_);
    local *OUT = shift(@_);
    my @fields = @_;

    if ($#fields != $#{$self->{field_names}}) {
        croak("ERROR: DataDescriptor::write_data_entry: "
                 ."number of fields does not match with description.\n");
    }

    my $entry = join($self->{data_separator_tag}, @fields);

    print OUT ( "$self->{entry_open_tag} ");
    print OUT ( "$entry ");
    print OUT ( "$self->{entry_close_tag}\n");
}
#-------------------------------------------------------------------------------
sub read_next_data_entry() {
    my $self   = shift(@_);
    local *IN  = shift(@_);

    my $o_ent = quotemeta($self->{entry_open_tag});
    my $c_ent = quotemeta($self->{entry_close_tag});
    my $d_sep = quotemeta($self->{data_separator_tag});

    my %data_entry = ();

    while (my $line = <IN>) {
        if ($line =~ m/^\s*$o_ent(.*)$c_ent/) {

            my $entry = $1;
            $entry =~ s/^\s*//;
            $entry =~ s/\s*$//;

            my @values = split(/$d_sep/, $entry);

            if ($#values != $#{$self->{field_names}}) {
                croak("ERROR: DataDescriptor::write_data_entry: number of"
                      ."fields does not match with description.\n");
            }
            foreach my $i (0..$#{$self->{field_names}}) {
                my $field_name = ${$self->{field_names}}[$i];
                $data_entry{$field_name} = $values[$i];
            }
            return %data_entry;
        }
    }
    return %data_entry;
}
#-------------------------------------------------------------------------------

1;
__END__

#-------------------------------------------------------------------------------
# Stub documentation for the DataDescriptor module.
#-------------------------------------------------------------------------------

=head1 NAME

Tifa::DataDescriptor - Write and read data description/entries in a file.

=head1 SYNOPSIS

  use Tifa::DataDescriptor;
  $descriptor = new Tifa::DataDescriptor();

=head1 REQUIRE

Perl 5.006002, Carp and Exporter.

=head1 DESCRIPTION

This Perl module is part of the TIFA (Tools for Integer FActorization) library.
Its goal is to provide a consistent interface to document and write a text data
file and to retrieve these documentation and data entries from a data file
in order to process it.

=head2 Available methods

This module provides the following methods:

    Tifa::DataDescriptor()
    set_comment_style($style)
    set_field_separator($separator)
    get_field_separator()
    add_field_description($field_name, $description)
    get_field_description($field_name)
    get_all_field_names()
    write_descriptions($file_handle)
    load_descriptions($file_handle)
    write_data_entry($file_handle, @data_values)
    read_next_data_entry($file_handle)

=head2 Methods description

    Tifa::DataDescriptor()
        Basic constructor allocating a Tifa::DataDescriptor object.

    set_comment_style($style)
        Sets the comment style used for the data description.
        Allowed values: "Perl", "C" and "C++".
        Default value : "Perl".

    set_field_separator($separator)
        Sets the motif used to separate the data fields.
        Default value: ":".

    get_field_separator()
        Returns the motif used to separate the data fields.

    add_field_description($field_name, $description)
        Adds in the data description a field $field_name with its
        description given by $description.

    get_field_description($field_name)
        Gets the description of the field $field_name.

    get_all_field_names()
        Returns an array containing all of the data fields' names.

    write_descriptions($file_handle)
        Writes the data description to an already opened file given
        by $file_handle.

    load_descriptions($file_handle)
        Loads the data description from the opened file referenced
        by $file_handle.

    write_data_entry($file_handle, @data_values)
        Writes a data entry with values given by the array @data_values
        in the opened file referenced by $file_hanfle.

    read_next_data_entry($file_handle)
        Reads the next data entry in the file referenced by $file_handle
        and returns it as a hashtable mapping the field names to their
        values. If no entry is found, returns an empty hashtable.

=head1 EXAMPLES

=head2 Example 1: write data description to a file

This following code snippet gives an example of how to use the
Tifa::DataDescriptor module to embed data description in a file.

    my $descr    = new Tifa::DataDescriptor();
    my $filename = "os_mascots.txt";
    open(my $handle, ">$filename") or die("Cannot open $filename\n");
    #
    # Fill the descriptor with the data description.
    #
    $descr->set_comment_style("Perl");
    $descr->set_field_separator(":");
    $descr->add_field_description("Name", "Name of mascot");
    $descr->add_field_description("Species", "Species of mascot");
    $descr->add_field_description("FavOS", "Favorite OS of mascot");
    #
    # Write the data description in the data file.
    #
    $descr->write_descriptions($handle);
    #
    # Write some data entries in the data file.
    #
    my @data = ('Tux', 'penguin', 'Linux');
    $descr->write_data_entry($handle, @data);

    @data = ('Hexley', 'platypus', 'Darwin');
    $descr->write_data_entry($handle, @data);

    @data = ('Beastie', 'daemon', 'FreeBSD');
    $descr->write_data_entry($handle, @data);

    @data = ('Puffy', 'pufferfish', 'OpenBSD');
    $descr->write_data_entry($handle, @data);

    close(my $handle);


The following description and entries are written at the current position in
the file:

    # <DataDescription>
    #
    #    <format> Name:Species:FavOS </format>
    #
    #    <field> Name : Name of mascot </field>
    #    <field> Species : Species of mascot </field>
    #    <field> FavOS : Favorite OS of mascot </field>
    #
    # </DataDescription>

    <data> Tux:penguin:Linux </data>
    <data> Hexley:platypus:Darwin </data>
    <data> Beastie:daemon:FreeBSD </data>
    <data> Puffy:pufferfish:OpenBSD </data>

=head2 Example 2: read data description in a file

This following code snippet gives an example of how to use the
Tifa::DataDescriptor module to retrieve the data description and entries
embedded within a file.

    my $descr    = new Tifa::DataDescriptor();
    my $filename = "data.txt";
    open(my $handle, "<$filename") or die("Cannot open $filename\n");
    #
    # Assume that the data file contains comments given in a Perl style.
    #
    $descr->set_comment_style("Perl");
    #
    # Read the data description embedded in the data file.
    #
    $descr->load_descriptions($filename);
    #
    # Read all the data entries and print them on stdout.
    #
    while (%entry = $descr->read_next_data_entry($handle)) {

        foreach my $key (keys %entry) {
            print("\$entry{$key} = $entry{$key}\n");
        }
        print("\n");
    }

    close($handle);

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

