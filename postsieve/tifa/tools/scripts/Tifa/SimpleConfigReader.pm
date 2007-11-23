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

package Tifa::SimpleConfigReader;

#-------------------------------------------------------------------------------
#                  Truly minimalist configuration file reader.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : SimpleConfigReader.pm
# Author        : Jerome Milan
# Created on    : circa late July / early August 2006
# Last modified : Wed Jan 31 2007
#
# Version : 0.1.2
# License : GNU Lesser General Public License (LGPL) v2.1 or later
#           Copyright (C) 2006, 2007 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
# 0.1.2: Wed Jan 31 2007 by JM
#       - Added $accept_unknown_param field and the methods
#         accept_unknown_param() and dont_accept_unknown_param().
#       - Multiple line declarations allowed (line must end with ' \')
#       - Now accept hashtable parameters
#       - Some cleanup.
#
# 0.1.1: Tue Oct 4 2006 by JM
#       - Added $warnings_are_fatal field and the methods use_fatal_warnings()
#         and use_non_fatal_warnings().
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# This is a truly minimalist configuration file reader used by the
# scripts of the TIFA library. As such, the needs it fits are precisely
# defined and consequently its feature set and tolerance to syntax errors
# are very limited. This is certainly not a general purpose configuration
# file parser, so its use outside of the TIFA libray is really discouraged.
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

    $self->{comment_open_tag}   = '#';

    $self->{begin_scalar_tag} = '';
    $self->{begin_list_tag}   = 'list';
    $self->{begin_range_tag}  = 'range';
    $self->{begin_hash_tag}   = 'hash';

    $self->{parameters}         = ();
    $self->{scalar_param_names} = ();
    $self->{array_param_names}  = ();
    $self->{hash_param_names}   = ();

    $self->{warnings_are_fatal} = 1;

    $self->{accept_unknown_param} = 0;

    bless $self, $type;

    return $self;
}
#-------------------------------------------------------------------------------
sub use_fatal_warnings {
    my $self = shift;

    $self->{warnings_are_fatal} = 1;
}
#-------------------------------------------------------------------------------
sub use_non_fatal_warnings {
    my $self = shift;

    $self->{warnings_are_fatal} = 0;
}
#-------------------------------------------------------------------------------
sub accept_unknown_param {
    my $self = shift;

    $self->{accept_unknown_param} = 1;
}
#-------------------------------------------------------------------------------
sub dont_accept_unknown_param {
    my $self = shift;

    $self->{accept_unknown_param} = 0;
}
#-------------------------------------------------------------------------------
sub print_all_param_values {
    my $self = shift;

    my %param_hash    = %{$self->{parameters}};

    foreach my $scalar_params (keys %{$self->{scalar_param_names}}) {
        print("$scalar_params =  $param_hash{$scalar_params}\n");
    }
    foreach my $array_params (keys %{$self->{array_param_names}}) {
        print("$array_params = {\n");
        foreach my $array_vals (@{$param_hash{$array_params}}) {
            print("    $array_vals\n");
        }
        print("}\n");
    }
    foreach my $hash_params (keys %{$self->{hash_param_names}}) {
        print("$hash_params = {\n");
        foreach my $hash_key (keys %{$param_hash{$hash_params}}) {
            print("    $hash_key => ".${$param_hash{$hash_params}}{$hash_key});
            print("\n");
        }
        print("}\n");
    }
}
#-------------------------------------------------------------------------------
sub get_parameter_hash {
    my $self = shift;

    return %{$self->{parameters}};
}
#-------------------------------------------------------------------------------
sub set_scalar_param_names {
     my $self   = shift(@_);

     my $error = "ERROR: SimpleConfigReader::set_scalar_param_names:";

     undef %{$self->{scalar_param_names}};

     foreach my $name (@_) {
         $name = scalar $name;

         if ($name !~ m/^[a-zA-Z]\w*/) {
             croak("$error Invalid parameter name syntax $name.\n");
         }
         if (defined ${$self->{array_param_names}}{$name}) {
             croak("$error $name is already defined as an array.\n");
         }
         if (defined ${$self->{hash_param_names}}{$name}) {
              croak("$error $name is already defined as a hashtable.\n");
         }
         ${$self->{scalar_param_names}}{$name} = 1;
     }
}
#-------------------------------------------------------------------------------
sub set_array_param_names {
    my $self   = shift(@_);

    my $error = "ERROR: SimpleConfigReader::set_array_param_names:";

    undef %{$self->{array_param_names}};

    foreach my $name (@_) {
        $name = scalar $name;

        if ($name !~ m/^[a-zA-Z]\w*/) {
             croak("$error Invalid parameter name syntax $name.\n");
        }
        if (defined ${$self->{scalar_param_names}}{$name}) {
             croak("$error $name is already defined as a scalar.\n");
        }
        if (defined ${$self->{hash_param_names}}{$name}) {
             croak("$error $name is already defined as a hashtable.\n");
        }
        @{$self->{array_param_names}}{$name} = 1;
    }
}
#-------------------------------------------------------------------------------
sub set_hash_param_names {
    my $self   = shift(@_);

    my $error = "ERROR: SimpleConfigReader::set_hash_param_names:";

    undef %{$self->{hash_param_names}};

    foreach my $name (@_) {
        $name = scalar $name;

        if ($name !~ m/^[a-zA-Z]\w*/) {
             croak("$error Invalid parameter name syntax $name.\n");
        }
        if (defined ${$self->{scalar_param_names}}{$name}) {
             croak("$error $name is already defined as a scalar.\n");
        }
        if (defined ${$self->{array_param_names}}{$name}) {
             croak("$error $name is already defined as a array.\n");
        }
        @{$self->{hash_param_names}}{$name} = 1;
    }
}
#-------------------------------------------------------------------------------
sub get_scalar_param_names {
    my $self   = shift(@_);

    return (keys %{$self->{scalar_param_names}});
}
#-------------------------------------------------------------------------------
sub get_array_param_names {
    my $self   = shift(@_);

    return (keys %{$self->{array_param_names}});
}
#-------------------------------------------------------------------------------
sub get_hash_param_names {
    my $self   = shift(@_);

    return (keys %{$self->{hash_param_names}});
}
#-------------------------------------------------------------------------------
sub read_config_file {

    my $self     = shift(@_);
    my $filename = shift(@_);

    my $error   = "ERROR: SimpleConfigReader::parse_config_file:";
    my $warning = "WARNING: SimpleConfigReader::parse_config_file:";

    open(CF, "<$filename") or croak("$error Cannot open $filename.\n");

    my $cmnt_tag = $self->{comment_open_tag};
    my $sclr_tag = $self->{begin_scalar_tag};
    my $list_tag = $self->{begin_list_tag};
    my $rng_tag  = $self->{begin_range_tag};
    my $hash_tag = $self->{begin_hash_tag};

    my $name_re = '\s*([a-zA-Z_]\w*)\s*';
    my $eq_re   = '\s*=\s*';
    my $sclr_re = '\s*'.quotemeta($sclr_tag).'\s*(.*)\s*';
    my $list_re = '\s*'.quotemeta($list_tag).'\s*\((.*)\)\s*';
    my $rng_re  = '\s*'.quotemeta($rng_tag).'\s*\((.*)\)\s*';
    my $hash_re = '\s*'.quotemeta($hash_tag).'\s*\((.*)\)\s*';

    my $rng_cnt_re  = '\s*from\s*([+-]?\d+\.?\d*)\s*to\s*([+-]?\d+\.?\d*)';
    $rng_cnt_re    .= '\s*increment\s*([+-]?\d+\.?\d*)\s*';

    my $hash_cnt_re = '(\s*\w*\s*=>\s*\w*\s*,)*\s*\w*\s*=>\s*\w*\s*';

    my $info = "";

    while (my $line = <CF>) {

        next if ($line =~ m/^\s*$cmnt_tag/);
        next if ($line =~ m/^\s*$/);

        $line =~ s/\s*$//g;
        $line =~ s/^\s*//g;

        if ($line =~ m/\ \\$/) {
            #
            # Multiple line declarations
            #
            my $nextline = $line;

            $line =~ s/\ \\$//;
            $line =~ s/\s+/\ /;

            while ($nextline =~ m/\ \\$/) {

                if (eof(CF)) {
                    croak("$error end of file $filename reached before end",
                          " of line!\n");
                }
                $nextline = <CF>;

                $nextline =~ s/\s*$//g;

                if (    ($nextline =~ m/^\s*$cmnt_tag/)
                     || ($nextline =~ m/^\s*$/))  {

                    $nextline = " \\";
                    next;
                }
                $line .= $nextline;

                $line =~ s/\ \\$//;
                $line =~ s/\s+/\ /g;
            }
        }
        $info = "line $. in file $filename";

        #
        # Matches hashtables
        #
        if ($line =~ m/$name_re$eq_re$hash_re/) {

            my $name    = $1;
            my $content = $2;

            if ($content =~ m/$hash_cnt_re/) {
                my @keyvals = split(/\s*,\s*/, $content);
                my @keys    = ();
                my @vals    = ();

                foreach my $keyval (@keyvals) {
                    (my $key, my $val) = split(/\s*=>\s*/, $keyval, 2);
                    $key =~ s/^\s*//;
                    $key =~ s/\s*$//;
                    $val =~ s/^\s*//;
                    $val =~ s/\s*$//;

                    push(@keys, $key);
                    push(@vals, $val);
                }
                $self->__add_hash_param__($info, $name, \@keys, \@vals);
                next;

            } else {
                croak("$error Invalid syntax at line $. in file $filename\n");
            }
        }
        #
        # Matches lists
        #
        if ($line =~ m/$name_re$eq_re$list_re/) {
            my $name    = $1;
            my $content = $2;
            $self->__add_array_param__($info, $name,
                                       split(/\s*,\s*/, $content));
            next;
        }
        #
        # Matches ranges
        #
        if ($line =~ m/$name_re$eq_re$rng_re/) {
            my $name    = $1;
            my $content = $2;

            if ($content =~ m/$rng_cnt_re/) {
                my $from = $1;
                my $to   = $2;
                my $incr = $3;
                my @data = ();
                for (my $i = $from; $i <= $to; $i += $incr) {
                    push(@data, $i);
                }
                $self->__add_array_param__($info, $name, @data);
            } else {
                croak("$error Invalid syntax at line $. in file $filename\n");
            }
            next;
        }
        #
        # Matches scalars
        #
        if ($line =~ m/$name_re$eq_re$sclr_re/) {
            my $name    = $1;
            my $content = $2;
            $self->__add_scalar_param__($info, $name, $content);
            next;
        }
    }
    close(CF);
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  "Private" functions
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
sub __add_scalar_param__ {
    my $self  = shift;
    my $info  = shift;
    my $name  = shift;
    my $value = shift;

    my $is_defined = eval { exists ${$self->{scalar_param_names}}{$name} };
    my $add_param  = eval { $self->{accept_unknown_param} || $is_defined };

    my $error = "ERROR: SimpleConfigReader::__add_scalar_param__:";

    if (defined ${$self->{array_param_names}}{$name}) {
        croak("$error $name is already defined as an array.\n");
    }
    if (defined ${$self->{hash_param_names}}{$name}) {
         croak("$error $name is already defined as a hashtable.\n");
    }

    if ($add_param) {
        if ($name !~ m/^[a-zA-Z_]\w*/) {
            croak("$error Invalid parameter name syntax $name.\n");
        }
        ${$self->{scalar_param_names}}{$name} = 1;
        ${$self->{parameters}}{$name}         = $value;

    } else {
        my $warning = "ERROR: SimpleConfigReader::__add_scalar_param__:";
        my $mess    =  "$name is an unknown scalar parameter name at $info."
                      ." Possible typo?\n";
        if ($self->{warnings_are_fatal}) {
            croak("$error $mess");
        } else {
            carp("$warning $mess");
        }
    }
}
#-------------------------------------------------------------------------------
sub __add_array_param__ {
    my $self  = shift;
    my $info  = shift;
    my $name  = shift;
    my @value = @_;

    my $is_defined = eval { exists ${$self->{array_param_names}}{$name} };
    my $add_param  = eval { $self->{accept_unknown_param} || $is_defined };

    my $error = "ERROR: SimpleConfigReader::__add_array_param__:";

    if (defined ${$self->{scalar_param_names}}{$name}) {
        croak("$error $name is already defined as a scalar.\n");
    }
    if (defined ${$self->{hash_param_names}}{$name}) {
         croak("$error $name is already defined as a hashtable.\n");
    }
    if ($add_param) {

        if ($name !~ m/^[a-zA-Z_]\w*/) {
            croak("$error Invalid parameter name syntax $name.\n");
        }
        ${$self->{array_param_names}}{$name} = 1;
        if (! @{$self->{parameters}}{$name}) {
            @{$self->{parameters}{$name}} = ();
        }
        push(@{$self->{parameters}{$name}}, @value);

    } else {
        my $warning = "ERROR: SimpleConfigReader::__add_array_param__:";
        my $mess    =  "$name is an unknown array parameter name at $info."
                      ." Possible typo?\n";
        if ($self->{warnings_are_fatal}) {
            croak("$error $mess");
        } else {
            carp("$warning $mess");
        }
    }
}
#-------------------------------------------------------------------------------
sub __add_hash_param__ {
    my $self   = shift;
    my $info   = shift;
    my $name   = shift;
    my $keysref   = shift;
    my $valuesref = shift;

    my @keys   = @$keysref;
    my @values = @$valuesref;

    my $is_defined = eval { exists ${$self->{hash_param_names}}{$name} };
    my $add_param  = eval { $self->{accept_unknown_param} || $is_defined };

    my $error = "ERROR: SimpleConfigReader::__add_hash_param__:";

    if (defined ${$self->{scalar_param_names}}{$name}) {
        croak("$error $name is already defined as a scalar.\n");
    }
    if (defined ${$self->{array_param_names}}{$name}) {
         croak("$error $name is already defined as an array.\n");
    }
    if ($add_param) {
        if ($name !~ m/^[a-zA-Z_]\w*/) {
            croak("$error Invalid parameter name syntax $name.\n");
        }
        ${$self->{hash_param_names}}{$name} = 1;
        if (exists $self->{parameters}{$name}) {
            %{$self->{parameters}{$name}} = ();
        }
        while (@keys) {
            ${$self->{parameters}{$name}}{shift(@keys)} = shift(@values);
        }

    } else {
        my $warning = "ERROR: SimpleConfigReader::__add_hash_param__:";
        my $mess    =  "$name is an unknown hash parameter name at $info."
                      ." Possible typo?\n";
        if ($self->{warnings_are_fatal}) {
            croak("$error $mess");
        } else {
            carp("$warning $mess");
        }
    }
}
#-------------------------------------------------------------------------------

1;

__END__


#-------------------------------------------------------------------------------
# Stub documentation for this module.
#-------------------------------------------------------------------------------

=head1 NAME

Tifa::SimpleConfigReader - Truly minimalist configuration file reader

=head1 SYNOPSIS

  use Tifa::SimpleConfigReader;
  $reader = new Tifa::SimpleConfigReader();

=head1 REQUIRE

Perl 5.006002, Carp and Exporter.

=head1 SUMMARY

This is a truly minimalist configuration file reader used by the scripts of the
TIFA library. As such, the needs it fits are precisely defined and consequently
its feature set and tolerance to syntax errors are very limited. This is
certainly not a general purpose configuration file parser, so its use outside
of the TIFA libray is really discouraged.

=head1 DESCRIPTION

Full-featured, general purpose configuration file parsers abound on the
Comprehensive Perl Archive Network, from the most simple readers (e.g.
ConfigReader::Simple) to the most sophiticated parsers (e.g. Config::Scoped)
understanding blocks and interpolations. The objective of the
Tifa::SimpleConfigReader module is extremely modest: this module is not
trying to re-invent (a worse version of) the wheel but to fit a very specific
need encountered during the development of the TIFA library. So, yes, its
feature set is very restricted and next to no syntaxic sugar is allowed in
the configuration files (welcome back to the 70s!).

=head2 Configuration file format

Comments beginning by a '#' are allowed in the configuration files provided
that they stand on their own line.

Three kinds of parameters are allowed: scalars, arrays and hashtables.

Scalar parameter values are given using the most basic notation:

    <scalar_param> = <value>

Array parameter values can be given either by providing a list, or by providing
a range. A list is introduced by the "list" keyword followed by a
comma-delimited list of values in brackets:

    <array_param> = list(<value1>, <value2>, <value3>)

A range is described by the "range" keyword and its three values given by the
"from", "to" and "increment" keywords:

    <array_param> = range(from <value1> to <value2> increment <value3>)

Finally, hashtable parameters are introduced by the "hash" keyword:

    <hash_param> = hash(<key1> => <value1>, ..., <keyN> => <valueN>)

Directives can span multiples lines, as long as line breaks are explicitely
marked by ' \' (a space followed by a backslash). For example, the following
code will be interpreted correctly:

    #
    # This is correct: line breaks are explicitely indicated
    #
    my_important_parameter = range(             \
                                from 10         \
                                to   120        \
                                increment 5     \
                             )

However, the following snippet will produce a syntax error:

    #
    # This is incorrect: line breaks are not explicitely indicated
    #
    my_important_parameter = range(
                                from 10
                                to   120
                                increment 5
                             )

Another restriction lies in the format of the parameter names: they should
imperatively begin by a letter or by an underscore, although digits are
allowed in the other positions.

Generally speaking, Tifa::SimpleConfigReader is quite leniant with whitespaces,
be it in the beginning of a line, before of after the = sign, in a list of
values or in a range definition.

If several directives for the same scalar parameter are given in the
configuration file, previous values are discarded. However, if several
directives for the same array or hash parameter are given, all the values will
be kept in the array or hash. For example, let's consider the following
configuration snippet:

    location = Stanford
    location = Hamburg
    location = Geneva
    #
    # The final value for the "location" key is Geneva.
    #
    years = list(1985, 1995)
    years = range(from 2010 to 2016 increment 2)
    #
    # The final "years" array will contain 1985, 1995, 2010, 2012, 2014
    # and 2016.
    #
    positions = hash(CMS   => Octant5)
    positions = hash(Atlas => Octant7)
    positions = hash(Atlas => Octant1)
    #
    # The final "positions" hashtable will contain the mapping
    # CMS => Octant5 and Atlas => Octant1.
    #

=head2 Available methods

This module provides the following methods:

    new()
    get_parameter_hash()
    set_scalar_param_names(@names)
    set_array_param_names(@names)
    set_hash_param_names(@names)
    get_scalar_param_names()
    get_array_param_names()
    get_hash_param_names()
    print_all_param_values()
    read_config_file($filename);
    accept_unknown_param()
    dont_accept_unknown_param()
    use_fatal_warnings()
    use_non_fatal_warnings()

=head2 Methods description

    new()
        Basic constructor allocating a Tifa::SimpleConfigReader object.

    get_parameter_hash()
        Returns a hashtable mapping the names of the parameters read in
        the configuration file to their respective values.

    set_scalar_param_names(@names)
        Sets the names of the parameters supposed to be scalars. Each
        entry in the @names array should be the name of a parameter
        expected to be read in the configuration file.

    set_array_param_names(@names)
        Sets the names of the parameters supposed to be arrays. Each
        entry in the @names array should be the name of a parameter
        expected to be read in the configuration file.

    set_hash_param_names(@names)
        Sets the names of the parameters supposed to be hashtables. Each
        entry in the @names array should be the name of a parameter
        expected to be read in the configuration file.

    get_scalar_param_names()
        Returns the names of the scalar parameters as an array.

    get_array_param_names()
        Returns the names of the array parameters as an array.

    get_hash_param_names()
        Returns the names of the hash parameters as an array.

    print_all_param_values()
        Prints all of the parameter names and their respective values.

    read_config_file($filename);
        Reads a configuration file $filename in the current context and
        stores the values of the parameters read.

    accept_unknown_param()
        Accept unknown parameter names to appear in the configuration
        file and treat them as if their were declared via the
        set_*_param_names() functions. This can be used if the
        parameters likely to appear in the configuration file are not
        known in advance.

    dont_accept_unknown_param()
        Do not allowed undeclared parameters names for a strict check
        of the configuration file syntax.
        This is the default behaviour.

    use_fatal_warnings()
        Exit if an unknown parameter name if found. (This setting has no
        effect if unknown parameter names were explicitely allowed by
        calling the accept_unknown_param() function.)
        This is the default behaviour.

    use_non_fatal_warnings()
        Warn if an unknown parameter name if found and ignore it, but
        do not exit. (This setting has no effect if unknown parameter
        names were explicitely allowed by calling the
        accept_unknown_param() function.)

=head1 EXAMPLE

This following code snippet gives an example of how to use the
Tifa::SimpleConfigReader module to parse a simple configuration file. Let the
configuration file be given by:

 #
 # Write some meaningful comments here...
 #
 accelerator = Large Hadron Collider
 location    = Geneva

 #
 # Now, some array definitions...
 #
 experiments = list(Alice, Atlas, CMS, LHCb, TOTEM)
 bunch_nb    = range(from 0 to 2834 increment 1)
 objectives  = list(everything)

 #
 # And finally, a hashtable
 #
 positions = hash(Atlas => Octant1, CMS => Octant5)

The following code reads this configuration file and prints the values of the
parameters on the standard output:

 $reader = new Tifa::SimpleConfigReader();
 #
 # Inform the SimpleConfigReader object about the parameters defined
 # in the configuration file. Be warned that any parameter not
 # specified via the set_*_param_names methods could be treated as a
 # syntax error if the function accept_unknown_param() has not been
 # called.
 #
 @scalar_names = ("accelerator", "location");
 @array_names = ("experiments", "bunch_nb", "objectives");

 $reader->set_scalar_param_names(@scalar_names);
 $reader->set_array_param_names(@array_names);

 $filename = "config.txt";

 $reader->read_config_file($filename);
 $reader->print_all_param_values();
 #
 # Parameter values are now accessible via the get_parameter_hash()
 # method
 #
 %all_params = $reader->get_parameter_hash();
 $where      = $all_params{"location"};
 @detectors  = @{$all_params{"experiments"}};

 print("$detectors[1] is located near $where.\n");

The (truncated) output is:

    location =  Geneva
    accelerator =  Large Hadron Collider
    bunch_nb = {
        0
        1
        2
        [...]
        2833
        2834
    }
    experiments = {
        Alice
        Atlas
        CMS
        LHCb
        TOTEM
    }
    objectives = {
        everything
    }
    positions = {
        CMS => Octant5
        Atlas => Octant1
    }
    Atlas is located near Geneva.

=head1 EXPORT

No functions are exported from this package by default.

=head1 SEE ALSO

=head1 AUTHOR

Jerome Milan, E<lt>milanj@lix.polytechnique.frE<gt>

=head1 VERSION

0.1.2

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
