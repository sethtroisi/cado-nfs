#!/usr/bin/perl
#  This file is part of the gf2x library.
# 
#  Copyright 2007, 2008, 2009, 2010
#  Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
# 
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 2.1 of the License, or (at
#  your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License for more details.
#  
#  You should have received a copy of the GNU Lesser General Public
#  License along with CADO-NFS; see the file COPYING.  If not, write to
#  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
#  Boston, MA 02110-1301, USA.


use warnings;
use strict;

my $fragment = '';
while (<>) {
    if (/^# -- begin generated code --$/) {
        print;
        print $fragment;
        while (<>) {
            if (/^# -- end generated code --$/) {
                print;
                last;
            }
        }
        next;
    }
    if (!/^#\*/) {
        print;
        next;
    }
    /^#\*\s*(.*)$/ or die;
    my $cmd = $1;
    $cmd =~ s/(\S*)\{(\d+)\.\.(\d+)\}(\S*)/join(" ", map { "$1$_$4" } ($2..$3))/gex;

    if ($cmd =~ /^test_sizes ([\d\s]+)$/) {
        my @sizes = split(' ', $1);
        for my $x (@sizes) {
            $fragment .= <<EOF;
tuning_undefs_$x.h: tuning_undefs_pre.h ; sed -e s/\@\@SIZE\@\@/$x/g < \$< > \$\@
tuneup_$x.c: \$(srcdir)/tuneup_pre.c ; sed -e s/\@\@SIZE\@\@/$x/g < \$< > \$\@
tuneup_$x.\$(OBJEXT): tuning_undefs_$x.h
BUILT_SOURCES +=tuneup_$x.c tuning_undefs_$x.h
EXTRA_LTLIBRARIES+=libtuneup-s$x.la
nodist_libtuneup_s${x}_la_SOURCES=tuneup_$x.c tuning_undefs_$x.h
libtuneup_s${x}_la_CPPFLAGS=-I\$(top_builddir) -I\$(top_srcdir)
CLEANFILES+=tuning_undefs_$x.h
DISTCLEANFILES+=tuning_undefs_$x.h
EOF
            $fragment .= "\n";
        }
    } elsif ($cmd =~ /^test\s+(.*)/) {
        my $funcs = $1;
        for my $f (split(' ', $funcs)) {
            $f =~ /^((?:gen_)?mul(\d+)\w+)$/ or die;
            my $code=$1;
            my $size=$2;
            if ($code =~ /^gen/) {
                $code =~ /w(\d+)k(\d+)/ or die;
                my $w = $1;
                my $k = $2;
                $fragment .= <<EOF;
LOWLEVEL+=tune_$code\$(EXEEXT)
tune_${code}_CPPFLAGS=-I\$(top_builddir) -I\$(top_srcdir) -DTUNING=$size
nodist_tune_${code}_SOURCES=$code.c tuning_undefs_$size.h
$code.\$(OBJEXT): tuning_undefs_$size.h
tune_${code}_LDADD=libtuneup-s$size.la libtiming.la ../libgf2x.la
$code.c: gen_bb_mul_code\$(EXEEXT) ; ./gen_bb_mul_code $w $k > \$@
BUILT_SOURCES+=$code.c
EOF
            } else {
                $fragment .= <<EOF;
LOWLEVEL+=tune_${code}\$(EXEEXT)
tune_${code}_CPPFLAGS=-I\$(top_builddir) -I\$(top_srcdir) -DTUNING=$size
tune_${code}_SOURCES=$code.c
nodist_tune_${code}_SOURCES=tuning_undefs_$size.h
tune_${code}_LDADD=libtuneup-s$size.la libtiming.la ../libgf2x.la
$code.\$(OBJEXT): tuning_undefs_$size.h
EOF
            }
            $fragment .= "\n";
            $fragment .= "\n";
        }
    } else {
        $fragment .= "$cmd\n";
    }
    print;
}
