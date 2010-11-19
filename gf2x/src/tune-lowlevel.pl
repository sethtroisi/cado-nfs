#!/usr/bin/perl -w
# This file is part of the gf2x library.
#
# Copyright 2007, 2008, 2009, 2010
# Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 2.1 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with CADO-NFS; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA 02110-1301, USA.


use warnings;
use strict;

sub usage {
    die "Usage: ./tune-lowlevel.pl <list of test programs>\n";
}

my %sizes=();

for my $x (@ARGV) {
    $x =~ /^tune_/ or die "Bad test program $x";
    $x =~ /mul(\d)/ or die "Bad test program $x";
    my $s=$1;
    if (!exists $sizes{$s}) { $sizes{$s}=[]; }
    push @{$sizes{$s}}, $x;
}

sub mysys {
    print STDERR "@_\n";
    system @_;
}



my $make;
if (!defined($make=$ENV{'MAKE'})) {
    $make="make";
}

# make sure the upper build is complete. It might be bad to do this
# check, but since this script is bound to call make anyway... 
mysys "cd .. ; $make";

my @summary = ();

for my $s (sort { $a <=> $b } keys %sizes) {
    print STDERR "Benching for $s words\n";
    # Now we check everything, always.
#    if (scalar @{$sizes{$s}} == 1 && !defined($ENV{'BENCH'})) {
#        (my $x = $sizes{$s}->[0]) =~ s/^tune_//;
#        print STDERR "Only one possibility ($x) -- check skipped\n";
#        push @summary, "mul$s -> $x.c (only choice)\n";
#        next;
#    }
    my @results;
    for my $p (@{$sizes{$s}}) {
        mysys "make $p";
        my $r = `./$p`;
        chomp($r);
        print STDERR "$r\n";
        $r =~ /^(.*)\s:\s([\d\.]+)\sns$/;
        push @results, [$1, $2];
    }
    @results = sort { $a->[1] <=> $b->[1] } @results;
    my $best = $results[0];
    $best->[0] =~ /tune_(.*)$/ or die;
    my $selected="$1.c";
    my $msg = "mul$s -> $selected [ $best->[1] ns ]";
    if (! -f $selected) {
        my $e;
        if (defined($e=$ENV{'srcdir'}) && -f "$e/$selected") {
            $selected="$e/$selected";
        } else {
            die "Cannot find $selected anywhere !"
        }
    }
    print STDERR "Selected $selected\n";
    my $ltarget="already_tuned/tuned/gf2x_mul$s.h";
    my $slot="gf2x/gf2x_mul$s.h";
    # my $prepared="ready_gf2x_mul$s.c";
    # mysys "sed -e s///g $selected > $prepared";
    my $rc=system "diff ../$slot $selected > /dev/null";

    if ($rc == 0) {
        print STDERR "Choice identical to the selected file\n";
        push @summary, "$msg (unchanged)\n";
        # unlink $prepared;
        next;
    }
    push @summary, "$msg\n";
    mkdir "../already_tuned" unless -d "../already_tuned";
    mkdir "../already_tuned/tuned" unless -d "../already_tuned/tuned";
    # Show the commands at the same time as we execute them.
    mysys "rm ../$slot";
    mysys "rm ../$ltarget";
    if ($selected =~ /gen_/) {
        # generated file: do a copy, not a link.
        mysys "cp -f $selected ../$ltarget";
    } else {
        mysys "ln -sf ../../src/$selected ../$ltarget";
    }
    mysys "ln -sf ../$ltarget ../$slot";

    print STDERR "Library source has changed -- rebuilding\n";
    mysys "cd .. ; $make";
}

print STDERR "Summary of tune-lowlevel results\n";
for my $m (@summary) {
    print STDERR $m;
}
