#!/usr/bin/env perl

# This script makes sure that for all compilation units found in the
# source tree, we have '#include "cado.h"' near the top, as we should.
#
# Of course we have exceptions. See @path_exceptions below.
#       gf2x/
#       linalg/bwc/mpfq/
#       linalg/bwc/flint-fft/
# This script must be called from the top of the tree.

use strict;
use warnings;

die "Please call $0 from the top of the source tree" unless -f "cado.h";

my @path_exceptions=qw|
        gf2x/
        linalg/bwc/mpfq/
        linalg/bwc/flint-fft/
        config/
        misc/
        |;


sub add_file {
    my @res = ();
    open my $fh, shift or return;
    while (defined($_=<$fh>)) {
        next unless /^[^#]/;
        chomp($_);
        next if m{/$};
        push @res, $_;
    }
    close $fh;
    return @res;
}

my @all_files;
push @all_files, &add_file('files.dist');
push @all_files, &add_file('files.nodist');

my @errors;
FILE: for my $f (@all_files) {
    # Check whether this matches in any way.
    $f =~ m{\.c(?:pp)?$} or next;
    for my $p (@path_exceptions) {
        next FILE if $f =~ /^$p/;
    }
    open my $fh, $f or die "$f: $!";
    while (defined($_=<$fh>)) {
        if (/\#include/) {
            push @errors, $f unless /"cado\.h"/;
            last;
        }
    }
    close $fh;
}
print "$_\n" for @errors;
if (@errors) {
    my $n = scalar @errors;
    print STDERR "$n errors found. Please fix\n";
    exit 1;
}
