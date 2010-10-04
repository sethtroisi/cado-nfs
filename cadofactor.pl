#!/usr/bin/perl -w

# General script for factoring integers with Cado-NFS.
#
# Copyright 2008 Pierrick Gaudry, Emmanuel Thome, Paul Zimmermann,
#                Jeremie Detrey
#
# This file is part of CADO-NFS.
#
# CADO-NFS is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CADO-NFS; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301, USA.



# Usage:
#    cadofactor.pl param=<paramfile> wdir=<...> ...
# Parameters passed in arguments *after* param=... will override choices
# that are made in paramfile.
#
# See params.c59 for an example of parameter file.
#
# If the parameter n=<n> is given, then n is factored. Otherwise, it is
# taken from stdin.
#
# NB: all the shell commands that are run by this script are reported in
# the $wdir/$name.cmd file, with full arguments, so that it is easy to
# reproduce part of the computation.



# TODO-list:
#  - Bench polynomials with some sieve before selecting the best one.
#  - Enable a 'lowmem' option
#  - To be able to import the msieve relations in CADO-NFS.
#  - Use new Murphy value in polynomial selection (C code only).
#  - Print in each kjout the number of real roots of polynomial
#  - Automate the parameters setting.

use Cwd qw(abs_path);
use File::Basename;
use lib abs_path(dirname($0));
use cadofct;
use strict;
use warnings;

select(STDERR); $| = 1; # always flush stderr
select(STDOUT); $| = 1; # always flush stdout

print "$0 @ARGV\n";

do_init();
do_task("sqrt");

banner("All done!");

open FILE, "$param{'prefix'}.allfactors"
or die "Cannot open `$param{'prefix'}.allfactors' for reading: $!.\n";
my @l = <FILE>;
close FILE;
chomp for @l;
print "@l\n";

if (defined(my $e = $param{'expected_factorization'})) {
	my @exp=split(',',$e);
	@exp = sort @exp;
	@l = sort @l;
	my $ok = "@exp" eq "@l";
	if ($ok) {
		print "Factorization matches expected result\n";
	} else {
		die "Factorization does not match expected result\n";
	}
}

