#!/usr/bin/perl -w
#
# General script for factoring integers with Cado-NFS.
#
# Copyright 2008, 2009, 2010 Pierrick Gaudry, Emmanuel Thome, Paul Zimmermann,
#                            Jeremie Detrey, Lionel Muller
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
#
#
# Usage:
# ======
#
# cadofactor.pl param=<paramfile> wdir=<...> ...
#
# Parameters passed in arguments *after* param=... will override choices
# that are made in paramfile.
#
# See params/params.c91 for an example of parameter file.
#
# If the parameter n=<n> is given, then n is factored. Otherwise, it is
# taken from stdin.
#
# NB: all the shell commands that are run by this script are reported in
# the $wdir/$name.cmd file, with full arguments, so that it is easy to
# reproduce part of the computation.
#
# cadofactor.pl is resistant to crash of the host or remote machines.
#  
# For a small factorization, you can use the factor.sh script (easier
# to use) who run the factorization only on your machine.
#
#
# Example for factorize an integer of 155 digits (on several machines):
# ====================================================================
# 
# Before to start the factorization, you must have to configure ssh (cf
# README) for all the machines.
#
# $mkdir $HOME/c155
# $cd $HOME/c155
# $cp $CADO_DIR/params/mach_desc .
#   Edit the file mach_desc and configure it.
# $cp $CADO_DIR/params/params.c155 .
#   It's not require to edit the paramfile but you can choose to modify
#   some parameters in this file or on line command of cadofactor.pl.
#   If the paramfile does not exist then you must to create it.
# $CADO_DIR/cadofactor.pl params.c155 wdir=$HOME/c155 name=rsa155 \
#    machines=mach_desc \
#    n=10941738641570527421809707322040357612003732945449205990913842131476349984288934784717997257891267332497625752899781833797076537244027146743531593354333897
#
#
# To factorize an integer with SNFS (on several machines):
# ==============================================================================
# 
# Before to start the factorization, you must have to configure ssh (cf
# README) for all the machines.
#
# $mkdir $HOME/snfs
# $cd $HOME/snfs
# $cp $CADO_DIR/params/mach_desc .
#   Edit the file mach_desc and configure it.
# $cp $CADO_DIR/params/params.c<size> .  
#   [substitute <size> by the integer size]
#   It's not require to edit the paramfile but you can choose to modify
#   some parameters in this file or on line command of cadofactor.pl.
#   If the paramfile does not exist then you must to create it.
# $cp <polynomial choosen in the cado format> snfs<size>.poly
#   The polynomial file must also contain the sieve parameters.
# $CADO_DIR/cadofactor.pl params.c<size> wdir=$HOME/snfs name=snfs<size> \
#    machines=mach_desc n=<n>
# $touch snfs<size>.polysel_done
# $CADO_DIR/cadofactor.pl params.c<size> wdir=$HOME/snfs name=snfs<size> \
#    machines=mach_desc n=<n>

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

