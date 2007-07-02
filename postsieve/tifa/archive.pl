#!/usr/bin/perl -w

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

#-------------------------------------------------------------------------------
#              Clean problematic .tgz archives created by SCons
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General Information
#-------------------------------------------------------------------------------
# File          : archive.pl
# Author        : Jerome Milan
# Created on    : Sometimes circa January (Februray?) 2007
# Last modified : Mon Mar 5 2007
#
# Version : 0.1.0
# License : GNU Lesser General Public License (LGPL)
#           Copyright (C) 2006, 2007 INRIA
#-------------------------------------------------------------------------------
# History
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Short description
#-------------------------------------------------------------------------------
# This script opens a .tgz file archive, cleans its content and repackages it
# in a new archive. Also, outputs the MD5 and SHA1 hashes of the produced
# archive.
#-------------------------------------------------------------------------------

use strict;
use Getopt::Long;
use File::Basename;
use File::Copy;
use File::Find;

my @times = localtime();
my $day   = sprintf("%02d", $times[3]);
my $month = sprintf("%02d", $times[4] + 1);
my $year  = $times[5] + 1900;

my $input_file  = '';
my $output_file = "tifa-$year-$month-$day.tgz";
my $help;

GetOptions(
    "in=s"  => \$input_file,
    "out=s" => \$output_file,
    "help"  => \$help
);

if (defined $help ) {
    print_help();
    exit();
}

check_options_validity();

my $subdir =  $output_file;
$subdir    =~ s/\.tgz$//;

#
# Perform all work in a randomly named subdirectory
#
my $tmpdir =  sprintf("tmp_%x", int rand(999999999));

system("mkdir -p ./$tmpdir/$subdir");
system("cp $input_file ./$tmpdir/$subdir/$input_file");
chdir("./$tmpdir/$subdir");
system("tar -xzvf $input_file");
system("rm $input_file");
chdir("..");

my @files = ();

find(sub {push(@files, $File::Find::name)}, "./$subdir");

foreach my $file (@files) {
    if (basename($file) =~ m/^\./) {
        system("rm $file\n");
        next;
    }
    if (basename($file) =~ m/\.(\w*)\.\w\w\w$/) {
        system("chmod u+wr $file\n");
        system("rm $file\n");
        next;
    }
}

system("tar -czvf $output_file $subdir");
system("mv $output_file ../");

chdir("..");
system("rm -fr ./$tmpdir");

my $hash;
my $line = '-' x (length($output_file) + 1);

my $hash_file = $output_file;
$hash_file =~ s/\.tgz/\.hash/;

open(OUT, ">$hash_file") or die("Cannot open $hash_file");

print(OUT "\n$output_file:\n");
print(OUT "$line\n");

$hash = `openssl md5 $output_file`;
$hash =~ s/\)=/\)\ =/;
print(OUT "        $hash");

$hash = `openssl sha1 $output_file`;
$hash =~ s/\)=/\)\ =/;
print(OUT "       $hash");

$hash = `openssl dgst -ripemd160 $output_file`;
$hash =~ s/\)=/\)\ =/;
print(OUT "  $hash\n");

close(OUT);

#-------------------------------------------------------------------------------
#   Subroutines
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
sub print_help {
    my $name = basename($0);

    my $title = "$name - Clean problematic .tgz archives created by SCons";
    my $line  = '-' x length($title);
    print << "EOF";

$title
$line
EOF

print << "EOF";

This script opens a .tgz file archive, cleans its content (garbage files and
files beginning by a dot) and repackages it in a new archive.

This is useful with the archives created by SCons since they do not expand in
a subdirectory but put everything in the current directory (not quite smart!).
The .tgz files produced by this script does not have this problem.

Also, outputs the MD5, SHA1 and RIPEMD hashes of the produced archive in an
output file ending with the extension .hash and using as prefix the archive
file name (minus its .tgz extension).

Usage: $name --in <input_file> --out <output_file> [--help]

  --in=s
      (mandatory)
      Name of the input .tgz archive file to clean.

  --out=s
      (optional, default value used if none provided)
      Name of the output cleaned .tgz archive file.
      Default: tifa-year-month-day.tgz where year, month and day are the
               current year (on 4 digits), month and day (on 2 digits).

  --help
      Prints this help.

EOF
  exit();
}
#-------------------------------------------------------------------------------
sub check_options_validity {
    if (! defined $input_file) {
        print("ERROR: No input file provided\n");
        exit();
    }
    if (! -e $input_file) {
        print("ERROR: Cannot find input file $input_file\n");
        exit();
    }
    if (($input_file !~ m/\.tgz$/) && ($input_file !~ m/\.tar\.gz$/)){
        print("ERROR: Input file $input_file is not a .tgz archive\n");
        exit();
    }
    if (-e $output_file) {
        print("ERROR: Output file $output_file already exists\n");
        exit();
    }
}
#-------------------------------------------------------------------------------
