#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

sub format_dhms {
    my $sec = shift;
    my ($d, $h, $m);
    $d = int ( $sec / 86400 ); $sec = $sec % 86400;
    $h = int ($sec / 3600 ); $sec = $sec % 3600;
    $m = int ($sec / 60 ); $sec = $sec % 60;
    return "$d"."d:$h"."h:$m"."m:$sec"."s";
}

my $wdir = abs_path(dirname(dirname($0)));
opendir DIR, "$wdir/merge_rels/"
    or die "Cannot open directory `$wdir/merge_rels/': $!\n";
my @files = grep /\.rels\.\d+-\d+\.gz$/, readdir DIR;
closedir DIR;

my $compl = 0;
my $trunc = 0;
my $cpu_time = 0;
my $count;
my $line;
my $rfile = 1;
my $nfiles = scalar (@files);
for (@files) {
    print "Reading $_... ($rfile/$nfiles)\r";
    $rfile++;
    $count = 0;
    open FILE, "gzip -dc $wdir/merge_rels/$_ |" or die "Cannot open `$_' for reading: $!.\n";
    while ($line = <FILE>) {
        if ($line =~ /cpu time (\S+)s \[/) {
            $cpu_time += $1;
            $compl += $count;
            $count = 0;
        } elsif ($line =~ /truncated/) {
            $trunc += $count;
            $count = 0;
        }
        $count++ unless ($line =~ /^#/);
    }
    close FILE;
}

my $total = $compl + $trunc;
my $est = int ( $cpu_time * $total / $compl );
print "completed: $compl, truncated: $trunc, total: $total.\n".
      "completed cpu time       : $cpu_time s\n".
      "                           ".format_dhms($cpu_time).".\n".
      "estimated total cpu time : $est s\n".
      "                           ".format_dhms($est).".\n";
