#!/usr/bin/env perl

use strict;
use warnings;
use File::Copy;


my $tmpdir="/tmp";
my $bindir="$HOME/cado-nfs/build/$HOSTNAME";
my $las="$bindir/sieve/las -checknorms";
my $makefb="$bindir/sieve/makefb";

my $polyfile="c155.poly";
my $qrange=1000;

my @list_I=(12, 13, 14);
my @list_qmin=(30000000, 100000000, 400000000);
my @list_arlim=([16e6,32e6], [32e6,64e6], [64e6,150e6]);
my @list_lpb=(28, 29);

for my $arlim (@list_arlim) {
    copy($polyfile, "$tmpdir/poly");
    open FILE, ">> $tmpdir/poly";
    print FILE "rlim: " . ${$arlim}[0] . "\n";
    print FILE "alim: " . ${$arlim}[1] . "\n";
    print FILE "lpbr: " . $list_lpb[0]. "\n";
    print FILE "lpba: " . $list_lpb[0]. "\n";
    print FILE "mfbr: " . 2*$list_lpb[0] . "\n";
    print FILE "mfba: " . 2*$list_lpb[0] . "\n";
    print FILE "rlambda: 2.3\n";
    print FILE "alambda: 2.3\n";
    close FILE;

    system("$makefb -poly $tmpdir/poly > $tmpdir/roots");

    for my $I (@list_I) {
        for my $qmin (@list_qmin) {
            for my $lpb (@list_lpb) {
                copy($polyfile, "$tmpdir/poly");
                open FILE, ">> $tmpdir/poly";
                print FILE "rlim: " . ${$arlim}[0] . "\n";
                print FILE "alim: " . ${$arlim}[1] . "\n";
                print FILE "lpbr: " . $lpb . "\n";
                print FILE "lpba: " . $lpb . "\n";
                print FILE "mfbr: " . 2*$lpb . "\n";
                print FILE "mfba: " . 2*$lpb . "\n";
                print FILE "rlambda: 2.3\n";
                print FILE "alambda: 2.3\n";
                close FILE;

                my $qmax = $qmin+$qrange;
                my $cmd = "$las -I $I -poly $tmpdir/poly" .
                 " -fb $tmpdir/roots -q0 $qmin -q1 $qmax" .
                 " >& $tmpdir/rels";
                system("$cmd");
                open FILE, "<$tmpdir/rels";
                my $lastline;
                my $rep = -1;
                my $rate = -1;
                my $nb_q = -1;
                while (<FILE>) {
                    if (/# Total (\d+) reports \[(\d+.\d+)s\/r\]/) {
                        $rep=$1;
                        $rate=$2;
                        next;
                    } 
                    if (/# Average J=\d+ for (\d+) special-q/) {
                        $nb_q=$1;
                        next;
                    }
                }
                close FILE;
                if ($rep==-1 || $rate==-1 || $nb_q==-1) {
                    die "Could not get info from $tmpdir/rels";
                }
                my $repq = $rep/$nb_q;
                $repq = sprintf("%.1f", $repq);
                print "I=$I rlim=". ${$arlim}[0] . 
                  " alim=" . ${$arlim}[1] .
                  " lpb=$lpb qmin=$qmin:\t rel/q=$repq rate=$rate s/r\n";
            }
        }
    }
}


