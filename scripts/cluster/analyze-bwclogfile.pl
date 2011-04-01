#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw/min max/;

sub analyze_file {
    my $f = shift ;
    my $res = {};
    my @fs=($f);
    if ($f =~ m{OAR\.(?:[\w-]+\.)?(\d+)\.stdout}) {
        $res->{'jobid'} = $1;
        my $g = $f;
        $g =~ s/stdout/stderr/;
        push @fs, $g;
    } else {
        $res->{'jobid'} = $f;
    }
    $res->{'tail'}="";
    $res->{'rev'}="???";
    my @tcpus=();
    my @tcomms=();
    for my $x (@fs) {
        open F, $x or die "$x: $!";
        while (<F>) {
            if (/^iteration \d+$/) {
                @tcpus=();
                @tcomms=();
            }
            if (/^\|/) {
                while (s/(\d+\.\d+)@\d+%\+(\d+\.\d+)@\d+%//) {
                    push @tcpus, $1;
                    push @tcomms, $2;
                }
            }
            if (/revision (?:svn \d+ -- )?git (\w+)$/) {
                $res->{'rev'}=$1;
            }
            if (/^# \((?:svn \d+ -- )?git (\w+)(\s\+mods)?\)/) {
                $res->{'rev'}=$1;
                $res->{'rev'}.='#M' if $2;
            }
            if (/^Total (\d+) rows (\d+) cols .* (\d+) coeffs$/) {
                my $N = int($1/1.0e6);
                my $w = int($3/$1 + 0.5);
                if ($1 != 55229198 || $2 != 55229006 || $3 != 4525926812) {
                    $res->{'tail'}=" ! ${N}M_$w";
                }
            }
            if (m{^(\d+) nodes on (\w+)}) {
                $res->{'cluster'}=sprintf('%-10s', $2);
                $res->{'nodes'}=$1;
            }
            if (m{^(\d+) nodes, no common name prefix}) {
                $res->{'cluster'}=sprintf('%-10s', "<no prefix>");
                $res->{'nodes'}=$1;
            }
            if (m{\bmpi=(\d+x\d+)}) { $res->{'mpi'}=$1; }
            if (m{\bthr=(\d+x\d+)}) { $res->{'thr'}=$1; }
            if (m{CPU:.*(\d+\.\d+)\ss/iter}) { $res->{'tcpu'}=$1; }
            if (m{COMM:.*(\d+\.\d+)\ss/iter}) { $res->{'tcomm'}=$1; }
            if (m{N=(\d+).*ETA\s.*\[(\d+\.\d+)\ss/iter\]}) {
                $res->{'niter'}=$1;
                $res->{'ttot'}=$2;
            }
            if (m{CPU:.*(\d+\.\d+)\sms/iter}) { $res->{'tcpu'}=$1/1000.0; }
            if (m{COMM:.*(\d+\.\d+)\sms/iter}) { $res->{'tcomm'}=$1/1000.0; }
            if (m{N=(\d+).*ETA\s.*\[(\d+\.\d+)\sms/iter\]}) {
                $res->{'niter'}=$1;
                $res->{'ttot'}=$2/1000.0;
            }
        }
        close F;
    }
    # if ((!defined($res->{'tcpu'}) || !defined($res->{'tcomm'}) ||
    if (scalar @tcpus && scalar @tcomms)
    {
        my $tcpu=0; $tcpu+=$_ for (@tcpus); $tcpu /= scalar @tcpus;
        my $tcomm=0; $tcomm+=$_ for (@tcomms); $tcomm /= scalar @tcomms;
        my $ttot=$tcpu+$tcomm;
        my $min_cpu = int(100.0 * (1-(min(@tcpus) / $tcpu)));
        my $max_cpu = int(100.0 * ((max(@tcpus) / $tcpu)-1));
        my $dev_cpu = max($min_cpu, $max_cpu);
        my $min_comm = int(100.0 * (1-(min(@tcomms) / $tcomm)));
        my $max_comm = int(100.0 * ((max(@tcomms) / $tcomm)-1));
        my $dev_comm = max($min_comm, $max_comm);
        $res->{'tcpu'} = sprintf('%.2f±%d', $tcpu, $dev_cpu);
        $res->{'tcomm'} = sprintf('%.2f±%d', $tcomm, $dev_comm);
        $res->{'ttot'} = sprintf('%.2f', $ttot);
    }

    for my $k qw/cluster nodes mpi thr/ {
        if (!defined($res->{$k})) {
            print STDERR "$f : missing info for \"$k\"\n";
            exit;
        }
    }
    for my $k qw/jobid tcpu tcomm ttot niter/ {
        if (!defined($res->{$k})) {
            warn "$f : missing info for \"$k\"\n";
            $res->{$k}='???';
        }
    }
    return $res;
}

my @a=();
for my $f (@ARGV) {
    my $res = analyze_file($f);
    next unless $res;
    my $text="$res->{'nodes'}" .
            " $res->{'cluster'}\t($res->{'mpi'} $res->{'thr'})" .
            "\t$res->{'tcpu'}+$res->{'tcomm'}=$res->{'ttot'}" .
            "  [$res->{'rev'} \@$res->{'niter'} $res->{'jobid'}]$res->{'tail'}\n";
    my $score;
    if ($res->{'ttot'} ne '???') {
        $score=$res->{'nodes'}*$res->{'ttot'};
    } else {
        $score=9999;
    }
    push @a, [$score, $text];
}
@a = sort { $a->[0] <=> $b->[0] } @a;
for (@a) {
    print "$_->[1]";
}
