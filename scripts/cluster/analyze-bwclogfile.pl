#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw/min max/;
use DateTime::Format::Strptime;
use Data::Dumper;

my $dateparser = new DateTime::Format::Strptime(
    pattern => '%a %b %d %T %Z %Y',
    time_zone => 'Europe/Paris',
    on_error        => 'croak',
) or die;

sub get_nodelist {
    my ($fh, $res) = @_;
    my $mpi = $res->{'mpi'} or die "mpi=$res->{'mpi'} ?";
    my $thr = $res->{'thr'} or die "thr=$res->{'thr'} ?";
    $mpi =~ /^(\d+)x(\d+)$/ or die; $mpi=[$1,$2];
    $thr =~ /^(\d+)x(\d+)$/ or die; $thr=[$1,$2];
    my @rows=();
    for(my $i = 0 ; $i < $mpi->[0] ; $i++) {
        my $row = <$fh>;
        my @r=();
        for(my $j = 0 ; $j < $mpi->[1] ; $j++) {
            $row =~ s/^\s*(\d+)// or die "unexpected: $row";
            push @r, $1;
            for(my $jt = 1 ; $jt < $thr->[1] ; $jt++) {
                $row =~ s/^\s*\.// or die "unexpected: $row";
            }
        }
        push @rows, \@r;
        for(my $it = 1; $it < $thr->[0] ; $it++) {
            my $row = <$fh>;
            for(my $x = 0 ; $x < $mpi->[1] * $thr->[1] ; $x++) {
                $row =~ s/^\s*\.// or die "unexpected: $row";
            }
        }
    }
    my @hrz = ();
    push @hrz, @$_ for @rows;
    my @vrt = ();
    for(my $j = 0 ; $j < $mpi->[1] ; $j++) {
        for(my $i = 0 ; $i < $mpi->[0] ; $i++) {
            push @vrt, $rows[$i]->[$j];
        }
    }

    return \@vrt;
}

sub analyze_file {
    my $f = shift ;
    return undef if -z $f;
    my $res = {};
    my @fs=($f);
    if ($f =~ m{OAR\.(?:[\w_-]+\.)?(\d+)\.stdout}) {
        $res->{'jobid'} = $1;
        my $g = $f;
        $g =~ s/stdout/stderr/;
        push @fs, $g;
    } else {
        $res->{'jobid'} = $f;
    }
    $res->{'tail'}="";
    $res->{'rev'}="???";
    $res->{'steps'}=[];
    my @tcpus=();
    my @tcomms=();
    my $time;
    for my $x (@fs) {
        open F, $x or die "$x: $!";
        while (<F>) {
            if (/Command\sline.*--mksol/) {
                $res->{'mksol'}=1;
            }
            if (/Abort.*error/) {
                $res->{'error'}=$_;
            }
            if (/mpi(?:exec|run).*(?:krylov|mksol)/) {
                $res->{'t1'}=$time;
            }
            if (/^Started: (.*)$/) {
                $time = $dateparser->parse_datetime($1)->epoch();
                if (!defined($res->{'t0'})) {
                    $res->{'t0'} = $time;
                }
            }
            if (/^(...\s...\s\d+\s\d+:\d+:\d+\s\w+\s\d+)$/) {
                $time = $dateparser->parse_datetime($1)->epoch();
                if (!defined($res->{'t0'})) {
                    $res->{'t0'} = $time;
                }
            }
            if (/^iteration \d+$/) {
                @tcpus=();
                @tcomms=();
            }
            if (/^\|/) {
                while (s/(\d+\.\d+)@\d+%\+(\d+\.\d+)@\d+%?//) {
                    push @tcpus, $1;
                    push @tcomms, $2;
                }
            }
            if (/revision (?:svn \d+ -- )?git (\w+)$/) {
                $res->{'rev'}=$1;
            }
            if (/^# \((?:svn \d+ -- )?([0-9a-f]+)(\+)?\)/) {
                $res->{'rev'}=$1;
                $res->{'rev'}.='#M' if $2;
            }
            if (/^(?:Total|.*?:) (\d+) rows (\d+) cols .*?\s*(\d+) coeffs$/) {
                my $N = int($1/1.0e6);
                my $w = int($3/$1 + 0.5);
                $res->{'matdesc'}="${N}M_$w";
            }
            if (m{\bmpi=(\d+x\d+)}) { $res->{'mpi'}=$1; }
            if (m{\bthr=(\d+x\d+)}) { $res->{'thr'}=$1; }
            if (m{^(\d+) nodes on (\w+)}) {
                $res->{'cluster'}=sprintf('%-10s', $2);
                $res->{'nodes'}=$1;
                $res->{'nodelist'} = get_nodelist(*F{IO}, $res);
            }
            if (m{^\s*Running on (\w+), (\d+) nodes}) {
                $res->{'cluster'}=sprintf('%-10s', $1);
                $res->{'nodes'}=$2;
            }
            if (m{^(\d+) nodes, no common name prefix}) {
                $res->{'cluster'}=sprintf('%-10s', "<no prefix>");
                $res->{'nodes'}=$1;
                $res->{'nodelist'} = get_nodelist(*F{IO}, $res);
            }
            if (m{CPU:.*(\d+\.\d+)\ss/iter}) { $res->{'tcpu'}=$1; }
            if (m{COMM:.*(\d+\.\d+)\ss/iter}) { $res->{'tcomm'}=$1; }
            if (m{N=(\d+).*ETA\s.*\[(\d+\.\d+)\ss/iter\]}) {
                die "$f: weird" unless defined $res->{'t1'};
                $res->{'niter'}=$1;
                $res->{'ttot'}=$2;
                push @{$res->{'steps'}}, $1;
            }
            if (m{CPU:.*(\d+\.\d+)\sms/iter}) { $res->{'tcpu'}=$1/1000.0; }
            if (m{COMM:.*(\d+\.\d+)\sms/iter}) { $res->{'tcomm'}=$1/1000.0; }
            if (m{N=(\d+).*ETA\s.*\[(\d+\.\d+)\sms/iter\]}) {
                $res->{'niter'}=$1;
                $res->{'ttot'}=$2/1000.0;
            }
            if (m{^iteration (\d+)$}) {
                $res->{'niter'}=$1;
            }
            if (m{start=(\d+).*end=\1}) {
                $res->{'stupid'}=1;
            }
        }
        close F;
    }
    if (defined(my $e = $res->{'error'})) {
        if ($e =~ /dest rank=(\d+)/) {
            $res->{'culprit'}=$res->{'nodelist'}->[$1];
        }
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
    $res->{'tx'} = $time;
    die "$f: no t0" unless defined $res->{'t0'};
    if (defined($res->{'t1'})) {
        $res->{'import'} = $res->{'t1'} - $res->{'t0'};
        $res->{'useful'} = $res->{'tx'} - $res->{'t1'};
    } else {
        $res->{'import'} = $res->{'tx'} - $res->{'t0'};
        $res->{'useful'} = 0;
    }

    for my $k (qw/cluster nodes mpi thr/) {
        if (!defined($res->{$k})) {
            print STDERR "$f : missing info for \"$k\"\n";
            return;
        }
    }

    if ($res->{'t1'}) {
        # We should have some data
        my @missing=();
        for my $k (qw/jobid tcpu tcomm ttot niter/) {
            next if defined $res->{$k};
            push @missing, $k;
            $res->{$k}='???';
        }
        if (@missing) {
            warn "$f : missing " . join(" ", @missing) . "\n" unless $res->{'stupid'};
        }
    }

    return $res;
}


my $categories = {
    error => [],
    worked => [],
    invain => [],
    stupid => [],
};

my $matrices={};
my @a=();
for my $f (@ARGV) {
    my $res = analyze_file($f);
    next unless $res;
    if ($res->{'error'}) {
        push @{$categories->{'error'}}, $res;
        next;
    # } elsif ($res->{'stupid'}) {
    # push @{$categories->{'stupid'}}, $res;
    # next;
    } elsif (!defined($res->{'t1'})) {
        push @{$categories->{'invain'}}, $res;
        next;
    }
    push @{$categories->{'worked'}}, $res;
    my $score;
    if ($res->{'ttot'} ne '???') {
        $score=$res->{'nodes'}*$res->{'ttot'};
    } else {
        $score=9999;
    }
    my $slot = $res->{'matdesc'};
    if ($res->{'mksol'}) {
        $slot .= " (mksol)";
    }
    $matrices->{$slot}++;   # autovivify
    push @a, [$score, $res];
}
@a = sort { $a->[0] <=> $b->[0] } @a;


# Count, for each cluster: the number of steps, the total I/O time, and
# the total cpu time. Exclude buggy jobs, and stupid jobs as well.

{
my $global=[0,0,0];
my $tally={};
for my $res (@{$categories->{'worked'}}) {
    my $cluster=$res->{'cluster'};
    my $n = scalar @{$res->{'steps'}};
    my $i = $res->{'import'};
    my $u = $res->{'useful'};
    $tally->{$cluster}->[0]+=$n; 
    $tally->{$cluster}->[1]+=$i; 
    $tally->{$cluster}->[2]+=$u; 
    $global->[0] += $n;
    $global->[1] += $i;
    $global->[2] += $u;
}

for my $res (@{$categories->{'invain'}}) {
    my $cluster=$res->{'cluster'};
    my $i = $res->{'import'};
    $tally->{$cluster}->[1]+=$i; 
    $global->[1] += $i;
}

print "Per-cluster:\n";
for my $c (keys %$tally) {
    my $x = $tally->{$c};
    printf "$c: %d iters, %d I/O, %d useful (%.1f%%)\n",
        $x->[0], $x->[1], $x->[2], 100.0 * $x->[2]/($x->[1]+$x->[2]);
}
my $x = $global;
printf "global: %d iters, %d I/O, %d useful (%.1f%%)\n",
    $x->[0], $x->[1], $x->[2], 100.0 * $x->[2]/($x->[1]+$x->[2]);
}

for my $m (sort { $matrices->{$b} <=> $matrices->{$a} } (keys %$matrices)) {

    my $t = "Matrix: $m ($matrices->{$m})";
    my $l = int((80-length($t)-4)/2);
    print '-' x $l, " $t ", '-' x $l, "\n";

    for (@a) {
        my $res = $_->[1];
        next if $res->{'matdesc'} ne $m;
        my $text="$res->{'nodes'}" .
        " $res->{'cluster'}\t($res->{'mpi'} $res->{'thr'})" .
        "\t$res->{'tcpu'}+$res->{'tcomm'}=$res->{'ttot'}" .
        "  [$res->{'rev'} \@$res->{'niter'} $res->{'jobid'}]$res->{'tail'}";
        print "$text\n";
    }
}
