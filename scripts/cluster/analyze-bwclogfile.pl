#!/usr/bin/perl

sub analyze_file {
    my $f = shift ;
    my $res = {};
    if ($f =~ m{OAR\.\w+\.(\d+)\.stdout}) {
        $res->{'jobid'} = $1;
    } else {
        $res->{'jobid'} = $f;
    }
    $res->{'tail'}="";
    $res->{'rev'}="???";
    open F, $f or die "$f: $!";
    while (<F>) {
        if (/revision (?:svn \d+ -- )?git (\w+)$/) {
            $res->{'rev'}=$1;
        }
        if (/^Total (\d+) rows (\d+) cols .* (\d+) coeffs$/) {
            $N = int($1/1.0e6);
            $w = int($3/$1 + 0.5);
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
    for my $k qw/jobid cluster nodes mpi thr tcpu tcomm ttot niter/ {
        if (!defined($res->{$k})) {
            print STDERR "$f : missing info for \"$k\"\n";
            return;
        }
    }
    return $res;
}

my @a=();
for my $f (@ARGV) {
    my $res = analyze_file($f);
    next unless $res;
    my $text="$res->{'cluster'}" .
            " $res->{'nodes'} ($res->{'mpi'} $res->{'thr'})" .
            "\t$res->{'tcpu'}+$res->{'tcomm'}=$res->{'ttot'}" .
            "\t[$res->{'rev'} \@$res->{'niter'} $res->{'jobid'}]$res->{'tail'}\n";
    my $score=$res->{'nodes'}*$res->{'ttot'};
    push @a, [$score, $text];
}
@a = sort { $a->[0] <=> $b->[0] } @a;
for (@a) {
    print "$_->[1]";
}
