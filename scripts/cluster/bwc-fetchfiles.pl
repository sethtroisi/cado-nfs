#!/usr/bin/perl -w

use warnings;
use strict;
use File::Spec;
use Time::HiRes qw/time/;
use Data::Dumper;
use List::Util 'shuffle';

my $ssh = $ENV{'SSH'};
if (!defined($ssh)) {
    if (-x "/usr/bin/oarsh") {
        $ssh='oarsh';
    } else {
        $ssh='ssh';
    }
}

select(STDOUT); $|=1;

my @storehosts = qw/
                   rsa768.rennes.grid5000.fr
                   /;

# nfs.nancy.grid5000.fr


my $server_module = 'cado';

my $timeout_min = 2;
my $timeout_max = 60;
my $timeout_max_on_friend_failure = 0; # typically less than $timeout_max
my $quiet=0;
my $put=0;

my $dispatch = {};

my $md5db={};

my $ldir="/tmp";
my $port=8873;

my $start_time = time;

my $bytes_per_iteration_a = 512;   # m*n/8
my $bytes_vector = 441833584;      # beware, this varies per splitting
my $bytes_eval = 441833584;        # n/64 times more.
my $mksol_coarsegrain = 100000;    #

my $verbose = 0;

sub local_progress_files {
    my $res=[];
    my %vput = ();
    my %sput = ();
    my $vre = qr/^V(\d+-\d+)\.(\d+)\.([\da-f]+)$/;
    my $are = qr/^A(\d+-\d+)\.(\d+-\d+)$/;
    my $ere = qr/^S(\d+-\d+)\.(\d+)\.([\da-f]+)$/;

    opendir DIR, $ldir;
    my @all = readdir DIR;
    my @vecs = grep(/$vre/, @all);
    my @scps = grep(/$are/, @all);
    my @evals = grep(/$ere/, @all);
    closedir DIR;

    my @x=();
    for my $f (@scps) { $f =~ /$are/ or next; push @x, [ $1, $f ]; }
    for my $f (@evals) { $f =~ /$ere/ or next; push @x, [ $1, $f ]; }
    for my $f (@vecs) { $f =~ /$vre/ or next; push @x, [ $1, $f ]; }

    @x = sort { $a->[0] <=> $b->[0] } @x;

    for my $xf (@x) {
        my $f = $xf->[1];
        if ($f =~ /$are/) {
            my $j = $1;
            my $n0 = $2;
            my $n1 = $3;
            my @s = stat("$ldir/$f");
            if (!-f _ ) { do print STDERR "\n$f: vanished\n"; next; };
            my $esize = $bytes_per_iteration_a * ($n1-$n0);
            if ($s[7] != $esize) {
                print STDERR "\n$f: size $s[7], waits $esize before transferring\n"
                unless $quiet;
            } else {
                push @$res, "$j/$f";
            }
        } elsif ($f =~ /$ere/) {
            my $j = $1;
            my $n0 = $2;
            my $n1 = $3;
            my $sv = int($2/$mksol_coarsegrain);
            my @s = stat("$ldir/$f");
            if (!-f _ ) { do print STDERR "\n$f: vanished\n"; next; };
            my $esize = $bytes_eval;
            if ($s[7] < $esize) {
                print STDERR "\n$f: size $s[7], waits $esize before transferring\n"
                unless $quiet;
            } else {
                push @$res, "$j,$sv/$f";
            }
        } elsif ($f =~ /$vre/) {
            my $j = $1;
            my $k = $2;
            my @s = stat("$ldir/$f");
            if (!-f _ ) { do print STDERR "\n$f: vanished\n"; next; };
            my $esize = $bytes_vector;
            if ($s[7] < $esize) {
                print STDERR "\n$f: size $s[7], waits $esize before transferring\n";
            } else {
                push @$res, "$j/$f";
            }
        } else {
            print STDERR "$f --> what is this ???\n";
            next;
        }
    }
    return $res;
}

my @orig_argv=@ARGV;

while (defined($_=shift(@ARGV))) {
    if (/^--dispatch-list$/) {
        my $cachelist = shift @ARGV or die "$_ needs 1 arg";
        my @hlist;
        open F, $cachelist or die "$cachelist: $!";
        while (defined(my $x = <F>)) {
            my ($host, @parts) = split ' ', $x;
            push @hlist, $host;
            $dispatch->{$host}=[] unless exists $dispatch->{$host};
            push @{$dispatch->{$host}}, @parts;
        }
        close F;
        next;
    }
    if (/^--localdir$/) {
        $ldir = shift @ARGV or die "$_ needs 1 arg";
        next;
    }
    if (/^--fetch$/) {
        my $h = shift @ARGV or die "$_ needs 2 args";
        my $f = shift @ARGV or die "$_ needs 2 args";
        $h =~ s/\..*$//;
        $dispatch->{$h}=[$f];
        next;
    }
    if (/^--put$/) { $put=1; next; }
    if (/^-q$/) { $quiet=1; next; }
    if (/^-v$/) { $verbose=1; next; }
    if (/^--timeout$/) {
        my $x = shift @ARGV or die "$_ needs 1 arg";
        $x -= $timeout_min;
        if ($x > 0) {
            $timeout_max=$x;
        }
        next;
    }
    if (/^--md5db$/) {
        my $f = shift @ARGV or die "$_ needs 1 arg";
        open F, $f or die "$f: $!";
        while (defined($_=<F>)) {
            chomp($_);
            my ($sum, $path) = split(' ', $_);
            $md5db->{$path}=$sum;
        }
        close F;
        next;
    }
    if (/^--local-progress-files$/) {
        my $res = &local_progress_files();
        $port=8874;
        $dispatch->{'localhost'}=$res;
        ## print Dumper($dispatch);
        next;
    }
    die "Usage: $0 [--put] [--dispatch-list <file> | --fetch <host> <path> | --local-progress-files]";
}

unless ($quiet) {
    print "## Calling $0 ", join(" ", @orig_argv), "\n";
}

my %ch=();
for my $h (keys %$dispatch) {
    my $child = fork;
    if (!defined($child)) { die "fork(): $!"; }
    if ($child == 0) {
        # child.
        select(STDOUT); $|=1;
        my $nfiles = scalar @{$dispatch->{$h}};
        MYFILES: for my $i (shuffle (0..$nfiles-1)) {
            my $z = $dispatch->{$h}->[$i];
            my $pos = sprintf "[%d/%d]", (1+$i), $nfiles;
            my ($base,$subdir);
            $z =~ /^((?:.*\/)?)([^\/]+)$/ or die "$z: bad pattern";
            $base = $2;
            $subdir = $1;
            my $server;
            # print "$cmd\n";
            my $t0 = time;
            my ($u0, $u1);
            my $tw=0;
            my $header = "[$h] $base $pos";
            my $dt;
            if (!$put && scalar keys %$md5db) {
                if (system("$ssh -q -n $h test -f $ldir/$base") != 0) {
                    # print "$header not there yet\n";
                } elsif (!exists($md5db->{$base})) {
                    $dt = sprintf('%.1f',  time-$start_time);
                    print "$dt $header md5 unknown\n";
                } else {
                    my $text = `$ssh -q -n $h md5sum $ldir/$base`;
                    chomp($text);
                    my $digest;
                    if ($text =~ m{^([0-9a-f]{32})\s+$ldir/$base$}) {
                        $digest=$1;
                    }
                    if (!defined($digest)) {
                        $dt = sprintf('%.1f',  time-$start_time);
                        print "$dt $header cannot check md5 [[$text]]\n";
                    } elsif ($digest eq $md5db->{$base}) {
                        $dt = sprintf('%.1f',  time-$start_time);
                        print "$dt $header md5 ok ($digest)\n";
                        next;
                    } else {
                        $dt = sprintf('%.1f',  time-$start_time);
                        my $argh="$dt $header BAD MD5 ! ($digest != $md5db->{$base})\n";
                        print $argh;
                        print STDERR $argh;
                    }
                }
            }
            my $rsync_output;
            while (1) {
                my $cachefile = "$ENV{'HOME'}/.cache.$base";
                my $cache_server;
                $server = undef;
                my $spath="$server_module/$subdir";
                if ($put) {
                    my $cmd="$ssh -q -n $h test -f $ldir/$base";
                    system $cmd;
                    if ($? != 0) {
                        print STDERR "$base: vanished [$cmd]\n";
                        exit 1;
                    }
                    system "$ssh -q -n $h test -f $ldir/done.$base";
                    if ($? == 0) {
                        $rsync_output='';
                        last;
                    }
                } else {
                    # We rely on the availability of an nfs-shared $HOME.
                    CACHE: {
                        my $chost = readlink $cachefile or last CACHE;
                        $server=$chost; 
                        $dt = sprintf('%.1f',  time-$start_time);
                        # printf("$dt $header Using $chost as alternative server\n");
                        $cache_server = $chost;
                        $spath="tmp/";
                    }
                }
                if (!defined($server)) {
                    $server = $storehosts[int(rand(scalar @storehosts))];
                }
                my $cmd = "";
                $cmd = "$ssh -q -n $h " unless $h eq 'localhost';
                $cmd .= "rsync -av --port=$port ";
                # $cmd .= join(" ", map {
                # "nfs.lyon.grid5000.fr::rsa768/${nh}x${nv}/rsa768.comp.local.${nh}x${nv}.$_"
                # } @{$dispatch->{$h}});
                if ($put) {
                    $cmd .= "$ldir/$base";
                    $cmd .= " ${server}::$spath";
                } else {
                    $cmd .= "${server}::$spath$base";
                    $cmd .= " $ldir/";
                }
                $cmd .= " 2>/dev/null";
                # print "$header <-- $server:$port\n";
                $u0 = time;
                print "$cmd\n" if $verbose;
                $rsync_output=`$cmd`;
                my $res = $?;
                $u1 = time;
                if ($res == 0) {
                    if ($put) {
                        system "$ssh -q -n $h touch $ldir/done.$base";
                    } else {
                        symlink $h, $cachefile;
                    }
                    last;
                }
                my $s = $timeout_min + int(rand($timeout_max-$timeout_min));
                $dt = sprintf('%.1f',  time-$start_time);
                if (defined($cache_server)) {
                    print "$dt $header removing $cachefile (error on $server)\n";
                    unlink $cachefile;
                    # sleep less in this case, as it's very probably an
                    # easy miss.
                    if ($timeout_max_on_friend_failure) {
                        $s = $timeout_min + int(rand($timeout_max_on_friend_failure-$timeout_min));
                    } else {
                        $s = 0;
                    }
                }
                if ($s) {
                    print "$dt $header retrying in $s seconds [$server]\n";
                    $tw+=$s;
                    sleep($s);
                }
            }
            my $t1 = time;
            my $kb = `$ssh -q -n $h du -k $ldir/$base`;
            chomp($kb);
            $kb =~ /^(\d+)/;
            $kb=$1;
            my $mb = $kb >> 10;
            if ($rsync_output =~ /$base/) {
                if ($port == 8874) {
                    printf "\n";
                }
                $dt = sprintf('%.1f',  time-$start_time);
                printf("$dt $header ($mb MB): %.1f seconds [%.1f MB/s], %.1f total, %.1f waittime [$server]\n",
                    ($u1-$u0), $mb/($u1-$u0), ($t1-$t0), $tw);
            }
            # Why ? Because it yields before other processes, thereby
            # improving fairness of the scheduling. It's mild, but
            # otherwise we would be too aggressive in grabbing the server
            # attention for as long as _we_ like.
            sleep(1);
        }
        exit 0;
   }
   $ch{$child}=$h;
}

my $remaining = scalar keys %ch;
my $finished = 0;

while ((my $pid=wait) > 0) {
    my $who = $ch{$pid};
    $finished++;
    $remaining--;
    delete $ch{$pid};
    my @rem = values %ch;
    @rem = map { s/^[\D-]+(\d+).*$/$1/; $_ } @rem;
    if (scalar keys %$dispatch > 1) {
        my $dt = sprintf('%.1f',  time-$start_time);
        print "$dt [$who] ready ($finished ready, $remaining remaining: ", join(' ', @rem), ")\n";
    }
}
