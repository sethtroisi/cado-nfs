#!/usr/bin/env perl

use warnings;
use strict;
use File::Spec;
use Time::HiRes qw/time/;
use Data::Dumper;
use List::Util 'shuffle';
use IPC::Open3;


my $ssh = $ENV{'SSH'};
if (!defined($ssh)) {
    if (-x "/usr/bin/oarsh") {
        $ssh='oarsh';
    } else {
        $ssh='ssh';
    }
}

select(STDOUT); $|=1;

my %storehosts = ();
my %failed_storehosts = ();


# nfs.nancy.grid5000.fr


my $server_module = 'cado';

my $timeout_min = 2;
my $timeout_max = 60;
my $quiet=0;
my $put=0;

my $dispatch = {};

my $md5db={};

my $ldir="/tmp";

my $start_time = time;

my $bytes_per_iteration_a = 512;   # m*n/8
my $bytes_vector = 441833584;      # beware, this varies per splitting
my $bytes_eval = 441833584;        # n/64 times more.
my $mksol_coarsegrain = 100000;    #

my $verbose = 0;

sub node_server {
    # returns the URI of the locally running rsync server on compute node
    # $h
    my $h = shift;
    return "rsync://$h:8873/tmp";
}

sub backtick {
    my $maxtime;
    if ($_[0] eq '--timeout') {
        shift;
        $maxtime = shift;
    }
    my $cmd = "@_";
    if ($verbose) {
        my $msg = "";
        $msg = "[timeout $maxtime] " if defined($maxtime);
        $msg .= $cmd;
        print STDERR "## $msg\n";
    }
    open NULL, "</dev/null" or die;
    my $pid = open3(*NULL,\*CHLD_OUT, \*CHLD_ERR, $cmd);
    close NULL;
    
    my $fds = {
        'out' => [ *CHLD_OUT{IO}, "", 1, "" ],
        'err' => [ *CHLD_ERR{IO}, "", 1, "" ],
    };
    
    my $t0 = time;
    while (scalar grep { $_->[2] } values %$fds) {
        my $rin = '';
        for my $v (values %$fds) {
            next if !$v->[2]; 
            vec($rin, fileno($v->[0]), 1) = 1;
        }
        
        my $rc;
        my $rout;
        while ($rc = select($rout=$rin, undef, undef, 1.0) == 0) {
            if (defined($maxtime) && time - $t0 > $maxtime) {
                kill 9, $pid;
                waitpid($pid,0) or warn "Can't wait() for children, weird";
                my $res = {
                    'rc' => -1,
                    'wait'=> -1,
                    'out' => $fds->{'out'}->[3],
                    'err' => $fds->{'err'}->[3],
                };
                return $res;
            }
        }
        if ($rc < 0) {
            print STDERR "Left select with $!\n";
            last;
        }
    
        for my $k (keys %$fds) {
            my $v = $fds->{$k};
            next unless vec($rout, fileno($v->[0]), 1);
            if (sysread $v->[0], my $x, 1024) { 
                $v->[1] .= $x;
                while ($v->[1]=~s/^(.*(?:\n|\r))//) {
                    # print $1 if $verbose;
                    # print "$k $1";
                    $v->[3] .= $1;
                }
            } else {
                if (length($v->[1])) {
                    $v->[3] .= $v->[1];
                    $v->[1] =~ s/^/$k /gm;
                    # print $v->[1] if $verbose;
                    $v->[1] = '';
                }
                $v->[2]=0;
            }
        }
    }
    waitpid($pid,0) or warn "Can't wait() for children, weird";
    my $res = {
        'rc' => $?>>8,
        'wait'=> $?,
        'out' => $fds->{'out'}->[3],
        'err' => $fds->{'err'}->[3],
    };
    return $res;
}

sub local_progress_files {
    my $res=[];
    my %vput = ();
    my %sput = ();
    my $vre = qr/^V(\d+-\d+)\.(\d+)$/;
    my $are = qr/^A(\d+-\d+)\.(\d+-\d+)$/;
    my $ere = qr/^S(\d+-\d+)\.(\d+)$/;

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

sub usage {
    die "Usage: $0 [--put] [--dispatch-list <file> | --on <host> <path> | --local-progress-files]";
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
    if (/^(--server|-s)$/) {
        my $x = shift @ARGV or die "$_ needs 1 arg";
        $x =~ m{(?:rsync://)?([^/:]+)(?:(\d+))?/(.*)$} or die "bad URI: $x";
        $storehosts{$x}=1;
        next;
    }
    if (/^--localdir$/) {
        $ldir = shift @ARGV or die "$_ needs 1 arg";
        my $x = $ldir;
        $x =~ s{/$}{};
        my @dirs=();
        while ($x =~ s{(/[^/]+)$}{}) {
            push @dirs, "$x$1";
        }
        for my $x (@dirs) {
            next if -d $x;
            print "Need to create directory $x !";
            mkdir $x;
        }
        next;
    }
    if (/^--on$/) {
        my $h = shift @ARGV or die "$_ needs 2 args";
        my $f = shift @ARGV or die "$_ needs 2 args";
        $h =~ s/\..*$//;
        $dispatch->{$h}=[] unless $dispatch->{$h};
        push @{$dispatch->{$h}}, split(' ', $f);
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
            next unless $path;
            next unless $sum;
            $path =~ s{^.*/}{};
            $md5db->{$path}=$sum;
        }
        close F;
        next;
    }
    if (/^--local-progress-files$/) {
        my $res = &local_progress_files();
        $dispatch->{'localhost'}=$res;
        ## print Dumper($dispatch);
        next;
    }
    usage;
}
usage unless scalar keys %storehosts;

sub fqdn {
    my $x=`host $_[0]`;
    $x=~/^(\S*)/;
    return $1;
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
        my @myfiles = shuffle @{$dispatch->{$h}};
        MYFILES: for my $i (0..$nfiles-1) {
            my $z = $myfiles[$i];
            my $pos = sprintf "[%d/%d]", (1+$i), $nfiles;
            my ($base,$subdir);
            $z =~ /^((?:.*\/)?)([^\/]+)$/ or die "$z: bad pattern";
            $base = $2;
            $subdir = $1;
            my $server;
            # print "$cmd\n";
            my $t0 = time;
            my ($u0, $u1, $filesize);
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
                    my $res = backtick("$ssh -q -n $h md5sum $ldir/$base");
                    die if $res->{'rc'} != 0;
                    my $text = $res->{'out'};
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
                my $cachedir = "$ENV{'HOME'}/.cache";
                mkdir $cachedir unless -d $cachedir;
                my $cachefile = "$cachedir/cache.$base";
                my $cache_server;
                $server = undef;
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
                        my $ok = $server =~
                        m{^(?:rsync://)([^/:]+)(?::(\d+))/(.*)$};
                        if (!$ok) {
                            print "$dt $header $cachefile contains bogus server location\n";
                            unlink $cachefile;
                            $server = undef;
                        }
                        if ($1 eq $h) {
                            print "$dt $header $cachefile points locally, ignoring\n";
                            $server = undef;
                        } else {
                            # print "$dt $header Using $chost as alternative server\n";
                            $cache_server = $chost;
                        }
                    }
                }
                if (!defined($server)) {
                    my @x = keys %storehosts;
                    if (scalar @x == 0) {
                        print STDERR "All servers failed\n";
                        for my $k (keys %failed_storehosts) {
                            my $v = $failed_storehosts{$k};
                            my $out = $v->[1] || '';
                            $out =~ s/^/\t/mg;
                            print STDERR "$k : status $v->[0]:\n$out";
                        }
                        exit(1);
                    }
                    $server = $x[int(rand(scalar @x))];
                }
                $server =~ m{(?:rsync://)?([^/:]+)(?:(\d+))?/(.*)$} or die "bad URI: $server";
                if ($put) {
                    $header = "[${h}->${server}] $base $pos";
                } else {
                    $header = "[${h}<-${server}] $base $pos";
                }
                my $cmd;
                my $cmd0 = "";
                $cmd0 = "$ssh -q -n $h " unless $h eq 'localhost';
                system($cmd0 . "mkdir -p $ldir");

                # Compute an a priori sleep value in case of failure
                my $s = $timeout_min + int(rand($timeout_max-$timeout_min));
                my $fatal;
                $u0 = time;

                my ($orig, $dest);
                if ($put) {
                    $orig = " $ldir/$base";
                    $dest = " $server/$subdir";
                } else {
                    if (defined($cache_server)) {
                        $orig = " $server/$base";
                    } else {
                        $orig = " $server/$subdir$base";
                    }
                    $dest = " $ldir/";
                }
                my $rc;
                # try to get the file size.
                $cmd = $cmd0 . "rsync -a $orig";
                my $getsize = backtick('--timeout', 2, $cmd);
                if ($getsize->{'rc'} == -1 || $getsize->{'rc'} == 5) {
                    $fatal = defined($cache_server);
                    $rc = "$getsize->{'rc'} (while reading size)";
                    goto FAILED_COMMAND;
                } elsif ($getsize->{'rc'} != 0) {
                    $fatal = 1;
                    $rc = "$getsize->{'rc'} (while reading size)";
                    goto FAILED_COMMAND;
                }
                my @gs = split(' ', $getsize->{'out'});
                $filesize = $gs[1];
                if (!defined($filesize) || $filesize !~ /^\d+/) {
                    $rc = "-1 (size format mismatch)";
                    $fatal = 1;
                    goto FAILED_COMMAND;
                }
                my $timeout_val = 1 + int($filesize / (4 * 1024 * 1024));

                $cmd = $cmd0 . "rsync -av $orig $dest";
                # We used to discard stderr.
                # $cmd .= " 2>&1";
                # print "$header <-- $server\n";
                print "$cmd\n" if $verbose;
                my $rsync_res=backtick('--timeout', $timeout_val, $cmd);
                $rsync_output = $rsync_res->{'out'};
                my $res = $rsync_res->{'rc'};
                if ($res == -1) {
                    $fatal = defined($cache_server);
                    $rc = "rsync timeout after $timeout_val s";
                    goto FAILED_COMMAND;
                }
                $res = $rsync_res->{'wait'};
                $u1 = time;
                if ($res == 0) {
                    if ($put) {
                        system "$ssh -q -n $h touch $ldir/done.$base";
                    } else {
                        unlink $cachefile;
                        symlink node_server(fqdn($h)), $cachefile or warn "$!";
                    }
                    last;
                }
                if ($res & 127) {
                    printf "$cmd: died with signal %d, %s coredump\n",
                    ($res & 127),  ($res & 128) ? 'with' : 'without';
                    die;
                }
                $rc = $res >> 8;
                if ($rc == 5) {
                    # This error is not fatal for the current server,
                    # because it is the error message we expect to get
                    # when the number of connections exceeds the
                    # configured limit.
                    $fatal = 0;
                } else {
                    $fatal = 1;
                }

FAILED_COMMAND:
                die unless defined $fatal;

                $dt = sprintf('%.1f',  time-$start_time);

                if  ($fatal) {
                    print "$cmd: exit code $rc\n";
                    # Assume it's fatal. This means we will try another
                    # server, which implies that we don't need to bother
                    # with sleeping.
                    $s = 0;
                    if (defined($cache_server)) {
                        print "$dt $header removing $cachefile\n";
                        unlink $cachefile;
                    } else {
                        delete $storehosts{$server};
                        $failed_storehosts{$server}=[$rc, $rsync_output];
                    }
                }
                if ($s) {
                    print "$cmd: non-fatal exit code $rc\n" if $verbose;
                    print "$dt $header retrying in ${s} s\n";
                    $tw+=$s;
                    sleep($s);
                }
            }
            my $t1 = time;
            my $mb;
            if (defined($filesize)) {
                $mb = $filesize / 1048576;
            } else {
                my $kb = backtick("$ssh -q -n $h du -k $ldir/$base");
                die unless $kb->{'rc'} == 0;
                $kb = $kb->{'out'};
                chomp($kb);
                $kb =~ /^(\d+)/ or die "backtick says: $kb";
                $kb=$1;
                $mb = $kb / 1024;
            }
            if ($rsync_output =~ /$base/) {
                # if ($port == 8874) { printf "\n"; }
                $dt = sprintf('%.1f',  time-$start_time);
                printf("$dt $header (%.1f MB): %.1f s [%.1f MB/s], %.1f total, %.1f w\n",
                    $mb, ($u1-$u0), $mb/($u1-$u0), ($t1-$t0), $tw);
            }
            # Why ? Because it yields before other processes, thereby
            # improving fairness of the scheduling. It's mild, but
            # otherwise we would be too aggressive in grabbing the server
            # attention for as long as _we_ like.
            sleep 0.25;
        }
        exit 0;
   }
   $ch{$child}=$h;
}

my $remaining = scalar keys %ch;
my $finished = 0;

my $all_ok=1;

while ((my $pid=wait) > 0) {
    my $who = $ch{$pid};
    my $status = $?;
    $finished++;
    $remaining--;
    delete $ch{$pid};
    if ($status & 127) {
        printf "$who: died with signal %d, %s coredump\n",
        ($status & 127),  ($status & 128) ? 'with' : 'without';
        $all_ok=0;
    } elsif ($status >> 8) {
        printf "$who: status %d\n", ($status>>8);
        $all_ok=0;
    }
    my @rem = values %ch;
    @rem = map { s/^[\D-]+(\d+).*$/$1/; $_ } @rem;
    if (scalar keys %$dispatch > 1) {
        my $dt = sprintf('%.1f',  time-$start_time);
        print "$dt [$who] ready ($finished ready, $remaining remaining: ", join(' ', @rem), ")\n";
    }
}

exit ($all_ok ? 0 : 1);
