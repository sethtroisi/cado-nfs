#!/usr/bin/perl -w

# General script for factoring integers with Cado-NFS.
#
# Copyright 2008 Pierrick Gaudry, Emmanuel Thome, Paul Zimmermann, Jeremie Detrey
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



# usage:
#    cadofactor.pl param=<paramfile> wdir=<...> ...
# Parameters passed in arguments *after* param=... will override choices
# that are made in paramfile.

# See params.c59 for an example of parameter file.

# If the parameter n=<n> is given, then n is factored. Otherwise, it is
# taken from stdin.

# NB: all the shell commands that are run by this script are reported in
# the $wdir/$name.cmd file, with full arguments, so that it is easy to
# reproduce part of the computation.


# TODO-list:
#  - When we discover that we have enough relations, kill jobs everywhere
#    and import partial files (having more rels can not hurt)
#  - Do not restart sieving tasks if we happen to have enough relations
#  - Create some "scan-holes" mechanism to fill in the blanks.
#  - Bench polynomials with some sieve before selecting the best one.
#  - Enable a 'lowmem' option
#  - Use compressed files (cf duplicates -out nodup.gz). Should we implement
#    such a mechanism in las, or should we simply use las | gzip?
#  - Recognize when the linear algebra was already started, and in such a case
#    don't restart the sieving/pruning/filtering/merging steps, but resume
#    bw with resume=1.
#  - Use free relations (they are computed, but never used...)

use strict;
use warnings;

use Data::Dumper;

use File::Copy;
use File::Basename;
use Cwd qw(abs_path);
use POSIX qw(ceil floor);

# Default parameters

# This list gives:
#  - The preferred ordering for parameters.
#  - The default values (if any).
my @parameter_defaults = (
    # global
    wdir=>undef,
    cadodir=>undef,
    name=>undef,
    machines=>undef,
    n=>undef,
    parallel=>0,

    # polyselect using Kleinjung (kj) algorithm
    degree=>5,
    kjkeep=>100,
    kjkmax=>10,
    kjincr=>60,
    kjl=>7,
    kjM=>1e25,
    kjpb=>256,
    kjp0max=>100000,
    kjadmin=>undef,
    kjadmax=>undef,
    kjadrange=>1e7,
    kjdelay=>120,
    selectnice=>10,

    # sieve
    rlim=>8000000,
    alim=>8000000,
    lpbr=>29,
    lpba=>29,
    mfbr=>58,
    mfba=>58,
    rlambda=>2.3,
    alambda=>2.3,
    I=>13,
    excess=>100,
    qmin=>12000000,
    qrange=>1000000,
    checkrange=>1000000,

    delay=>120,
    sievenice=>19,
    keeprelfiles=>0,

    # filtering
    prune=>1.0,
    keep=>160, # should be 128+skip
    keeppurge=>100000,
    maxlevel=>15,
    cwmax=>200,
    rwmax=>200,
    ratio=>1.5,
    bwstrat=>1,
    skip=>32,

    # linalg
    linalg=>'bw',
    bwmt=>2,

    # characters
    nkermax=>30,
    nchar=>50,
);

my @official_param_list=();
# Build the ordered list of parameters.
for(my $i = 0 ; $i < $#parameter_defaults ; $i+=2) {
    push @official_param_list, $parameter_defaults[$i];
}

my $default_param = {};
# Build the list of default parameters
%$default_param = @parameter_defaults;

sub read_param {
    my @args = @_;
    my $param = $default_param;
    READ_ARGS: while (defined($_=shift @args)) {
        for my $p (@official_param_list) {
            if (/^$p=(.*)$/) { $param->{$p}=$1; next READ_ARGS; }
        }
        if (/^params?=(.*)$/) {
            my $file = $1;
            open FILE, "<$file" or die "$file: $!\n";
            my @newargs;
            my $line;
            while (<FILE>) {
                $line = $_;
                if (/^\s*#/) { next; }
                if (/^\s*$/) { next; }
                if (/^\s*(\w+)=([\w\.\/\$]*)\s*(?:#.*)?$/) {
                    push @newargs, "$1=$2";
                    next;
                }
                die "Could not parse line: $_ in file $file\n";
            }
            close FILE;
            unshift @args, @newargs;
            next;
        }
        if (-f $_) { unshift @args, "param=$_"; next; }
        die "Unknown argument: $_\n";
    }
    # sanity check ; old config files may still have true/false values
    while (my ($k,$v) = each %$param) {
        next unless defined $v;
        if ($v eq 'false') { $param->{$k}=0; }
        if ($v eq 'true') { $param->{$k}=1; }
    }
    # checking mandatory parameters:
    if (!$param->{'wdir'}) { die "I need a wdir argument!\n"; }
    if (!$param->{'name'}) { die "I need a name argument!\n"; }
    if (!$param->{'machines'} && $param->{'parallel'}) {
        die "With parallel, I need a machines argument!\n"; }
    if (!$param->{'cadodir'}) {
        print STDERR "Warning: taking current script's basedir as cadodir\n";
        $param->{'cadodir'} = abs_path(dirname($0));
    }
    return $param;
}

sub print_param {
    my $param = shift @_;
    my $fh = *STDOUT{IO};
    if (scalar @_) {
        $fh = shift @_;
    }
    for my $k (@official_param_list) {
        if (my $v = $param->{$k}) {
            print $fh "$k=$v\n";
        }
    }
}

# This is a filename, or possibly undef.
my $cmdlog;

# Run $cmd
# In case the return code of the command is not 0, we die, unless 
# a 2rd argument 'no_kill' is passed.
sub my_system {
    my $cmd = shift @_;
    my $nokill=0;
    if ($_[0] && $_[0] eq 'no_kill') {
        $nokill=1;
    }
    system "echo $cmd >> $cmdlog" if $cmdlog;
    my $ret = system($cmd);
    if ($ret != 0 && !$nokill) {
        die "$cmd exited unexpectedly: $!\n";
    }
}

# Execute the given $cmd using backticks, allowing $timeout time for it
# to finish. 
# Returns a triple ( $t, $ret, $status )
#   where $t is a 0 or 1 depending whether the timout was reached (0
#   means timeout, 1 means success),
#   $ret is the returned value (empty string if timeout),
#   and $status is the exit status of the command (-1 if timeout).
#
# For timeout, we use the alarm mechanism. See perldoc -f alarm.
sub my_system_timeout {
    my ($cmd, $timeout) = @_;
    my $ret = "";
    my $status = 0;
    eval {
        local $SIG{'ALRM'} = sub { die "alarm\n" }; # NB: \n required
        alarm $timeout;
        $ret = `$cmd`;
        $status = $? >> 8;
        alarm 0;
    };
    if ($@) {
        die unless $@ eq "alarm\n";   # propagate unexpected errors
        print "WARNING! The command $cmd timed out.\n";
        return (0, "", -1);
    }
    else {
        return (1, $ret, $status);
    }
}

sub banner {
    print STDERR "#" x 70, "\n";
    print STDERR "## $_[0]\n";
}

# returns -1, 0 or 1, depending is the last modification date of 
# $f1 is less, equal or more than the modif date of $f2.
sub cmp_filedate {
    my ($f1, $f2) = @_;
    my $d1 = (stat($f1))[9];
    my $d2 = (stat($f2))[9];
    if ($d1 < $d2) {
        return -1;
    } elsif ($d1 == $d2) {
        return 0;
    } else {
        return 1;
    }
}



###########################################################################
### Some functions to manipulate a list of ranges and detect holes
###########################################################################

# Invariant: A good-formed list of ranges looks like
#     [ [a0,b0] , [a1,b1], ...   [ak,bk] ]
# with bi < a_{i+1}   (and of course ai < bi)
# All the inegalities are *strict* (otherwise, merge ranges!)

sub check_wellformed_rangelist {
    my $T = shift @_;
    my $curr_value = -1;

    foreach my $range (@$T) {
        my ($a, $b) = @$range;
        if ($a >= $b) { return 0; }
        if ($a <= $curr_value) { return 0; }
        $curr_value = $b;
    }
    return 1;
}


# Take a list of ranges and returns a good-formed list
sub filter_ranges {
    my @ranges = @_;

    @ranges = sort { $a->[0] <=> $b->[0] } @ranges;

    my @brokendown = ();
    my @overlap = ();

    my $current = shift @ranges or return @ranges;
    my ($cs,$ce) = @$current;

    while (defined(my $x = shift @ranges)) {
        my ($s, $e) = @$x;
        if ($e == $s) {
            print STDERR "stupid sub-range [${s}..${e}[\n";
            next;
        }
        die "Uh : $e <= $s" if $e <= $s;

        if ($s == $ce) {
            $ce = $e;
            @overlap=($x);
        } elsif ($s < $ce) {
            # fetch the viewed ranges that overlap with this range.
            my @y=();
            for my $z (@overlap) {
                push @y, $z if $z->[1] > $s;
                my ($os, $oe) = @{$z};
            }
            @overlap = @y;
            push @overlap, $x;
            @overlap = sort { $a->[1] <=> $b->[1] } @overlap;
            if ($e > $ce) {
                $ce = $e;
            }
        } else {
            push @brokendown, [$cs, $ce];
            $cs = $s;
            $ce = $e;
            @overlap=($x);
        }
    }
    push @brokendown, [$cs, $ce];
    return @brokendown;
}


# Return first hole in list of ranges, with given rangemax, and amin
# The list of ranges is updated.
sub get_next_hole {
  my $ranges = shift @_;
  my $rangemax = shift @_;
  my $amin = shift @_;

  if (scalar @$ranges == 0) {
      my @r = ($amin, $amin+$rangemax);
      return \@r, [\@r];
  }

  my $r0 = $$ranges[0];
  # Case where there is some room at the begining
  if ($$r0[0] > $amin) {
      my $x = $amin+$rangemax;
      if ($x >= $$r0[0]) {
          $$ranges[0] = [$amin, $$r0[1]];
          my @r = ($amin, $$r0[0]);
          return \@r, $ranges;
      } else {
          my @r = ($amin, $x);
          unshift @$ranges, \@r;
          return \@r, $ranges;
      }
  }
  # Case where there is no second range
  if (scalar @$ranges == 1) {
      my $x = $$r0[1];
      my $y = $x+$rangemax;
      my @r = ($x, $y);
      $$ranges[0] = [$$r0[0], $y];
      return \@r, $ranges;
  }
  # last case: a hole between first and second range
  my $r1 = $$ranges[1];
  my $x = $$r0[1];
  my $y = $x + $rangemax;
  if ($y >= $$r1[0]) { 
      # the hole will be filled in
      shift @$ranges;
      $$ranges[0] = [$$r0[0], $$r1[1]];
      my @r = ($$r0[1], $$r1[0]);
      return \@r, $ranges;
  } else {
      $$ranges[0] = [$$r0[0], $y];
      my @r = ($x, $y);
      return \@r, $ranges;
  }
}
       
###########################################################################
## Parallel polynomial selection
###########################################################################

# Read a selection status file:
# Each line is of the form
#   hostname admin admax
# Any empty line is ignored.
# any line that starts with # is ignored
sub read_select_status_file {
    my $param = shift @_;
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $statfile = "$prefix.selectstatus";
    if (! -f $statfile) {
        print "No status file found. Creating an empty one...\n";
        open GRR, ">$statfile" or die "$statfile: $!\n";
        close GRR;
        return ();
    }
    my @status;
    open SF, "<$statfile" or die "$statfile: $!\n";
    while (<SF>) {
        if (/^\s*\n$/) { next; }
        if (/^#/) { next; }
        if (/^\s*([\w\.]*)\s+([\w\.]+)\s+([\w\.]+)(\s+([\w\.]+))?/) {
            push @status, [$1, $2, $3, $5];
            next;
        }
        die "Could not parse status file: $_\n";
    }
    close SF;
    return @status;
}

sub write_select_status_file {
    my $param = shift @_;
    my $status = shift @_;
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $statfile = "$prefix.selectstatus";
    open GRR, ">$statfile" or die "$statfile: $!\n";
    foreach my $t (@$status) {
        my @tt = @$t;
        print GRR "" . $tt[0] . " " . $tt[1] . " " . $tt[2] . " " . $tt[3] . "\n";
    }
    close GRR;
}

# A task is a quadruple (hostname admin admax pid)
# This function ssh to hostname, and tests if the task is still alive,
# dead or finished.
# 1 means running
# 0 means finished
# -1 means (probably dead)
# -2 means binary for polynomial selection was not found
sub check_running_select_task {
    my $param = shift @_;
    my $name = $param->{'name'};
    my $mach_desc = shift @_;
    my ($host, $admin, $admax, $pid) = @_;
    print "    $host \t$admin\t$admax:\n";
    my %host_data=%{$mach_desc->{$host}};
    my $wdir = $host_data{'tmpdir'};
    my $cadodir = $host_data{'cadodir'};

    # First check if the last line corresponds to finished job:
    my ($t, $lastline) = my_system_timeout(
        "ssh $host tail -1 $wdir/$name.kjout.$admin-$admax", 30);
    if (! $t) {
        print "      Assume task is still running...\n";
        return 1;
    }
    $_ = $lastline;
    if (/# generated/) {
        print "      finished!\n";
        return 0;
    }
    if (/No polynomial found/) {
        print "      finished but found no polynomial!\n";
        return 0;
    }
    if (/No such file or directory/) {
        print "      program for polynomial selection was not found!\n";
        return -2;
    }

    ## If file is partial, check its last modification time:
    #my $date;
    #($t, $date) = my_system_timeout("ssh $host date +%s", 30);
    #if (! $t) {
    #    print "      Assume task is still running...\n";
    #    return 1;
    #}
    #my $modifdate;
    #($t, $modifdate) = my_system_timeout(
    #    "ssh $host stat -c %Y $wdir/$name.kjout.$admin-$admax", 30);
    #if (! $t) {
    #    print "      Assume task is still running...\n";
    #    return 1;
    #}
    ## If didn't move for 10 minutes, assume it's dead
    #if ($date > ($modifdate + 600)) {
    #    print "      dead ?!?\n";
    #    return -1;
    #}
    #

    # TODO: this is still not good enough: if the process dies and another
    # one jumps in and takes its PID, we won't detect it.
    # Possible solutions: check the /proc/$pid/cmdline (Linux-specific)
    # or keep a temporary file regularly touched by the program.

    if ($pid) {
        my $ret;
        my $status;
        ($t, $ret, $status) = my_system_timeout(
            "ssh $host \"/bin/kill -0 $pid\"", 30);
        if (! $t) {
            print "      Assume task is still running...\n";
            return 1;
        }
        if ($status != 0) {
            print "      dead ?!?\n";
            return -1;
        }
    }

    # otherwise it's running:
    print "      running...\n";
    return 1;
}

sub push_select_files {
    my ($mach, $wdir, $param) = @_;
    my $name = $param->{'name'};
    my $ldir = $param->{'wdir'};
    my $t;
    my $ret;
    my $status;
    ($t, $ret) = my_system_timeout("ssh $mach mkdir -p $wdir", 60);
    if (! $t) { die "Connection to $mach timeout\n"; }
    ($t, $ret, $status) = my_system_timeout("ssh $mach test -e $wdir/$name.n", 120);
    if (! $t) { die "Connection to $mach timeout\n"; }
    if ($status != 0) {
        print "    Pushing $name.n to $mach.\n";
        system("rsync --timeout=120 $ldir/$name.n $mach:$wdir/");
    }
}

sub restart_select_tasks {
    my $param = shift @_;
    my $running = shift @_;
    my $mach_desc = shift @_;
    my $name = $param->{'name'};
    my $nice = $param->{'selectnice'};

    my @ranges = ();
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $wdir = $param->{'wdir'};
    opendir(DIR, $wdir) or die "can't opendir $wdir: $!";
    my @files = readdir(DIR);
    close DIR;
    foreach my $f (@files) {
        if ($f =~ /$name\.kjout\.([\w\.]+)-([\w\.]+)/) {
            push @ranges, [$1, $2];
        }
    }
    foreach my $t (@$running) {
        push @ranges, [${$t}[1], ${$t}[2]];
    }

    @ranges = filter_ranges(@ranges);
    die unless check_wellformed_rangelist(\@ranges);

    ### Restart as many tasks as neeeded.
    for my $m (keys %$mach_desc) {
        my $cpt = 0;
        my $t;
        my %desc = %{$mach_desc->{$m}};
        for $t (@$running) {
            if (${$t}[0] eq $m) { $cpt++; }
        }
        if ($cpt >= $desc{'cores'}) {
            # already enough tasks running on $m
            next;
        }
        # Number of new tasks to run is:
        $cpt = $desc{'cores'} - $cpt;
        for (my $i = 0; $i < $cpt; $i++) {
            ## Get next task to start
            my ($xx, $yy) = get_next_hole(\@ranges, $param->{'kjadrange'},
                $param->{'kjadmin'});
            my ($admin, $admax) = @$xx;
            if ($admin >= $param->{'kjadmax'}) {
                return;
            }
            if ($admax >= $param->{'kjadmax'}) {
                $admax = $param->{'kjadmax'};
            }
            @ranges = @$yy;
            my $mach = $m;
            my $wdir = $desc{'tmpdir'};
            my $bindir = $desc{'cadodir'};
            push_select_files($mach, $wdir, $param);
            my $cmd = "nice -$nice $bindir/polyselect/polyselect" .
              " -keep " . $param->{'kjkeep'} .
              " -kmax " . $param->{'kjkmax'} .
              " -incr " . $param->{'kjincr'} .
              " -l " . $param->{'kjl'} .
              " -M " . $param->{'kjM'} .
              " -pb " . $param->{'kjpb'} .
              " -p0max " . $param->{'kjp0max'} .
              " -admin " . $admin .
              " -admax " . $admax .
              " -degree " . $param->{'degree'} .
              " < $wdir/$name.n";
            my $outfile = "$wdir/$name.kjout.$admin-$admax";
            my $pid = `ssh $mach "/bin/sh -c '$cmd > $outfile 2>&1& echo \\\$!'"`;
            chomp $pid;
            print "    Starting $mach $admin $admax.\n";
            open FH, ">> " . $param->{'wdir'} . "/$name.cmd";
            print FH "ssh $mach \"/bin/sh -c '$cmd > $outfile 2>&1& echo \\\$!'\"\n";
            close FH;
            push @$running, [$m, $admin, $admax, $pid];
        }
    }
}

sub get_logmualpha_value {
    my $filename = shift @_;
    open FILE, "$filename";
    my @lines = readline(FILE);
    close FILE;
    @lines = grep (/E=/, @lines);
    $_ = $lines[$#lines];
    if (/E=(\d+\.\d*)/) {
        my $E = $1;
        return $E;
    } else {
        print STDERR "Can not parse output of $filename for geting logmu+alpha\n";
        print STDERR "Putting an arbitrary (large) value for E\n";
        return 100000;
    }
}


# rsync a finished / dead task to the local working dir.
# return the logmu+alpha value.
sub import_select_task_result {
    my $param = shift @_;
    my $name = $param->{'name'};
    my $mach_desc = shift @_;
    my ($host, $admin, $admax) = @_;
    my %host_data=%{$mach_desc->{$host}};
    my $rwdir = $host_data{'tmpdir'};
    my $lwdir = $param->{'wdir'};
    my $outfile = "$name.kjout.$admin-$admax";
    my $cmd = "rsync --timeout=30 $host:$rwdir/$outfile $lwdir/";
    my $ret = system($cmd);
    if ($ret != 0) {
        print STDERR "Problem when importing file $outfile from $host.\n";
        return 0;
    }
    my_system_timeout("ssh $host /bin/rm $rwdir/$outfile", 30);

    return get_logmualpha_value("$lwdir/$outfile");
}


sub parallel_polyselect {
    my $param = shift @_;
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $name = $param->{'name'};
    my $t;
    my $finished = 0;
    while (! $finished) {
        my %mach_desc = read_machine_description($param);
        my @status = read_select_status_file($param);
        my @newstatus;
        print "  Check all running tasks:\n";
        foreach $t (@status) {
            my $running = check_running_select_task($param, \%mach_desc, @$t);
            if ($running == 1) {
                push @newstatus, $t;
            } elsif ($running == 0) {
                import_select_task_result($param, \%mach_desc, @$t);
            } elsif ($running == -2) {
                die "program for polynomial selection not found\n";
            } else { # dead task!
                # do nothing, this will be restarted soon...
            }
        }
        print "  Start new tasks.\n";
        restart_select_tasks($param, \@newstatus, \%mach_desc);
        write_select_status_file($param, \@newstatus);
        # check if finished
        my @ranges = ();
        my $wdir = $param->{'wdir'};
        opendir(DIR, $wdir) or die "can't opendir $wdir: $!";
        my @files = readdir(DIR);
        close DIR;
        foreach my $f (@files) {
            if ($f =~ /$name\.kjout\.([\w\.]+)-([\w\.]+)/) {
                push @ranges, [$1, $2];
            }
        }
        @ranges = filter_ranges(@ranges);
        die unless check_wellformed_rangelist(\@ranges);
        if (scalar @ranges > 0) {
            my $first_range=$ranges[0];
            my ($a, $b) = @$first_range;
            if ($a <= $param->{'kjadmin'} && $b >= $param->{'kjadmax'}) {
                $finished = 1;
            }
        }

        if (!$finished) {
            my $delay = $param->{'kjdelay'};
            print "Wait for $delay seconds before checking again.\n";
            sleep($delay);
        }
    }

    # A bit of cleaning on slaves
    my %mach_desc = read_machine_description($param);
    for my $m (keys %mach_desc) {
        my %desc = %{$mach_desc{$m}};
        my $wdir = $desc{'tmpdir'};
        my_system_timeout("ssh $m /bin/rm -f $wdir/$name.n", 30);
    }


    # Choose best according to logmu+alpha
    my $Emin;
    my $bestfile;

    my $wdir = $param->{'wdir'};
    opendir(DIR, $wdir) or die "can't opendir $wdir: $!";
    my @files = readdir(DIR);
    close DIR;
    foreach my $f (@files) {
        if ($f =~ /$name\.kjout\.([\w\.]+)-([\w\.]+)/) {
            my $E = get_logmualpha_value("$wdir/$f");
            if ((!$Emin) || $E < $Emin) {
                $Emin = $E;
                $bestfile = "$wdir/$f";
            }
        }
    }
    print "Best polynomial is from $bestfile\n";
    copy("$bestfile" , "$prefix.poly");
    system("echo cp $bestfile $prefix.poly >> $prefix.cmd");
}



sub polyselect {
    my $param = shift @_;
    banner("polynomial selection");

    if ($param->{'parallel'}) {
        parallel_polyselect($param, @_);
        return;
    }
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $b;
    my $bestb=0;
    my $Emin;
    for ($b = $param->{'bmin'}; $b <= $param->{'bmax'}; $b++) {
        print "Running polynomial selection with b=$b...\n";
        my $effort=$param->{'e'};
        my $degree=$param->{'degree'};
        if (-f "$prefix.poly.$b") {
            print "Result file is already there. Skip the computation!\n";
        } else {
            my $cmd = $param->{'cadodir'}."/polyselect/polyselect " .
              "-b $b -e $effort -degree $degree < $prefix.n";
            my_system "$cmd > $prefix.poly.$b";
        }

        my $E = get_logmualpha_value "$prefix.poly.$b";
        if ((!$Emin) || $E < $Emin) {
            $Emin = $E;
            $bestb = $b;
        }
    }
    # choose best.
    # For the moment, we do it according to logmu+alpha
    print "Best polynomial is for b = $bestb\n";
    copy("$prefix.poly.$bestb" , "$prefix.poly");
    system("echo cp $prefix.poly.$bestb $prefix.poly >> $prefix.cmd");
}

###########################################################################
## Parallel sieving
###########################################################################

sub try_singleton {
    my $param = shift @_;
    my $nrels = shift @_;

    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};

    my @files = @_;
    my $filestring = "$prefix.rels*";
    if (scalar @files) {
        $filestring = join(' ', @files);
    }

    banner("Removing duplicates");
    my $cmd = $param->{'cadodir'} .
      "/linalg/duplicates -nrels $nrels -out $prefix.nodup $filestring 2> $prefix.duplicates.stderr";
    my_system $cmd;
    my @grouik = split(/ /, `tail -1 $prefix.duplicates.stderr`);
    my $nrels_dup = $grouik[(scalar @grouik)-1];
    chomp $nrels_dup;
    print "We have $nrels_dup relations left after removing duplicates\n";

    banner("Singleton removal");
    my $keep = $param->{'keeppurge'};
    $cmd = $param->{'cadodir'} .
      "/linalg/purge -poly $prefix.poly -keep $keep -nrels $nrels_dup -purged $prefix.purged $prefix.nodup  2> $prefix.purge.stderr";
    my_system $cmd, "no_kill";
    if (-z "$prefix.purged" || ! -f "$prefix.purged") {
        printf `grep "expected" $prefix.purge.stderr`;
        printf "Not enough relations!\n";
        return 0;
    } else {
        open FH, "< $prefix.purged";
        my $line = <FH>;
        my ($nrows, $ncols) = split(/ /, $line);
        close FH;
        print "Nrows = $nrows , Ncols = $ncols";
        print "Excess = " . ($nrows - $ncols) . "\n";
        if ($nrows - $ncols <= $param->{'excess'}) {
            print "Not enough relations!\n";
            return 0;
        }
    }
    return 1;
}

# Read a status file:
# Each line is of the form
#   hostname q0 q1
# Any empty line is ignored.
# any line that starts with # is ignored
sub read_status_file {
    my $param = shift @_;
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $statfile = "$prefix.status";
    if (! -f $statfile) {
        print "No status file found. Creating an empty one...\n";
        open GRR, ">$statfile" or die "$statfile: $!\n";
        close GRR;
        return ();
    }
    my @status;
    open SF, "<$statfile" or die "$statfile: $!\n";
    while (<SF>) {
        if (/^\s*\n$/) { next; }
        if (/^#/) { next; }
        if (/^\s*(\w*)\s+(\d+)\s+(\d+)(\s+(\d+))?/) {
            push @status, [$1, $2, $3, $5];
            next;
        }
        die "Could not parse status file: $_\n";
    }
    close SF;
    return @status;
}

sub write_status_file {
    my $param = shift @_;
    my $status = shift @_;
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $statfile = "$prefix.status";
    open GRR, ">$statfile" or die "$statfile: $!\n";
    foreach my $t (@$status) {
        my @tt = @$t;
        print GRR "" . $tt[0] . " " . $tt[1] . " " . $tt[2] . " " . $tt[3] . "\n";
    }
    close GRR;
}


# A task is a quadruple (hostname q0 q1 pid)
# This function ssh to hostname, and tests if the task is still alive,
# dead or finished.
sub check_running_task {
    my $param = shift @_;
    my $name = $param->{'name'};
    my $mach_desc = shift @_;
    my ($host, $q0, $q1, $pid) = @_;
    print "    $host \t$q0-$q1:\n";
    my %host_data=%{$mach_desc->{$host}};
    my $wdir = $host_data{'tmpdir'};
    my $cadodir = $host_data{'cadodir'};
    # First check if the last line corresponds to finished job:
    my ($t, $lastline) = my_system_timeout(
        "ssh $host tail -1 $wdir/$name.rels.$q0-$q1", 30);
    if (! $t) {
        print "      Assume task is still running...\n";
        return 1;
    }
    $_ = $lastline;
    if (/# Total/) {
        print "      finished!\n";
        return 0;
    }
    # If file is partial check its last modification time:
#    my $date;
#    ($t, $date) = my_system_timeout("ssh $host date +%s", 30);
#    if (! $t) {
#        print "      Assume task is still running...\n";
#        return 1;
#    }
#    my $modifdate;
#    ($t, $modifdate) = my_system_timeout(
#        "ssh $host stat -c %Y $wdir/$name.rels.$q0-$q1", 30);
#    if (! $t) {
#        print "      Assume task is still running...\n";
#        return 1;
#    }
#    ## If didn't move for 10 minutes, assume it's dead
#    if ($date > ($modifdate + 600)) {
#        print "      dead ?!?\n";
#        return 0;
#    }

    # TODO: this is still not good enough: if the process dies and another
    # one jumps in and takes its PID, we won't detect it.
    # Possible solutions: check the /proc/$pid/cmdline (Linux-specific)
    # or keep a temporary file regularly touched by the program.
    if ($pid) {
        my $ret;
        my $status;
        ($t, $ret, $status) = my_system_timeout(
            "ssh $host \"/bin/kill -0 $pid\"", 30);
        if (! $t) {
            print "      Assume task is still running...\n";
            return 1;
        }
        if ($status != 0) {
            print "      dead ?!?\n";
            # TODO: this is no good! We should return a different status code
            # to distinguish the "dead" status from the "finished" status.
            return 0;
        }
    }

    # otherwise it's running:
    print "      running...\n";
    return 1;
}


#
# rsync a finished / dead task to the local working dir.
# Return the number of relations in the imported file.
sub import_task_result {
    my $param = shift @_;
    my $name = $param->{'name'};
    my $mach_desc = shift @_;
    my ($host, $q0, $q1) = @_;
    my %host_data=%{$mach_desc->{$host}};
    my $rwdir = $host_data{'tmpdir'};
    my $lwdir = $param->{'wdir'};
    my $fname="$name.rels.$q0-$q1";
    my $lfname="$lwdir/$fname";
    my $rfname="$rwdir/$fname";
    my $cmd = "rsync --timeout=30 $host:$rfname $lwdir/";
    my $ret = system($cmd);

    if ($ret != 0) {
        print STDERR "Problem when importing file $fname from $host.\n";
        return 0;
    }
    $cmd = $param->{'cadodir'} . "/utils/check_rels -poly $lwdir/$name.poly $lfname";
    $ret = system($cmd);
    if ($ret != 0) {
        print STDERR "Buggy relation in file $fname from $host.\n";
        print STDERR "Moving it to FIXME.$fname .\n";
        print STDERR "If you fix it, should remove the prefix and it will "
          . "be taken into account in the computation";
        rename "$lfname", "$lwdir/FIXME.$fname"
            or die "Failed rename: $!\n";
        return 0;
    }
    # Import succeeded, so we can remove the remote file.
    unless ($param->{'keeprelfiles'}) {
        my_system_timeout("ssh $host /bin/rm $rfname", 30);
    }

    my $n = `grep -v "^#" $lfname | wc -l`;
    chomp($n);
    return ($n, $lfname);
}

sub read_machine_description {
    my $param = shift @_;
    my $machine_file = $param->{'machines'};
    my %mach_desc;
    my %vars = ();
    open MACH, "< $machine_file";
    while (<MACH>) {
        next if /^\s*#/ || /^\s*$/;

        s/^\s*//;

        if (/^(\w+)=(\S*)\s*$/) {
            $vars{$1}=$2;
            next;
        }

        if (/^\[(\S+)\]\s*$/) {
            %vars = ( cluster=>$1 );
            next;
        }

        if (s/^([\w\.]+)\s*(.*)/$2/) {
            my $m = $1;
            my $opts = $2;
            my %h = %vars;
            while ($opts =~ s/^(\w+)=(\S*)\s*//) {
                $h{$1}=$2;
            }
            die "$opts" if $opts;
            # Check mandatory args and complete with defaults.
            if (! $h{'tmpdir'}) { die "No tmpdir given for machine $m\n"; }
            if (! $h{'cadodir'}) { die "No cadodir given for machine $m\n"; }
            if (! defined($h{'cores'})) { $h{'cores'}=1; }

            # Substitute $name if %s is found.
            if ($h{'tmpdir'} =~ /%s/) {
                $h{'tmpdir'}=sprintf($h{'tmpdir'}, $param->{'name'});
            }

            $mach_desc{$m}=\%h;
            next;
        }

        die "Could not parse: $_ in file $machine_file\n";
    }
    close MACH;

    return %mach_desc;
}



# Check whether .poly and .roots files are present in the remote working
# directory. Otherwise push them.
sub push_files {
    my ($mach, $wdir, $param) = @_;
    my $name = $param->{'name'};
    my $ldir = $param->{'wdir'};
    my $t;
    my $ret;
    my $status;
    ($t, $ret) = my_system_timeout("ssh $mach mkdir -p $wdir", 60);
    if (! $t) { die "Connection to $mach timeout\n"; }
    ($t, $ret, $status) = my_system_timeout("ssh $mach test -e $wdir/$name.poly", 120);
    if (! $t) { die "Connection to $mach timeout\n"; }
    if ($status != 0) {
        print "    Pushing $name.poly to $mach.\n";
        system("rsync --timeout=120 $ldir/$name.poly $mach:$wdir/");
    }
    ($t, $ret, $status) = my_system_timeout("ssh $mach test -e $wdir/$name.roots", 120);
    if (! $t) { die "Connection to $mach timeout\n"; }
    if ($status != 0) {
        print "    Pushing $name.roots to $mach.\n";
        system("rsync --timeout=120 $ldir/$name.roots $mach:$wdir/");
    }
}    


# Look for next special q to be launched.
# This takes into account filenames in the working directory and running
# tasks as given.
sub get_next_q {
    my ($param, $running) = @_;
    my $maxq = 0;
    my $wdir = $param->{'wdir'};
    my $name = $param->{'name'};
    opendir(DIR, $wdir) or die "can't opendir $wdir: $!";
    my @files = readdir(DIR);
    close DIR;
    foreach my $f (@files) {
        if ($f =~ /$name\.rels\.(\d+)-(\d+)/) {
            if ($2 > $maxq) {
                $maxq = $2;
            }
        }
    }
    foreach my $t (@$running) {
        if (${$t}[2] > $maxq) {
            $maxq = ${$t}[2];
        }
    }
    if ($maxq < $param->{'qmin'}) {
        $maxq = $param->{'qmin'};
    }
    return $maxq;
}

sub restart_tasks {
    my $param = shift @_;
    my $running = shift @_;
    my $mach_desc = shift @_;
    my $name = $param->{'name'};
    my $nice = $param->{'sievenice'};

    my $q_curr = get_next_q($param, $running);

    for my $m (keys %$mach_desc) {
        my $cpt = 0;
        my $t;
        my %desc = %{$mach_desc->{$m}};
        for $t (@$running) {
            if (${$t}[0] eq $m) { $cpt++; }
        }
        if ($cpt >= $desc{'cores'}) {
            # already enough tasks running on $m
            next;
        }
        # Number of new tasks to run is:
        $cpt = $desc{'cores'} - $cpt;
        for (my $i = 0; $i < $cpt; $i++) {
            my $qend = $q_curr + $param->{'qrange'};
            my $mach = $m;
            my $wdir = $desc{'tmpdir'};
            my $bindir = $desc{'cadodir'};
            push_files($mach, $wdir, $param);
            my $cmd = "nice -$nice $bindir/sieve/las" .
              " -I " . $param->{'I'} .
              " -poly $wdir/$name.poly" .
              " -fb $wdir/$name.roots" .
              " -q0 $q_curr -q1 $qend";
            my $pid = `ssh $mach "/bin/sh -c '$cmd > $wdir/$name.rels.$q_curr-$qend 2>&1& echo \\\$!'"`;
            chomp $pid;
            print "    Starting $mach $q_curr-$qend.\n";
            open FH, ">> " . $param->{'wdir'} . "/$name.cmd";
            print FH "ssh $mach \"/bin/sh -c '$cmd > $wdir/$name.rels.$q_curr-$qend 2>&1& echo \\\$!'\"\n";
            close FH;
            push @$running, [$m, $q_curr, $qend, $pid];
            $q_curr = $qend;
        }
    }
}



# returns whether some new completed file has been imported.
sub parallel_sieve_update {
    my $param = shift @_;
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    print "  Reading machine description file.\n";
    my %mach_desc = read_machine_description($param);
    print "  Reading status file.\n";
    my @status = read_status_file($param);
    my $new_rels=0;

    my $t;
    my @newstatus;
    my @files=();
    # Loop over all tasks that are listed in status file and
    # check whether they are finished.
    print "  Check all running tasks:\n";
    foreach $t (@status) {
        my $running = check_running_task($param, \%mach_desc, @$t);
        if ($running) {
            push @newstatus, $t;
        } else {
            # Task is finished or dead. 
            my ($trels,@tfiles) = import_task_result($param, \%mach_desc, @$t);
            $new_rels += $trels;
            push @files, @tfiles;
        }
    }
    # If some tasks have finished, restart them according to mach_desc
    # and try singleton removal
    print "  Restarting tasks if necessary.\n";
    restart_tasks($param, \@newstatus, \%mach_desc);
    print "  Writing new status file.\n";
    write_status_file($param, \@newstatus);
    return $new_rels, @files;
}

sub count_rels {
    my $param = shift @_;
    my $name = $param->{'name'};
    my $wdir = $param->{'wdir'};
    my $nrels = 0;
    opendir(DIR, $wdir) or return 0;
    my @files = readdir(DIR);
    close DIR;
    my @relfiles=();
    foreach my $f (@files) {
        if ($f =~ /$name\.rels\.(\d+)-(\d+)/ || $f =~ /$name\.freerels/) {
            $nrels+= `grep -v "^#" $wdir/$f | wc -l`;
            push @relfiles, "$wdir/$f";
        }
    }
    return $nrels, @relfiles;
}
 

sub parallel_sieve {
    my $param = shift @_;
    my $finished = 0;
    my ($nrels,@files) = count_rels($param);
    print "We have found $nrels relations in working dir\n";
    print "Let's start the main loop!\n";
    my $prev_check = 0;
    while (! $finished) {
        print "Check what's going on on different machines...\n";
        my ($new_rels,@newfiles) = parallel_sieve_update($param);
        if ($new_rels) {
            $nrels += $new_rels;
            push @files, @newfiles;
            print "We have now $nrels relations\n";
            if ($nrels-$prev_check > $param->{'checkrange'}) {
                print "Trying singleton removal...\n";
                $finished = try_singleton($param, $nrels, @files);
                $prev_check = $nrels;
            }
        }
        if (! $finished) { 
            my $delay = $param->{'delay'};
            print "Number of relations is $nrels.\n";
            print "Wait for $delay seconds before checking again.\n";
            sleep($delay);
        }
    }

    ## kill remaining siever processes
    my @status = read_status_file($param);
    foreach my $t (@status) {
        my $m   = ${$t}[0];
        my $pid = ${$t}[3];
        my $ret = `ssh $m "/bin/kill -KILL $pid"` if $pid;
    }

    ## clean remaining sieving-related files on slaves
    my %mach_desc = read_machine_description($param);
    my $name = $param->{'name'};
    for my $m (keys %mach_desc) {
        my %desc = %{$mach_desc{$m}};
        my $wdir = $desc{'tmpdir'};
        my $ret = `ssh $m "/bin/rm $wdir/$name.rels.* $wdir/$name.poly $wdir/$name.roots"`;
    }
    return $nrels;
}

sub sieve {
    my $param = shift @_;
    banner("Sieve");
    if ($param->{'parallel'}) {
        return parallel_sieve($param, @_);
    }
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $finished=0;
    my $qcurr = $param->{'qmin'};
    my $nrels = 0;
    # taking (pi(2^lpba) + pi(2^lpbr)) / 3   as limit
    my $wantedrels = exp($param->{'lpba'}*log(2)) / ($param->{'lpba'}*log(2))
      + exp($param->{'lpbr'}*log(2)) / ($param->{'lpbr'}*log(2));
    $wantedrels = ceil($wantedrels/3);

    print "Sieving for $wantedrels relations\n";

    my %ls_result=();
    {
        opendir D, $param->{'wdir'};
        for my $x (readdir D) {
            next if $x =~ /^\./;
            next if $x =~ /(?:FIXME|bak|~)$/;
            $ls_result{$x}=1;
        }
        closedir D;
    }

    # Look for finished files.
    my @done=();

    {
        my $filename = "$prefix.rels";
        my $basename = "$param->{'name'}.rels";
        if (-f $filename) {
            print "Found file $filename...";
            my $lines = `grep -v "^#" $filename | wc -l`;
            chomp $lines;
            print "$lines rels\n";
            delete $ls_result{$basename};
            push @done, $filename;
            $nrels += $lines;
        }
    }

    while (1) {
        my $qend = $qcurr+$param->{'qrange'};
        my $filename = "$prefix.rels.$qcurr-$qend";
        my $basename = "$param->{'name'}.rels.$qcurr-$qend";
        last unless -f $filename;
        $qcurr = $qend;
        print "Found file $filename...";
        my $lines = `grep -v "^#" $filename | wc -l`;
        chomp $lines;
        print "$lines rels\n";
        delete $ls_result{$basename};
        push @done, $filename;
        $nrels += $lines;
    }

    my @missed = grep(/^$param->{'name'}.rels.\d+-\d+$/,keys %ls_result);
    print STDERR "Uh, I'm not considering these files:\n",
        join("\n",@missed), "\n";

    while (1) {
        if ($nrels >= $wantedrels) {
            last if try_singleton($param, $nrels, @done);
        }
        # Sieve another range.
        my $qend = $qcurr+$param->{'qrange'};
        my $filename = "$prefix.rels.$qcurr-$qend";
        my $las = "$param->{'cadodir'}/sieve/las";
        if (! -f $filename) {
            my $cmd = $las .
              " -I $param->{'I'}" .
              " -poly $prefix.poly" .
              " -fb $prefix.roots" .
              " -q0 $qcurr -q1 $qend";
            print "Running lattice siever for q in [$qcurr,$qend]...\n";
            my_system "/bin/sh -c '$cmd > $filename 2>&1'";
        } else {
            print "Found an existing relation file!\n";
        }
        print "Counting relations...";
        my $lines = `grep -v "^#" $filename | wc -l`;
        chomp $lines;
        print "$lines\n";
        push @done, $filename;
        $nrels += $lines;
        print "We have now $nrels relations; need $wantedrels\n";
        $qcurr=$qend;
    }
    return $nrels;
}

MAIN: {
    select(STDOUT); $|=1; # flush stdout always.

    my $param = read_param(@ARGV);

    my $wdir = $param->{'wdir'};

    if ($param->{'wdir'} =~ /%s/) {
        # substitute $name.
        $param->{'wdir'}=sprintf($param->{'wdir'}, $param->{'name'});
        $wdir = $param->{'wdir'};
    }

    # Create working directory if not there
    if (!-d $wdir) {
        mkdir $wdir or die "Cannot create $wdir: $!\n";
    }
    
    # Check if there is already some stuff relative to $name in $wdir
    # First thing is $name.n. If it is not there, we consider that
    # everything is obsolete, anyway.
    if (-f "$wdir/$param->{'name'}.n") {
        print STDERR "Warning: there is already some data relative to " .
        "this name in this directory.\n";
    }

    # Read n if not given on command line
    if (!$param->{'n'}) {
        print "'n' is not given in parameters, please enter the number to factor:\n";
        $param->{'n'} = <STDIN>;
        chomp ($param->{'n'});
    }

    # Create $name.n and $name.param in wdir.
    my $prefix = "$wdir/$param->{'name'}";
    open FN, ">$prefix.n";
    print FN "n:$param->{'n'}\n";
    close FN;
    my $fh;
    open $fh, ">$prefix.param";
    print_param($param,$fh);

    $cmdlog = "$prefix.cmd";

    # Polynomial selection
    if (! -f "$prefix.poly") {
        polyselect($param);
        # Appending some parameters to the poly file
        open FILE, ">> $prefix.poly";
        print FILE "rlim: $param->{'rlim'}\n";
        print FILE "alim: $param->{'alim'}\n";
        print FILE "lpbr: $param->{'lpbr'}\n";
        print FILE "lpba: $param->{'lpba'}\n";
        print FILE "mfbr: $param->{'mfbr'}\n";
        print FILE "mfba: $param->{'mfba'}\n";
        print FILE "rlambda: $param->{'rlambda'}\n";
        print FILE "alambda: $param->{'alambda'}\n";
    }

    # Creating the factor base
    my $cmd;
    if (-e "$prefix.roots" && 
        cmp_filedate("$prefix.poly", "$prefix.roots") < 1) {
            print "A factor base file exists... Let's trust it!\n";
    } else {
        banner("Factor base");
        my $makefb = "$param->{'cadodir'}/sieve/makefb";
        $cmd = "$makefb -poly $prefix.poly > $prefix.roots";
        my_system $cmd;
    }

    # Computing free relations
    if (! -f "$prefix.freerels") {
        banner("Free relations");
        my $freerel = "$param->{'cadodir'}/linalg/freerel";
        $cmd = "$freerel -poly $prefix.poly -fb $prefix.roots > $prefix.freerels";
        my_system $cmd;
    }

    # Sieving, removing duplicates, singleton removal
    # This is the same command, since we continue sieving until it works.
    my $nrels;
    $nrels = sieve($param);

    # Merge
    print "Merging...\n";
    $cmd = "$param->{'cadodir'}/linalg/merge".
      " -out  $prefix.merge.his" .
      " -mat $prefix.purged" .
      " -forbw $param->{'bwstrat'}" .
      " -prune $param->{'prune'}" .
      " -keep "  . $param->{'keep'} .
      " -maxlevel $param->{'maxlevel'}" .
      " -cwmax $param->{'cwmax'}" .
      " -rwmax $param->{'rwmax'}" .
      " -ratio $param->{'ratio'}";
      my_system "$cmd 2> $prefix.merge.stderr";
    my $bwcostmin=`tail $prefix.merge.his | grep "BWCOSTMIN:" | awk '{print \$NF}'`;
    chomp $bwcostmin;
    my $replay="$param->{'cadodir'}/linalg/replay";
    $cmd = $replay .
        " -his $prefix.merge.his" .
        " -index $prefix.index" .
        " -purged $prefix.purged" .
        " -out $prefix.small";
    if (defined($bwcostmin)) {
        print "Bwcostmin = $bwcostmin\n";
        $cmd = $cmd . " -costmin $bwcostmin";
    }
    my_system "$cmd 2> $prefix.replay.stderr";

    # Linear algebra
    print "Transposing...\n";
    $cmd = "$param->{'cadodir'}/linalg/transpose" .
      " -T $param->{'wdir'}" .
      " -skip $param->{'skip'}" .
      " -in $prefix.small -out $prefix.small.tr";
    my_system $cmd;
    if ($param->{'linalg'} eq 'bw') { 
        print "Calling Block-Wiedemann...\n";
        $cmd = "$param->{'cadodir'}/linalg/bw/bw.pl" .
        " seed=1" . # For debugging purposes, we use a deterministic BW
        " mt=$param->{'bwmt'}" .
        " matrix=$prefix.small.tr" .
        " mn=64" .
        " vectoring=64" .
        " multisols=1" .
        " wdir=$wdir/bw" .
        " tidy=0" .
        " solution=$prefix.W";
        my_system "/bin/sh -c '$cmd > $prefix.bw.stderr 2>&1'";
    } else {
        if ($param->{'linalg'} ne 'bl') {
            print "WARNING: I don't know linalg=$param->{'linalg'}" .
              ". Use bl as default.\n";
        }
        print "Calling Block-Lanczos...\n";
        $cmd = "$param->{'cadodir'}/linalg/bl/bl.pl" .
        " matrix=$prefix.small.tr" .
        " wdir=$wdir/bl" .
        " solution=$prefix.W";
        my_system "/bin/sh -c '$cmd > $prefix.bl.stderr 2>&1'";
    }
    print "Converting dependencies to CADO format...\n";
    $cmd = "$param->{'cadodir'}/linalg/bw/mkbitstrings " .
      " $prefix.W";
    my_system "$cmd > $prefix.ker_raw";
    my $nker = `wc -l < $prefix.ker_raw`;
    chomp $nker;
    print "We have computed $nker vectors of the Kernel.\n";

    # Characters
    print "Adding characters...\n";
    $cmd = "$param->{'cadodir'}/linalg/characters" .
      " -poly $prefix.poly" .
      " -purged $prefix.purged" .
      " -ker $prefix.ker_raw" .
      " -index $prefix.index" .
      " -rel $prefix.nodup" .
      " -small $prefix.small" .
      " -nker $nker" .
      " -skip $param->{'skip'}" .
      " -nchar $param->{'nchar'}";
    my_system "$cmd > $prefix.ker 2> $prefix.characters.stderr";

    my $ndepmax=`wc -l $prefix.ker | awk '{print \$1}'`;
    chomp $ndepmax;
    print "We have $ndepmax remaining after characters.\n";
    if ($ndepmax > $param->{'nkermax'}) {
        $ndepmax = $param->{'nkermax'};
    }

    # Sqrt
    print "Preparing $ndepmax squareroots...\n";
    $cmd = "$param->{'cadodir'}/linalg/allsqrt" .
      " $prefix.nodup $prefix.purged $prefix.index $prefix.ker" .
      " $prefix.poly" .
      " 0 $ndepmax ar $prefix.dep";
    my_system $cmd;
    
    my $i;
    for ($i = 0; $i < $ndepmax; $i++) {
        my $suf = sprintf("%03d", $i);
        print "Testing dependency number $i...\n";
        $cmd = "$param->{'cadodir'}/sqrt/naive/algsqrt" .
          "  $prefix.dep.alg.$suf $prefix.dep.rat.$suf $prefix.poly";
        my_system "$cmd> $prefix.fact.$suf 2>> $prefix.algsqrt.stderr";
        open FH, "< $prefix.fact.$suf";
        my $line = <FH>; chomp $line;
        if ($line ne 'Failed') {
            print $line . " ";
            $line = <FH>;
            print $line;
            close FH;
            last;
        }
        close FH;
    }
}
