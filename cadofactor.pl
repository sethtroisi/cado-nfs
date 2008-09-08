#!/usr/bin/perl -w

# General script for factoring integers with Cado-NFS.
#
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
#  - Use Kleinjung program instead of polyselect
#  - Bench polynomials with some sieve before selecting the best one.
#  - Enable a 'lowmem' option
#  - How to remove a computer from mach_desc? (should work, but test and
#    document).


use strict;
use warnings;

use Data::Dumper;

use File::Copy;
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

    # polyselect
    degree=>5,
    bmin=>1,
    bmax=>20,
    e=>1e8,
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
    keep=>160,
    keeppurge=>100000,
    maxlevel=>15,
    cwmax=>200,
    rwmax=>200,
    ratio=>1.5,
    bwstrat=>1,

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

sub splitpath {
    # Returns dirname and basename. 
    $_[0] =~ m{^(?:(.*)/)?(.*)$};
    my $dir=defined($1) ? $1 : ".";
    my $base=$2;
    if ($dir eq '.') {
        $dir = `pwd`;
        chomp $dir;
    }
    return ($dir, $base);
}

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
        ($param->{'cadodir'}) =  splitpath($0);
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
# Returns a pair ( $t, $ret )
#   where $t is a 0 or 1 depending whether the timout was reached (0
#   means timeout, 1 means success).
#   and $ret is the returned value (empty string if timeout).
#
# For timeout, we use the alarm mechanism. See perldoc -f alarm.
sub my_system_timeout {
    my ($cmd, $timeout) = @_;
    my $ret = "";
    eval {
        local $SIG{'ALRM'} = sub { die "alarm\n" }; # NB: \n required
        alarm $timeout;
        $ret = `$cmd`;
        alarm 0;
    };
    if ($@) {
        die unless $@ eq "alarm\n";   # propagate unexpected errors
        print "WARNING! The command $cmd timed out.\n";
        return (0, "");
    }
    else {
        return (1, $ret);
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
# Read a selection status file:
# Each line is of the form
#   hostname b
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
        if (/^\s*(\w*)\s+(\d+)/) {
            push @status, [$1, $2];
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
        print GRR "" . $tt[0] . " " . $tt[1] . "\n";
    }
    close GRR;
}

# A task is a pair (hostname b)
# This function ssh to hostname, and tests if the task is still alive,
# dead or finished.
# 1 means running
# 0 means finished
# -1 means (probably dead)
sub check_running_select_task {
    my $param = shift @_;
    my $name = $param->{'name'};
    my $mach_desc = shift @_;
    my ($host, $b) = @_;
    print "    $host \t$b:\n";
    my %host_data=%{$mach_desc->{$host}};
    my $wdir = $host_data{'tmpdir'};
    my $cadodir = $host_data{'cadodir'};

    # First check if the last line corresponds to finished job:
    my ($t, $lastline) = my_system_timeout(
        "ssh $host tail -1 $wdir/$name.poly.$b", 30);
    if (! $t) {
        print "      Assume task is still running...\n";
        return 1;
    }
    $_ = $lastline;
    if (/# generated/) {
        print "      finished!\n";
        return 0;
    }
    # If file is partial check its last modification time:
    my $date;
    ($t, $date) = my_system_timeout("ssh $host date +%s", 30);
    if (! $t) {
        print "      Assume task is still running...\n";
        return 1;
    }
    my $modifdate;
    ($t, $modifdate) = my_system_timeout(
        "ssh $host stat -c %Y $wdir/$name.poly.$b", 30);
    if (! $t) {
        print "      Assume task is still running...\n";
        return 1;
    }
    ## If didn't move for 10 minutes, assume it's dead
    if ($date > ($modifdate + 600)) {
        print "      dead ?!?\n";
        return -1;
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
    ($t, $ret) = my_system_timeout("ssh $mach mkdir -p $wdir", 60);
    if (! $t) { die "Connection to $mach timeout\n"; }
    ($t, $ret) = my_system_timeout("ssh $mach ls $wdir/$name.n", 120);
    if (! $t) { die "Connection to $mach timeout\n"; }
    chomp $ret;
    if ($ret eq "") {
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

    my %doneb;
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $wdir = $param->{'wdir'};
    opendir(DIR, $wdir) or die "can't opendir $wdir: $!";
    my @files = readdir(DIR);
    close DIR;
    foreach my $f (@files) {
        if ($f =~ /$name\.poly\.(\d+)/) {
            $doneb{$1}=1;
        }
    }
    foreach my $t (@$running) {
        $doneb{${$t}[1]}=1;
    }

    my @setb;
    for ($b = $param->{'bmin'}; $b <= $param->{'bmax'}; $b++) {
        if (! $doneb{$b}) {
            push @setb, $b;
        }
    }
    if (scalar @$running == 0 && scalar @setb == 0) {
        # finished!
        return 1;
    }

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
            if (! scalar @setb) { last; }
            $b = shift @setb;
            $t = [$m, $b ];
            my $mach = $m;
            my $wdir = $desc{'tmpdir'};
            my $bindir = $desc{'cadodir'};
            push_select_files($mach, $wdir, $param);
            my $cmd = "/bin/nice -$nice $bindir/polyselect/polyselect" .
              " -b $b" .
              " -e " . $param->{'e'} .
              " -degree " . $param->{'degree'} .
              " < $wdir/$name.n";
            my $ret = `ssh $mach "$cmd >& $wdir/$name.poly.$b&"`;
            print "    Starting $mach $b.\n";
            open FH, ">> " . $param->{'wdir'} . "/$name.cmd";
            print FH "ssh $mach \"$cmd >& $wdir/$name.poly.$b&\"\n";
            close FH;
            push @$running, $t;
        }
    }
    return 0;
}

sub get_logmualpha_value {
    my $filename = shift @_;
    open FILE, "$filename";
    my @lines = readline(FILE);
    close FILE;
    @lines = grep (/E=/, @lines);
    $_ = $lines[0];
    if (/E=(\d+\.\d*)/) {
        my $E = $1;
        return $E;
    } else {
        print STDERR "Can not parse output of $filename for geting logmu+alpha\n";
        return 0;
    }
}


# rsync a finished / dead task to the local working dir.
# return the logmu+alpha value.
sub import_select_task_result {
    my $param = shift @_;
    my $name = $param->{'name'};
    my $mach_desc = shift @_;
    my ($host, $b) = @_;
    my %host_data=%{$mach_desc->{$host}};
    my $rwdir = $host_data{'tmpdir'};
    my $lwdir = $param->{'wdir'};
    my $cmd = "rsync --timeout=30 $host:$rwdir/$name.poly.$b $lwdir/";
    my $ret = system($cmd);
    if ($ret != 0) {
        print STDERR "Problem when importing file $name.poly.$b from $host.\n";
        return 0;
    }

    return get_logmualpha_value("$lwdir/$name.poly.$b");
}


sub parallel_polyselect {
    my $param = shift @_;
    my $prefix = $param->{'wdir'} . "/" . $param->{'name'};
    my $b;
    my $t;
    my $bestb=0;
    my $effort=$param->{'e'};
    my $degree=$param->{'degree'};
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
            } else { # dead task!
                # do nothing, this will be restarted soon...
            }
        }
        print "  Start new tasks.\n";
        $finished = restart_select_tasks($param, \@newstatus, \%mach_desc);
        write_select_status_file($param, \@newstatus);
        if (!$finished) {
            my $delay = $param->{'delay'};
            print "Wait for $delay seconds before checking again.\n";
            sleep($delay);
        }
    }
    # Choose best according to logmu+alpha
    my $Emin;
    for ($b = $param->{'bmin'}; $b <= $param->{'bmax'}; $b++) { 
        if (! -f "$prefix.poly.$b") {
            die "Hey! Where is $prefix.poly.$b ????\n";
        }
        my $E = get_logmualpha_value("$prefix.poly.$b");
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
      "/linalg/duplicates -nrels $nrels $filestring > $prefix.nodup 2> $prefix.duplicates.stderr";
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
        if (/^\s*(\w*)\s+(\d+)\s+(\d+)/) {
            push @status, [$1, $2, $3];
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
        print GRR "" . $tt[0] . " " . $tt[1] . " " . $tt[2] . "\n";
    }
    close GRR;
}


# A task is a triple (hostname q0 q1)
# This function ssh to hostname, and tests if the task is still alive,
# dead or finished.
sub check_running_task {
    my $param = shift @_;
    my $name = $param->{'name'};
    my $mach_desc = shift @_;
    my ($host, $q0, $q1) = @_;
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
    my $date;
    ($t, $date) = my_system_timeout("ssh $host date +%s", 30);
    if (! $t) {
        print "      Assume task is still running...\n";
        return 1;
    }
    my $modifdate;
    ($t, $modifdate) = my_system_timeout(
        "ssh $host stat -c %Y $wdir/$name.rels.$q0-$q1", 30);
    if (! $t) {
        print "      Assume task is still running...\n";
        return 1;
    }
    ## If didn't move for 10 minutes, assume it's dead
    if ($date > ($modifdate + 600)) {
        print "      dead ?!?\n";
        return 0;
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
            if (! $h{'cores'}) { $h{'cores'}=1; }

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
    ($t, $ret) = my_system_timeout("ssh $mach mkdir -p $wdir", 60);
    if (! $t) { die "Connection to $mach timeout\n"; }
    ($t, $ret) = my_system_timeout("ssh $mach ls $wdir/$name.poly", 120);
    if (! $t) { die "Connection to $mach timeout\n"; }
    chomp $ret;
    if ($ret eq "") {
        print "    Pushing $name.poly to $mach.\n";
        system("rsync --timeout=120 $ldir/$name.poly $mach:$wdir/");
    }
    ($t, $ret) = my_system_timeout("ssh $mach ls $wdir/$name.roots", 120);
    if (! $t) { die "Connection to $mach timeout\n"; }
    chomp $ret;
    if ($ret eq "") {
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
            $t = [$m, $q_curr, $qend ];
            my $mach = $m;
            my $wdir = $desc{'tmpdir'};
            my $bindir = $desc{'cadodir'};
            push_files($mach, $wdir, $param);
            my $cmd = "/bin/nice -$nice $bindir/sieve/las" .
              " -I " . $param->{'I'} .
              " -poly $wdir/$name.poly" .
              " -fb $wdir/$name.roots" .
              " -q0 $q_curr -q1 $qend";
            my $ret = `ssh $mach "$cmd >& $wdir/$name.rels.$q_curr-$qend&"`;
            print "    Starting $mach $q_curr-$qend.\n";
            open FH, ">> " . $param->{'wdir'} . "/$name.cmd";
            print FH "ssh $mach \"$cmd >& $wdir/$name.rels.$q_curr-$qend&\"\n";
            close FH;
            push @$running, $t;
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
            # Task is finished. 
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
    foreach my $f (@files) {
        if ($f =~ /$name\.rels\.(\d+)-(\d+)/) {
            $nrels+= `grep -v "^#" $wdir/$f | wc -l`;
        }
    }
    return $nrels;
}
 

sub parallel_sieve {
    my $param = shift @_;
    my $finished = 0;
    my $nrels = count_rels($param);
    print "We have found $nrels relations in working dir\n";
    print "Let's start the main loop!\n";
    my $prev_check = 0;
    my @files=();
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

    # Look for finished files.
    my @done=();

    {
        my $filename = "$prefix.rels";
        if (-f $filename) {
            print "Found file $filename...";
            my $lines = `grep -v "^#" $filename | wc -l`;
            chomp $lines;
            print "$lines rels\n";
            push @done, $filename;
            $nrels += $lines;
        }
    }

    while (1) {
        my $qend = $qcurr+$param->{'qrange'};
        my $filename = "$prefix.rels.$qcurr-$qend";
        last unless -f $filename;
        $qcurr = $qend;
        print "Found file $filename...";
        my $lines = `grep -v "^#" $filename | wc -l`;
        chomp $lines;
        print "$lines rels\n";
        push @done, $filename;
        $nrels += $lines;
    }

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
            my_system "$cmd >& $filename";
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
    print "Bwcostmin = $bwcostmin\n";
    my $replay="$param->{'cadodir'}/linalg/replay";
    $cmd = $replay .
        " -his $prefix.merge.his" .
        " -index $prefix.index" .
        " -costmin $bwcostmin" .
        " -purged $prefix.purged" .
        " -out $prefix.small";
    my_system "$cmd 2> $prefix.replay.stderr";

    # Linear algebra
    print "Transposing...\n";
    $cmd = "$param->{'cadodir'}/linalg/transpose" .
      " -T $param->{'wdir'}" .
      " -in $prefix.small -out $prefix.small.tr";
    my_system $cmd;
    if ($param->{'linalg'} eq 'bw') { 
        print "Calling Block-Wiedemann...\n";
        $cmd = "$param->{'cadodir'}/linalg/bw/bw.pl" .
        " mt=$param->{'bwmt'}" .
        " matrix=$prefix.small.tr" .
        " mn=64" .
        " vectoring=64" .
        " multisols=1" .
        " wdir=$wdir/bw" .
        " solution=$prefix.W";
        my_system "$cmd >& $prefix.bw.stderr";
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
        my_system "$cmd >& $prefix.bl.stderr";
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
        print "Testing dependency numnber $i...\n";
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
