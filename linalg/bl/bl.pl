#!/usr/bin/perl -w

# Wrapper script for block lanczos.

use strict;
use warnings;

use Time::HiRes qw(gettimeofday);

use Data::Dumper;

print "$0 @ARGV\n";

sub dirname  { $_[0] =~ m{^(.*)/[^/]*$}; return $1 || '.'; };
sub basename { $_[0] =~ m{([^/]*)$}; return $1; };

my $favorite_tmp = "$ENV{HOME}/Local";
$favorite_tmp = "/tmp" if ! -d $favorite_tmp;

# Parameters that will be used only if nothing sets them earlier.
my $weak = {
    tidy=>1,
    wdir=> "$favorite_tmp/bl-tmp",
};

my $srcdir = dirname $0;
my $bwbindir;
{
    my $cmd="make -s -C \"$srcdir\"/../bw variables";
    open F, "$cmd |";
    while (<F>) {
        /^BINARY_DIR=(.*)$/ && do { $bwbindir=$1; };
    }
    close F;

    if ($bwbindir eq '') {
        $bwbindir="$srcdir/../bw";
    }
}

my $param = {};
my @args = @ARGV;
while (defined(my $x = shift(@args))) {
    if ($x =~ m/^(.*)=(.*)$/) {
        $param->{$1} = $2;
    } else {
        die "bad argument $x";
    }
}

sub action {
    print scalar join(' ',@_), "\n";
    system @_;
}

# We could possibly do cp --symlink just as well, or maybe ln (cp --link
# is a gnu-ism). Not ln -s without caution, because this would require
# some knowledge about the right path back, which could be quite a bit of
# a hack. I'm happy with cp --link at the moment.
sub do_cp { my @x = @_; unshift @x, "cp", "--link"; action @x;
    if ($? >> 8 != 0) {
        shift @x;
        shift @x;
        unshift @x, "cp";
        action @x;
    }
}


# check parameters
my $dumped = '';
sub dumpvar {
    my $x = shift @_;
    if (defined(my $y = $param->{$x})) {
        $dumped .= "$x=$y\n";
    }
}

if (!defined($param->{'matrix'})) {
    die "Undefined matrix file\n";
}
if (!defined($param->{'solution'})) {
    die "Undefined solution file\n";
}


# The default working directory depends on the matrix.
{
    my $hash = `head -c 2048 $param->{'matrix'} | md5sum | cut -c1-8`;
    chomp($hash);
    $weak->{'wdir'} .= $hash;
}

for my $k (keys %$weak) {
    next if exists $param->{$k};
    $param->{$k}=$weak->{$k};
}

my $matrix =	$param->{'matrix'};	dumpvar 'matrix';
my $solution =	$param->{'solution'};	dumpvar 'solution';
my $machines =	$param->{'machines'};	dumpvar 'machines';
my $wdir =	$param->{'wdir'};	dumpvar 'wdir';
my $tidy =	$param->{'tidy'};	dumpvar 'tidy';

if ($param->{'dumpcfg'}) {
    print $dumped;
    exit 0;
}



# Prepare the directory.
print "wdir is $wdir\n";
system "rm -rf $wdir";
if (!-d $wdir) {
    mkdir $wdir or die
    "Cannot create $wdir: $! -- select something non-default with wdir=";
}

my @mlist = split ' ', $machines || '';
my $njobs = scalar @mlist;

open STDOUT, "| tee -a $wdir/solve.out";
open STDERR, "| tee -a $wdir/solve.err";

print "Solver started ", scalar localtime, "\nargs:\n$dumped\n";
print STDERR "Solver started ", scalar localtime, "\nargs:\n$dumped\n";
# my $version=`$srcdir/util/version.sh`;
# print "version info:\n$version";
# print STDERR "version info:\n$version";

my $tstart = gettimeofday;
my $tlast = $tstart;

my $times = { };

sub account {
    my $t = gettimeofday;
    my $k = $_[0];
    if (!exists($times->{$k})) {
        $times->{$k} = 0.0;
    }
    $times->{$_[0]} += $t - $tlast;
    $tlast = $t;
}


sub rsync_push {
    return unless $njobs;
    for my $m (@mlist) {
        action "rsync -av @_ $wdir/ $m:$wdir/";
    }
}

sub rsync_pull {
    return unless $njobs;
    for my $m (@mlist) {
        action "rsync -av @_ $m:$wdir/ $wdir/";
    }
}

do_cp $matrix, "$wdir/matrix.txt";

rsync_push "--delete";
my @cmdline=();
if (@mlist) {
    @cmdline=('orterun', '-np', $njobs, '--host', scalar join(',',@mlist));
}

if ($njobs) {
    action "${bwbindir}/bw-balance --matrix $matrix --nslices $njobs";
}

account 'io';

push @cmdline, "${srcdir}/MainLanczos", "$wdir/matrix.txt", "$wdir/ker";

action @cmdline;

account 'bl';

unlink $solution;
do_cp "$wdir/ker", $solution;


account 'io';

if ($tidy) {
    system "rm -rf $wdir";
}
account 'io';
print "-+" x 30, "\n";
print "times (in seconds):\n";
for my $k (keys %$times) {
    print "$k\t", sprintf("%.2f", $times->{$k}), "\n";
}
print "total\t", sprintf("%.2f",(gettimeofday-$tstart)), "\n";
print "-+" x 30, "\n";


