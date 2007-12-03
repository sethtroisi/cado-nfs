#!/usr/bin/perl -w

# Wrapper script for block wiedemann.
#
# This is an evolution of the previous shell script. The latter grew a
# bit cumbersome in the end, so it seems to become more manageable in
# perl now.

use strict;
use warnings;

my $param = {
	# example data
	modulus=>"531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127",
	msize => 5000,
	dens => 10,
	# everything can be overriden from the command line with switches
	# of the form ``param=xxx''
	m=>4,
	n=>4,
	vectoring=>4,
	method=>'ifft',
	multisols=>0,
};

# Parameters that will be used only if nothing sets them earlier.
my $weak = {
	tidy=>1,
	wdir=>"$ENV{HOME}/Local/testmat",
};

sub dirname  { $_[0] =~ m{^(.*)/[^/]*$}; return $1 || '.'; };
sub basename { $_[0] =~ m{([^/]*)$}; return $1; };

my $srcdir = dirname $0;
my $bindir;

{
	my $cmd="make -s -C \"$srcdir\" variables";
	open F, "$cmd |";
	while (<F>) {
		/^BINARY_DIR=(.*)$/ && do { $bindir=$1; };
	}
	close F;
}

my $slavebindir=$bindir;

my @args = @ARGV;
while (defined(my $x = shift(@args))) {
	if ($x =~ m/^mn=(.*)$/) {
		# special case.
		$param->{'m'} = $1;
		$param->{'n'} = $1;
	} elsif ($x =~ m/^(.*)=(.*)$/) {
		$param->{$1} = $2;
	} elsif ($x =~ /\.cfg$/) {
		# read a possible config file, and apply options right
		# now.
		open my $xh, $x or die "$x: $!";
		unshift @args, <$xh>;
		close $xh;
	} else {
		# Then it should be a matrix file.
		$param->{'matrix'} = $x;
	}
}

sub action { print "@_\n"; system @_; }

# If we already have a matrix, parse its header to obtain some
# information.
if ($param->{'matrix'}) {
	my $f = $param->{'matrix'};
	open my $mh, $f or die "$f: $!";
	my @header = ();
	my $hline;
	HEADER: {
		do {
			last HEADER unless defined(my $x = <$mh>);
			$hline = $x if $x =~ /ROWS/;
			push @header, $x;
		} until (scalar @header == 10);
	}
	if (!defined($hline)) {
		die "No header line in first ten lines of $f";
	}
	$hline =~ /(\d+)\s+ROWS/i or die "bad header line in $f : $hline";
	$param->{'msize'} = $1;
	$hline =~ /(\d+)\s+COLUMNS/i or die "bad header line in $f : $hline";
	if ($param->{'msize'} != $1) {
		die "Matrix is not square ($param->{'msize'}x$1)";
	}
	$hline =~ /MODULUS\s+(\d+)/i or die "bad header line in $f : $hline";
	$param->{'modulus'} = $1;

	my $hash = `head -c 2048 $f | md5sum | cut -c1-8`;
	chomp($hash);
	$weak->{'wdir'} .= ".$hash";
}

for my $k (keys %$weak) {
	next if exists $param->{$k};
	$param->{$k}=$weak->{$k};
}

# ok -- now we're fed up with all of the $param stuff, so short-circuit
# all this, and import the names in the main namespace.

my $dumped = '';
sub dumpvar {
	my $x = shift @_;
	if (defined(my $y = $param->{$x})) {
		$dumped .= "$x=$y\n";
	}
}

my $matrix =	$param->{'matrix'};	dumpvar 'matrix';
my $modulus =	$param->{'modulus'};	dumpvar 'modulus';
my $msize =	$param->{'msize'};	dumpvar 'msize';
my $m =		$param->{'m'};		dumpvar 'm';
my $n =		$param->{'n'};		dumpvar 'n';
my $wdir =	$param->{'wdir'};	dumpvar 'wdir';
my $resume =	$param->{'resume'};	dumpvar 'resume';
my $tidy =	$param->{'tidy'};	dumpvar 'tidy';
my $vectoring =	$param->{'vectoring'};	dumpvar 'vectoring';
my $dens =	$param->{'dens'};	dumpvar 'dens';
my $dump =	$param->{'dump'};	dumpvar 'dump';
my $mt =	$param->{'mt'};		dumpvar 'mt';
my $machines =	$param->{'machines'};	dumpvar 'machines';
my $rsync =	$param->{'rsync'};	dumpvar 'rsync';
my $method =	$param->{'method'};	dumpvar 'method';
my $threshold =	$param->{'threshold'};	dumpvar 'threshold';
my $multisols =	$param->{'multisols'};	dumpvar 'multisols';

if ($param->{'dumpcfg'}) {
	print $dumped;
	exit 0;
}

if ($vectoring && ($n % $vectoring != 0)) {
	die "vectoring parameter must divide n";
}

# Prepare the directory.
system "rm -rf $wdir" unless $resume;
mkdir $wdir unless -d $wdir;

# In case of multiple machines, make sure that all machines see the
# directory.
my @mlist = split ' ', $machines || '';

my $shared=1;

if (@mlist) {
	my $tms = time;
	system "touch $wdir/$tms";
	for my $m (@mlist) {
		open my $xh, "ssh $m ls $wdir/$tms 2>&1 |";
		while (defined(my $x=<$xh>)) {
			print "($m)\t$x";
			chomp($x);
			if ($x ne "$wdir/$tms") {
				$shared=0;
			}
		}
		close $xh;
	}
	unlink "$wdir/$tms";
	if (!$shared) {
		print STDERR "Some machines can not see $wdir\n";
		if (!$rsync) {
			print STDERR "Either use rsync=1 or an nfs dir\n";
			exit 1;
		}
	}
	# TODO: what about inhomogeneous computations ?
        action "cp ${bindir}bw-slave $wdir/";
	$slavebindir="$wdir/";
}

if ($matrix) {
	action "cp $matrix $wdir/matrix.txt";
} else {
	# If no matrix parameter is set at this moment, then surely it
	# means we're playing with a random sample: create it !
	action "${bindir}bw-random $msize $modulus $dens > $wdir/matrix.txt";
}

action "${bindir}bw-balance --subdir $wdir"
	unless ($resume && -f "$wdir/matrix.txt.orig");
action "${bindir}bw-secure --subdir $wdir"
	unless ($resume && -f "$wdir/X0-vector");

if ($resume && -f "$wdir/X00") {
	# Make sure there's the proper number of vectors.
	my $dh;
	opendir $dh, $wdir;
	my $gx = scalar grep { /^X\d+$/ } readdir $dh;
	closedir $dh;
	if (scalar $gx != $m) { die "Bad m param ($m != $gx)\n"; }
	opendir $dh, $wdir;
	my $gy = scalar grep { /^Y\d+$/ } readdir $dh;
	closedir $dh;
	if (scalar $gy != $n) { die "Bad m param ($n != $gy)\n"; }
} else {
	action "${bindir}bw-prep --subdir $wdir $m $n"
}

if ($dump) {
	action "${bindir}bw-printmagma --subdir $wdir > $wdir/m.m";
}

sub rsync_push {
	if (!$shared && $rsync) {
		for my $m (@mlist) {
			action "rsync -av @_ $wdir/ $m:$wdir/";
		}
	}
}

sub rsync_pull {
	if (!$shared && $rsync) {
		for my $m (@mlist) {
			action "rsync -av @_ $m:$wdir/ $wdir/";
		}
	}
}

rsync_push "--delete";

sub compute_spanned {
	my @jlist = @_;
	for my $job (@jlist) {
		print "$job\n";
	}
	if (@mlist) {
		for my $i (0..$#jlist) {
			my $job = $jlist[$i];
			my $machine = $mlist[$i % scalar @mlist];
			next if fork;
			# kid code
			print "($machine/job$i)\t$job\n";
			open my $ch, "ssh $machine $job 2>&1 |";
			while (defined(my $x = <$ch>)) {
				print "($machine/job$i)\t$x";
			}
			close $ch;
			exit 0;
			# end of kid code.
		}
		for my $i (0..$#jlist) { wait; }
	} else {
		for my $job (@jlist) {
			system "$job\n";
		}
	}
}

my $m1=$m-1;
my $n1=$n-1;

# Gather the list of slave commands.
SLAVE : {
	my $exe="${slavebindir}bw-slave";
	if ($mt) {
		$exe .= "-mt --nthreads $mt";
	}
	$exe .= " --task slave --subdir $wdir";


	my @slavejobs = map
		{
			my ($i,$ni1)=($vectoring*$_,$vectoring*($_+1)-1);
			"$exe 0,$i $m1,$ni1";
		}
		(0..int($n/$vectoring) - 1);


	compute_spanned @slavejobs;
}

rsync_pull;

my @sols;

MASTER : {
	my $exe = "${bindir}bw-master";
	if ($method =~ /^(?:q(?:uadratic)?|old)$/) {
		$exe .= '-old';
	} else {
		die "threshold needed" unless $threshold;
		# TODO: allow runtime selection of method or (probably
		# smarter) runtime checking that the proper binary is
		# being used.
		$exe .= "2 -t $threshold";
	}
	$exe .= " --subdir $wdir matrix.txt $m $n";
	open my $mlog, ">$wdir/master.log";
	print "$exe\n";
	open my $mh, "$exe |";
	while (defined(my $x=<$mh>)) {
		print $mlog $x;
		print $x;
		if ($x =~ /LOOK \[\s*([\d\s]*)\]$/) {
			@sols = split ' ',$1;
		}
	}
	if (!scalar @sols) {
		die "No solution found";
	}
}

rsync_push;

# Aggressively select only one solution.
if (!$multisols) {
	@sols = ($sols[0]);
}

my @solfiles=();

MKSOL : {
	my $exe="${slavebindir}bw-slave";
	if ($mt) {
		$exe .= "-mt --nthreads $mt";
	}
	$exe .= " --task mksol --subdir $wdir";

	for my $s (@sols) {
		my @mksoljobs = map
			{
			my ($i,$ni1)=($vectoring*$_,$vectoring*($_+1)-1);
			"$exe --sc $s 0,$i $m1,$ni1";
			}
			(0..int($n/$vectoring) - 1);

		compute_spanned @mksoljobs;
		rsync_pull;
		action "${bindir}bw-gather --subdir $wdir $s --nbys $vectoring";

		if ($matrix) {
			my $d = dirname $matrix;
			system "cp \"$wdir/W0$s\" \"$d\"";
			push @solfiles, "$d/W0$s";
		}
	}

	if (@solfiles) {
		print "solution files\n";
		for my $f (@solfiles) { print "$f\n"; }
	}

	if ($matrix && !$multisols) {
		(my $sf = $matrix) =~ s/matrix/solution/;
		if ($sf ne $matrix) {
			my $s = $sols[0];
			system "cp \"$wdir/W0$s\" \"$sf\"";
			print "also: $sf\n";
		}
	}
}

if ($tidy) {
	system "rm -rf $wdir";
}
