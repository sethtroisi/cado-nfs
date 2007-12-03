#!/usr/bin/perl -w

# Wrapper script for block wiedemann.
#
# This is an evolution of the previous shell script. The latter grew a
# bit cumbersome in the end, so it seems to become more manageable in
# perl now.

use strict;
use warnings;

use Time::HiRes qw(gettimeofday); # For seed values

use Data::Dumper;

print "$0 @ARGV\n";

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
	maxload=>1, # How many simultaneous jobs on one machine.
};


my $favorite_tmp = "$ENV{HOME}/Local";
$favorite_tmp = "/tmp" if ! -d $favorite_tmp;

# Parameters that will be used only if nothing sets them earlier.
my $weak = {
	tidy=>1,
	wdir=> "$favorite_tmp/testmat",
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

	if ($bindir eq '') {
		$bindir="$srcdir/";
	}
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

sub parse_matrix_header {
	my $f = shift @_;
	open my $mh, $f or die "$f: $!";
	my @header = ();
	my $h = {};
	my $hline;
	HEADER: {
		do {
			last HEADER unless defined(my $x = <$mh>);
			$hline = $x if $x =~ /ROWS/;
			if (!defined($hline)) {
				$hline = $x if $x =~ /^\d+\s\d+$/;
			}
			push @header, $x;
		} until (scalar @header == 10);
	}
	if (!defined($hline)) {
		die "No header line in first ten lines of $f";
	}
	if ($hline =~ /^(\d+)\s(\d+)$/) {
		$h->{'msize'} = $1;
		$h->{'modulus'} = 2;
		print "CADO format $1 rows $2 columns\n";
	} else {
		$hline =~ /(\d+)\s+ROWS/i
			or die "bad header line in $f : $hline";
		$h->{'msize'} = $1;
		$hline =~ /(\d+)\s+COLUMNS/i
			or die "bad header line in $f : $hline";
		if ($h->{'msize'} != $1) {
			die "Matrix is not square ($h->{'msize'}x$1)";
		}
		$hline =~ /MODULUS\s+(\d+)/i
			or die "bad header line in $f : $hline";
		$h->{'modulus'} = $1;
	}
	my $hash = `head -c 2048 $f | md5sum | cut -c1-8`;
	chomp($hash);
	$h->{'hash'} = $hash;
	return $h;
}

if ($param->{'matrix'}) {
	my $h = parse_matrix_header $param->{'matrix'};
	$param->{'msize'} = $h->{'msize'};
	$param->{'modulus'} = $h->{'modulus'};
	$weak->{'wdir'} .= ".$h->{'hash'}";
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
my $maxload =	$param->{'maxload'};	dumpvar 'maxload';
my $seed =	$param->{'seed'};	dumpvar 'seed';
my $precond =	$param->{'precond'};	dumpvar 'precond';

if ($param->{'dumpcfg'}) {
	print $dumped;
	exit 0;
}

my $seeding = '';
if ($seed) {
	$seeding = "--seed $seed";
} else {
	$seed = (gettimeofday)[1] % 999983;
	print "Selecting seed $seed\n";
	$seeding = "--seed $seed";
}



if ($vectoring && ($n % $vectoring != 0)) {
	die "vectoring parameter must divide n";
}

sub magmadump {
	if ($dump) {
		action "${bindir}bw-printmagma --subdir $wdir > $wdir/m.m";
	}
}

# Prepare the directory.
system "rm -rf $wdir" unless $resume;
if (!-d $wdir) {
	mkdir $wdir or die
		"Cannot create $wdir: $! -- select something non-default with wdir=";
}

if ($resume) {
	if ($precond && -f "$wdir/precond.txt") {
		print "There is already a $wdir/precond.txt file; let's hope it is the same as $precond !\n"
	} elsif ($precond && ! -f "$wdir/precond.txt")  {
		die "$wdir/precond.txt is missing";
	}
	# Having a precond file present but not required by cmdline is
	# conceivable if I introduce rhs someday.
	#
	# Having none is ok of course.
} else {
	if ($precond) {
		action "cp $precond $wdir/precond.txt";
	}
}

# In case of multiple machines, make sure that all machines see the
# directory.
my @mlist = split ' ', $machines || '';

my $shared=1;

if (@mlist) {
	my $tms = time;
	system "touch $wdir/$tms";
	for my $m (@mlist) {
		open my $xh, "ssh -n $m ls $wdir/$tms 2>&1 |";
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
        action "cp ${bindir}bw-slave-mt $wdir/";
	$slavebindir="$wdir/";
}

if ($resume) {
	die "Cannot resume: no matrix file" unless -f "$wdir/matrix.txt";
	my $h = parse_matrix_header "$wdir/matrix.txt";
	$param->{'msize'} = $h->{'msize'};
	$param->{'modulus'} = $h->{'modulus'};
} else {
	if ($matrix) {
		action "cp $matrix $wdir/matrix.txt";
	} else {
		# If no matrix parameter is set at this moment, then surely it
		# means we're playing with a random sample: create it !
		action "${bindir}bw-random $seeding $msize $modulus $dens > $wdir/matrix.txt";
	}
}

action "${bindir}bw-balance --subdir $wdir"
	unless ($resume && -f "$wdir/matrix.txt.old");
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
	action "${bindir}bw-prep $seeding --subdir $wdir $m $n"
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
	# for my $job (@jlist) { print "$job\n"; }
	if (@mlist) {
		my @kids=();
		my @assignments=();
		for my $j (1..($maxload * scalar @mlist)) {
			push @assignments, [];
		}
		for my $j (0..$#jlist) {
			my $job = $jlist[$j];
			push @{$assignments[$j % scalar @assignments]}, $j;
		}

		for my $i (0..$#assignments) {
			my $a = $assignments[$i];
			my $machine = $mlist[$i % scalar @mlist];
			my $pid=fork;
			if ($pid) {
				push @kids, $pid;
				next;
			}
			for my $j (@$a) {
				my $job = $jlist[$j];
				# kid code
				print "($machine/job$j)\t$job\n";
				open my $ch, "ssh -n $machine $job 2>&1 |";
				while (defined(my $x = <$ch>)) {
					print "($machine/job$j)\t$x";
				}
				close $ch;
			}
			# end of kid code.
			exit 0;
		}
		print "Child processes: ", join(' ', @kids), "\n";
		my @dead=();
		for my $i (0..$#assignments) { push @dead, wait; }
		print "Reaped processes: ", join(' ', @dead), "\n";
	} else {
		for my $job (@jlist) {
			action $job;
		}
	}
}

sub nlines {
	my $f = shift @_;
	open my $fh, "<$f" or return 0;
	my $l=0;
	while (<$fh>) { $l++; }
	close $fh;
	return $l;
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


	my @slavejobs = ();
	
	for my $vi (0..int($n/$vectoring) - 1) {
		my ($i,$ni1)=($vectoring*$vi,$vectoring*($vi+1)-1);
		# That's really a crude approximation.
		my $l_approx=$msize/$m + $msize/$n;
		my @linecounts=();
		for my $j (0..$m1) {
			my $x = nlines "$wdir/A-0$m1-00";
			if ($x < $l_approx) {
				# Then this one is not finished. We have
				# to go forward.
				last;
			}
			push @linecounts, $x;
		}
		if (scalar @linecounts == $m) {
			for my $j (0..$m1) {
				print "$wdir/A-0$m1-00 : $linecounts[$j] lines, over.\n";
			}
		} else {
			push @slavejobs, "$exe 0,$i $m1,$ni1";
		}
	}

	if (scalar @slavejobs) {
		print "--- slave jobs: ---\n";
		for my $j (@slavejobs) {
			print "// $j\n";
		}
		compute_spanned @slavejobs;
	}
}

rsync_pull;

my @sols;

MASTER : {
	RECYCLE_MASTER: {
		if ($resume) {
			my $allfiles=-f "$wdir/master.log";
			for my $x (0..$m+$n-1) {
				my $p = sprintf "%02d", $x;
				$allfiles = $allfiles && -f "$wdir-$p";
			}
			if ($allfiles) {
				print "master: all files found, reusing\n";
			} else {
				last RECYCLE_MASTER;
			}
			open my $ph, "tail -1 $wdir/master.log |";
			while (defined(my $x=<$ph>)) {
				if ($x =~ /LOOK \[\s*([\d\s]*)\]$/) {
					print "$x\n";
					@sols = split ' ',$1;
				}
			}
			if (!scalar @sols) {
				magmadump;
				die "No solution found";
			}
			last MASTER;
		}
	}

	my $exe = "${bindir}bw-master";
	if ($modulus eq '2') {
		$exe .= '-binary';
	} elsif ($method =~ /^(?:q(?:uadratic)?|old)$/) {
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
		magmadump;
		die "No solution found ; seed was $seed";
	}
}

rsync_push;

# Aggressively select only one solution.
if (!$multisols) {
	@sols = ($sols[0]);
}

my @sols_found  = ();

# within the matrix source directory
my @solfiles=();

MKSOL : {
	my $exe="${slavebindir}bw-slave";
	if ($mt) {
		$exe .= "-mt --nthreads $mt";
	}
	$exe .= " --task mksol --subdir $wdir";

	my @mksoljobs=();
	for my $s (@sols) {
		if ($resume && -f "$wdir/W0$s") {
			print "$wdir/W0$s already exists\n";
			next;
		}
		my $t = [ map
			{
			my ($i,$ni1)=($vectoring*$_,$vectoring*($_+1)-1);
			"$exe --sc $s 0,$i $m1,$ni1";
			}
			(0..int($n/$vectoring) - 1) ];

		push @mksoljobs, @$t;
	}

	if (scalar @mksoljobs) {
		print "--- mksol jobs: ---\n";
		for my $t (@mksoljobs) {
			#for my $j (@$t) {
			#print "// $j\n";
			#}
			print "// $t\n";
		}
		compute_spanned @mksoljobs;
			# for my $t (@mksoljobs) {
			#	compute_spanned @$t;
			#	rsync_pull;
			# }
	}

	if (!$multisols) {
		@sols = ($sols[0]);
		print "Restricting to solution @sols\n";
	}

	for my $s (@sols) {
		action "${bindir}bw-gather --subdir $wdir $s --nbys $vectoring";
		if (-f "$wdir/W0$s") {
			push @sols_found, "W0$s";
		}
		if ($matrix) {
			my $d = dirname $matrix;
			if (-f "$wdir/W0$s") {
				system "cp \"$wdir/W0$s\" \"$d\"";
				push @solfiles, "$d/W0$s";
			}
		}
	}

	if ($multisols) {
		print scalar @sols_found, " solutions found ",
			"(from ", scalar @sols, " out of master)\n";
		for my $f (@sols_found) { print "$f\n"; }
	}

	if (@solfiles) {
		print "solution files\n";
		for my $f (@solfiles) { print "$f\n"; }
	}

	if ($matrix) {
		if (scalar @solfiles != scalar @sols) {
			print "POSSIBLE BUG: not all wanted solutions were found\n";
			$tidy=0;
		}

		if (!$multisols && @solfiles) {
			(my $sf = $matrix) =~ s/matrix/solution/;
			if ($sf ne $matrix) {
				my $s = $solfiles[0];
				system "cp --link \"$s\" \"$sf\"";
				print "also: $sf\n";
			}
		}
	}
}

if ($dump) {
	magmadump;
}

if ($tidy) {
	system "rm -rf $wdir";
}

print "Seed value for this run was $seed\n";
