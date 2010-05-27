#!/usr/bin/perl -w

# Usage: cadosieve.pl
# 	read files params_* in dir "." and estimate time for sieve
# 	the params_* must contain: different name, machines, params sieve
# 							   and filter ...
#							   qrange=1000, 2000, 4000 ...
# 	( read one file *.poly in dir and generate the others and run makefb )
#
# TODO: generate input
# 		...

use Cwd qw(abs_path);
use File::Basename;
use lib abs_path(dirname(dirname($0)));
use cadofct;
use strict;
use warnings;

# input
opendir DIR, "."
	or die "Cannot open directory `.': $!\n";
my @files = grep /^params_/,
				readdir DIR;
closedir DIR;
@files = sort @files;

# recover
my @done = <*.sieve_done>;
my @recover = <*.sieve_jobs>;

my @poly = <*.poly>;
die unless (@poly);

my %link_name_params;

# lpbr => [jump to next interval, nrels max with 20% duplicates]
my %max_rels = ( 25 => [  500000,   5000000],
				 26 => [ 1000000,   8000000],
				 27 => [ 2000000,  14000000], 
				 28 => [ 5000000,  25000000],
				 29 => [10000000,  50000000],
				 30 => [20000000, 100000000],
				 31 => [40000000, 200000000] );
# sieve
my $file;
while (@files) {
	$file = shift @files;
	read_param(\%param, { strict => 1 }, "$file");
	$link_name_params{"$param{'name'}"} = $file;
	next if (@done);
	if (@recover) {
		next if ($recover[0] ne basename("$param{'prefix'}.sieve_jobs"));
		shift @recover;
	}
	if ($poly[0] ne basename("$param{'prefix'}.poly")) {
		cmd ("sed '/^rlim/d;/^alim/d;/^lpbr/d;/^lpba/d;/^mfbr/d;/^mfba/d;
				/^rlambda/d;/^alambda/d' $poly[0] > $param{'prefix'}.poly");
	    open FILE, ">> $param{'prefix'}.poly"
    	    or die "Cannot open `$param{'prefix'}.poly' for writing: $!.\n";
    	print FILE "$_: $param{$_}\n"
        	for qw(rlim alim lpbr lpba mfbr mfba rlambda alambda);
    	close FILE;
	}
	do_factbase() unless (-f "$param{'prefix'}.roots");
	
	read_machines();
	my @r;
	my @t = @{$max_rels{"$param{'lpbr'}"}};
	$r[0] = $param{'qmin'} + $param{'qrange'};
	$r[1] = $param{'qmin'} + $t[0];
	my @ranges;
	while ( $r[1] < 50 * $t[1] ) {
		cmd ("touch $param{'prefix'}.rels.$r[0]-$r[1]");
		$r[0] += $t[0];
		$r[1] += $t[0];
	}
	if (@files) {
		do_sieve_bench(\@t);
	} else {
		do_sieve_bench(\@t, 1);
		cmd ("touch $param{'prefix'}.sieve_done");
	}
}

# info rels
banner "Info rels";
my @rels_files;

foreach my $name (sort keys %link_name_params) {
    my $nrels = 0;
	my $reports;
	my $time_by_rels;
	my $estimate_time = 0;

	opendir DIR, "."
		or die "Cannot open directory `.': $!\n";
	@rels_files = grep /^$name\.rels\.[\de.]+-[\de.]+\.gz$/,
   						readdir DIR;
	closedir DIR;
	my @int_files = map { /^$name\.rels\.([\de.]+)-/; $1 } @rels_files;
	@int_files = sort ( {$a <=> $b } @int_files );

	$file = $link_name_params{"$name"};
	read_param(\%param, { strict => 1 }, "$file");
	sub plus_qrange {
		my $x = shift;
		return ( $x + $param{'qrange'} );
	}
	my @rels_files = map { /^(\d+)$/; "$name.rels.$1-".plus_qrange($1).".gz" } @int_files;

	foreach my $f (@rels_files) {
		if ( last_line($f) =~ /# Total (\d+) reports \[(\S+)s/ ) {
				info "$f: ".last_line($f)."\n";
				$reports = $1;
				$time_by_rels = $2; 
				my @t = @{$max_rels{"$param{'lpbr'}"}};
				my $rels = $t[0] * $reports / $param{'qrange'};
        		$nrels += $rels;
				if ($nrels < $t[1]) {
					$estimate_time += $rels * $time_by_rels;
				} else {
					$estimate_time += ( $t[1] - ($nrels - $rels) ) * $time_by_rels;
					last;
				}
		}
	}
	
	$tab_level++;
	info "Estimate time for configuration $name : \033[01;31m$estimate_time"."s, ".
								format_dhms($estimate_time)."\033[01;00m\n";
	$tab_level--;
}
	
