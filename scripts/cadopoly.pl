#!/usr/bin/perl -w

# Usage: cadopoly
# 	read files params_* in dir "." and run polyselect
# 	the params_* must contain: different name, n, machines ...
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

# recover (to do)

# polysel
my $old_prefix;
my $file;
while (@files) {
	$file = shift @files;
	read_param(\%param, { strict => 1 }, "$file");
	read_machines();
	open FILE, "> $param{'prefix'}.n"
		or die "Cannot open `$param{'prefix'}.n' for writing: $!.\n";
	print FILE "n: $param{'n'}\n";
	close FILE;
	cmd ("mv $old_prefix.polysel_jobs $param{'prefix'}.polysel_jobs") if ($old_prefix);
	if (@files) {
		do_polysel_bench();
	} else {
		do_polysel();
	}
	$old_prefix = $param{'prefix'};
}

# info kjout
banner "Info kjout";
opendir DIR, "."
	or die "Cannot open directory `.': $!\n";
my @kjout_files = grep /\.kjout\.[\de.]+-[\de.]+$/,
   					readdir DIR;
closedir DIR;

my $time_max = 1000000;
info "Time max for one phase: $time_max"."s";

my @names = map { /^(\S+)\.kjout/; $1 } @kjout_files;
my %names;
foreach (@names) {
	$names{$_}++;	
}

my $time;
my $phase;

NAME : foreach my $name (sort keys %names) {
	my $Emin;
	my $best;
	my ($min_first, $min_second, $min_third);
	my ($max_first, $max_second, $max_third);

	foreach my $f (@kjout_files) {
		next unless $f =~ /^$name/;
		#shift @kjout_files;
		open FILE, "< $f"
			or die "Cannot open `$f' for reading: $!.\n";
		my $last;
		while (<FILE>) {
			if ( /#\s+(\S+)\s+phase took\s+([\d.]+)s/ ) {
				$time = $2; 
				$phase = $1;
				if ($time > $time_max ) {
					warn "$phase phase of polyselect for configuration $name is too long!\n";
					next NAME;
				}
	
				if ( $phase eq "First" ) {
        			$min_first = $time if (!defined $min_first || $time < $min_first);
        			$max_first = $time if (!defined $max_first || $time > $max_first);
				} elsif ( $phase eq "Second" ) {
        			$min_second = $time if (!defined $min_second || $time < $min_second);
        			$max_second = $time if (!defined $max_second || $time > $max_second);
				} else {
        			$min_third = $time if (!defined $min_third || $time < $min_third);
        			$max_third = $time if (!defined $max_third || $time > $max_third);
				}
			}

			$last = $_ if /E=/;
		}
		close FILE;
		
		next unless $last;
        $last =~ /E=([\d.]+)/;
        if (!defined $Emin || $1 < $Emin) {
            $Emin = $1;
            $best = $f;
		}
	}
		
    die "No polynomial was found for configuration $name!\n"
    	unless defined $Emin;

    info "The best polynomial for configuration $name is from `".basename($best)."' (E = $Emin).\n";
	$tab_level++;
	info "Time min for first phase: $min_first"."s,  second phase: $min_second".
		 "s,  third phase: $min_third"."s\n".
		 "Time max for first phase: $max_first"."s,  second phase: $max_second".
		 "s,  third phase: $max_third"."s\n";
	$tab_level--;
}

			
