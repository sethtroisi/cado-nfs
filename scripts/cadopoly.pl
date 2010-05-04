#!/usr/bin/perl -w

# Usage: cadopoly.pl
# 	read files params_* in dir "." and run polyselect
# 	the params_* must contain: different name, n, machines, params polyselect ...
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
my @done = <*.polysel_done>;
my @recover = <*.polysel_jobs>;

my %link_name_params;

# polysel
my $old_prefix;
my $file;
while (@files) {
	$file = shift @files;
	read_param(\%param, { strict => 1 }, "$file");
	$link_name_params{"$param{'name'}"} = $file;
	next if (@done);
	if (@recover) {
		next if ($recover[0] ne basename("$param{'prefix'}.polysel_jobs"));
		shift @recover;
	}
	cmd ("mv $old_prefix.polysel_jobs $param{'prefix'}.polysel_jobs") if ($old_prefix);
	read_machines();
	open FILE, "> $param{'prefix'}.n"
		or die "Cannot open `$param{'prefix'}.n' for writing: $!.\n";
	print FILE "n: $param{'n'}\n";
	close FILE;
	if (@files) {
		do_polysel_bench();
	} else {
		do_polysel_bench(1);
		cmd ("mv $param{'prefix'}.polysel_jobs $param{'prefix'}.polysel_done");
	}
	$old_prefix = $param{'prefix'};
}

# info kjout
banner "Info kjout";
my @kjout_files;

my $time_max = 100000;
info "Time max for one phase: $time_max"."s";

my $time;
my $phase;

NAME : foreach my $name (sort keys %link_name_params) {
	my $Emoy;
	my $Emax;
	my $Emax2;
	my $Emax3;
	my $Emin;
	my $best;
	my ($min_first, $min_second, $min_third);
	my ($max_first, $max_second, $max_third);
	my ($total_first, $total_second, $total_third);
	my $total;

	opendir DIR, "."
		or die "Cannot open directory `.': $!\n";
	@kjout_files = grep /^$name\.kjout\.[\de.]+-[\de.]+$/,
   						readdir DIR;
	closedir DIR;
	my $size = scalar(@kjout_files);

	foreach my $f (@kjout_files) {
		open FILE, "< $f"
			or die "Cannot open `$f' for reading: $!.\n";
		my $murphy;
		while (<FILE>) {
			if ( /#\s+(\S+)\s+phase took\s+([\d.]+)s/ ) {
				$time = $2; 
				$total += $time;
				$phase = $1;
				if ($time > $time_max ) {
					warn "$phase phase of polyselect for configuration $name is too long!\n";
					next NAME;
				}
	
				if ( $phase eq "First" ) {
					$total_first += $time;
        			$min_first = $time if (!defined $min_first || $time < $min_first);
        			$max_first = $time if (!defined $max_first || $time > $max_first);
				} elsif ( $phase eq "Second" ) {
					$total_second += $time;
        			$min_second = $time if (!defined $min_second || $time < $min_second);
        			$max_second = $time if (!defined $max_second || $time > $max_second);
				} else {
					$total_third += $time;
        			$min_third = $time if (!defined $min_third || $time < $min_third);
        			$max_third = $time if (!defined $max_third || $time > $max_third);
				}
			}

			$murphy = $_ if /Murphy/;
		}
		close FILE;
		
		next unless $murphy;
        $murphy =~ /\)=(.+)$/;
		$Emoy += $1;
        if (!defined $Emax3 || $1 > $Emax3) {
        	if (!defined $Emax2 || $1 > $Emax2) {
        		if (!defined $Emax || $1 > $Emax) {
					$Emax3 = $Emax2;
					$Emax2 = $Emax;
					$Emax = $1;
            		$best = $f;
				} else {
					$Emax3 = $Emax2;
					$Emax2 = $1;
				}
			} else {
				$Emax3 = $1;
			}
		}
        if (!defined $Emin || $1 < $Emin) {
            $Emin = $1;
		}
	}
		
    die "No polynomial was found for configuration $name!\n"
    	unless defined $Emax;

	$Emoy = ($Emoy-$Emin-$Emax) / ($size-2);
	$file = $link_name_params{"$name"};
	read_param(\%param, { strict => 1 }, "$file");
	info "Params: M $param{'kjM'} l $param{'kjl'} kmax $param{'kjkmax'} ".
		 	"adrange $param{'kjadrange'} admin $param{'kjadmin'} admax $param{'kjadmax'}\n".
		 	"        degree $param{'degree'} p0max $param{'kjp0max'} ".
			"keep $param{'kjkeep'} incr $param{'kjincr'} pb $param{'kjpb'}\n";
	$tab_level++;
    info "The best polynomial for $name is from `".basename($best)."'\n";
	info "Emax = $Emax - ".
		 "Emax2 = $Emax2 - ".
		 "Emax3 = $Emax3 - ".
		 "Emoy = \033[01;31m$Emoy\033[01;00m  ";
	info "Time min for first phase: $min_first"."s,  second phase: $min_second".
		 "s,  third phase: $min_third"."s\n".
		 "Time max for first phase: $max_first"."s,  second phase: $max_second".
		 "s,  third phase: $max_third"."s\n";
	info "Time total for first phase: $total_first"."s,  second phase: ".
$total_second."s, third phase: $total_third"."s\n";
	info "Time total: $total"."s,  \033[01;31m".format_dhms($total)."s\033[01;00m\n";
	print "\n";
	$tab_level--;
}

			
