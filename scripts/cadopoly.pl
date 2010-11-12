#!/usr/bin/perl -w

# Usage: cadopoly.pl $paramfile1 $paramfile2 ...
#        read paramfiles and run polyselect
#        the paramfiles must contain: different name!, n, machines, params polyselect ...
#
# TODO: generate input

use Cwd qw(abs_path);
use File::Basename;
use lib abs_path(dirname(dirname($0)));
use cadofct;
use POSIX qw(ceil);
use strict;
use warnings;

# input
print "$0 @ARGV\n";
my @files;
while (@ARGV) {
  push @files, shift @ARGV;
}
@files = sort @files;

# recover
my @done = <*.polysel_done>;
my @recover = <*.polysel_jobs>;

# polyselect
my %link_name_params;
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
  rename ( "$old_prefix.polysel_jobs", "$param{'prefix'}.polysel_jobs" )
    if ($old_prefix);
  read_machines();
  open FILE, "> $param{'prefix'}.n"
    or die "Cannot open `$param{'prefix'}.n' for writing: $!.\n";
  print FILE "n: $param{'n'}\n";
  close FILE;
  if (@files) {
    do_polysel_bench();
  } else {
    do_polysel_bench(1);
    unlink "$param{'prefix'}.polysel_jobs";
    open FILE, "> $param{'prefix'}.polysel_done"
      or die "Cannot open `$param{'prefix'}.polysel_done' for writing: $!.\n";
    close FILE;
  }
  $old_prefix = $param{'prefix'};
}

# info kjout
banner "Info kjout";
my @kjout_files;
my $time;
my $phase;
my %h = ( First  => 0,
          Second => 1,
          Third  => 2 );
foreach my $name (sort keys %link_name_params) {
  my (@Emax, $Emin, $Emoy);
  my $best;
  my (@min, @max, @total);
  my $total;

  opendir DIR, "."
    or die "Cannot open directory `.': $!\n";
  @kjout_files = grep /^$name\.kjout\.[\de.]+-[\de.]+$/, readdir DIR;
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
        $total[$h{$phase}] += $time;
        $min[$h{$phase}] = $time if (!defined $min[$h{$phase}] || $time < $min[$h{$phase}]);
        $max[$h{$phase}] = $time if (!defined $max[$h{$phase}] || $time > $max[$h{$phase}]);
      }
      $murphy = $_ if /Murphy/;
    }
    close FILE;
		
    next unless $murphy;
    $murphy =~ /\)=(.+)$/;
    $Emoy += $1;
    if (!defined $Emax[2] || $1 > $Emax[2]) {
      if (!defined $Emax[1] || $1 > $Emax[1]) {
        $Emax[2] = $Emax[1];
        if (!defined $Emax[0] || $1 > $Emax[0]) {
          $Emax[1] = $Emax[0];
          $Emax[0] = $1;
          $best = $f;
        } else {
          $Emax[1] = $1;
        }
      } else {
        $Emax[2] = $1;
      }
    }
    $Emin = $1 if (!defined $Emin || $1 < $Emin);
  }
		
  die "No polynomial was found for configuration $name!\n"
    unless defined $Emax[0];

  $Emoy = $Emoy / $size;
  my $ampl = $Emax[0]-$Emin;
  $file = $link_name_params{"$name"};
  read_param(\%param, { strict => 1 }, "$file");
  my $size_expect = ceil (($param{'kjadmax'} -  $param{'kjadmin'} ) /
                          $param{'kjadrange'});
  die "kjout missing for configuration $name! ($size on $size_expect)\n"
    if ( $size != $size_expect );

  # print infos
  info "Params: M $param{'kjM'} l $param{'kjl'} kmax $param{'kjkmax'} ".
       "adrange $param{'kjadrange'} admin $param{'kjadmin'} admax $param{'kjadmax'}\n".
       "        degree $param{'degree'} p0max $param{'kjp0max'} ".
       "keep $param{'kjkeep'} incr $param{'kjincr'} pb $param{'kjpb'}\n";
  $tab_level++;
  info "The best polynomial for $name is from `".basename($best)."'\n".
       "Emax = $Emax[0] - Emax2 = $Emax[1] - Emax3 = $Emax[2]\n".
       "Emoy = \033[01;31m$Emoy\033[01;00m - ".
       "Emin = $Emin - ". "Ampl = $ampl\n".
       "Time min for first phase: $min[0]"."s,  second phase: $min[1]".
       "s,  third phase: $min[2]"."s\n".
       "Time max for first phase: $max[0]"."s,  second phase: $max[1]".
       "s,  third phase: $max[2]"."s\n".
       "Time total for first phase: $total[0]"."s,  second phase: ".
       $total[1]."s, third phase: $total[2]"."s\n".
       "Time total: $total"."s,  \033[01;31m".format_dhms($total)."\033[01;00m\n";
  print "\n";
  $tab_level--;
}
			
