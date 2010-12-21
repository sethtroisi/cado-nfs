#!/usr/bin/perl -w
# usage: $0 $number_files_clean

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

my $nfiles=shift;
my $is_gzip = 0;
my $wdir = abs_path(dirname(dirname($0)));

sub last_line {
    my ($f) = @_;
    my $last = "";
    if ($f =~ /\.gz$/) {
            $last = `gzip -dc $f 2> /dev/null | tail -n 1`;
            chomp $last;
            return $last;
    }
    open FILE, "< $f" or die "Cannot open `$f' for reading: $!.\n";

    # That should be enough to catch the last line
    seek FILE, -512, 2;
    $last = $_ while <FILE>;
    close FILE;
    chomp $last;
    return $last;
}

sub clean_1file {
  my $f = shift;
  my $start;
  my $length;

  if (last_line($f) !~ /^# (Total \d+ reports)/) {
    if ($is_gzip) {
        # avoid gzip -d, since it fails on truncated files.
        basename($f) =~ /^(.*)\.rels\.(\d+)-(\d+)\.gz$/;
        my $f_new = "$wdir/working_rels/$1.rels.$2-$3";
        $start = $2;
        $length = $3 - $start;
        my $nlines = `gzip -dc $f 2> /dev/null | wc -l`;
        if ( $nlines == 0 ) {
            print basename($f)." corrupted, cleaned file...\n";
            unlink $f;
            unlink "$wdir/inprogress/$start,$length";
            system "touch $wdir/inqueue/$start,$length";
            return -1;
        }
        $nlines = $nlines - 1;
        system "gzip -dc $f 2> /dev/null | head -$nlines > $f_new";
        $f = $f_new;
    } else {
        my $nlines = `wc -l $f | cut -d" " -f1`;
        $nlines = $nlines - 1;
        rename $f, "$f.tmp";
        system "head -$nlines $f.tmp 2> /dev/null > $f";
        unlink "$f.tmp";
    }
    open FILE, "+< $f"
        or die "Cannot open `$f' for update: $!.\n";

    
    # Since the file is truncated, we assume that the last
    # reported special q was not completely sieved, so we remove it.
    # It would not be a good idea to try to save it (incompletely
    # sieved ranges are a nuisance).
    my @lastq;
    my $pos = 0;
    while (<FILE>) {
        # Keep track of the last two special q's
        if (/^### q=(\d+): roots?/) {
            shift @lastq if scalar @lastq == 2;
            push @lastq, { q => $1, pos => $pos };
        }
        $pos = tell FILE;
    }

    # Less than two special q's in this file: nothing to recover
    if (scalar @lastq < 2) {
        warn "File `$f' contains no useable data. Deleting...\n";
        close FILE;
        unlink $f;
        unlink "$f.gz" if ($is_gzip);
        basename($f) =~ /^.*\.rels\.(\d+)-(\d+)$/;
        $start = $1;
        $length = $2 - $start;
        unlink "$wdir/inprogress/$start,$length";
        system "touch $wdir/inqueue/$start,$length";
        return -2;
    }

    # Truncate the file and add a marker at the end
    truncate FILE, $lastq[-1]->{'pos'};
    seek FILE, $lastq[-1]->{'pos'}, 0;
    print FILE "# Warning: truncated file\n";
    close FILE;

    # Rename the file to account for the truncated range
    basename($f) =~ /^(.*)\.rels\.(\d+)-(\d+)$/;
    my $name = $1;
    $length = $3 - $2;
    my @r = ($2, $lastq[0]->{'q'}+1);
    my @s = ($r[1], $3);
    print "Truncating `".basename($f)."' to range $r[0]-$r[1]...\n";
    rename $f, "$wdir/output_rels/$name.rels.$r[0]-$r[1]";

    # Cleanup and moving
    unlink $f;
    unlink "$f.gz" if ($is_gzip);
    unlink "$wdir/inprogress/$r[0],$length";
    $f = "$wdir/output_rels/$name.rels.$r[0]-$r[1]";
    system "gzip $f";
    $length = $s[1] - $s[0];
    system "touch $wdir/inqueue/$s[0],$length";
    return 0;
  }

  basename($f) =~ /\.rels\.(\d+)-(\d+)(|\.gz)$/;
  $length = $2 - $1;
  rename $f, "$wdir/output_rels/".basename($f);
  $f = "$wdir/output_rels/".basename($f);
  unless ($is_gzip) {
      print basename($f)." completed, compress and move to output_rels.\n";
      system "gzip $f";
  }
  unlink "$wdir/inprogress/$1,$length";
  return 1;
}

# clean old files
opendir DIR, "$wdir/working_rels/"
    or die "Cannot open directory `$wdir/working_rels/': $!\n";
my @files;
if ($is_gzip) {
    @files = grep /\.rels\.\d+-\d+\.gz$/, readdir DIR;
} else {
    @files = grep /\.rels\.\d+-\d+$/, readdir DIR;
}
my %cache;
$cache{$_} = -M "$wdir/working_rels/$_" for @files;
@files = sort { $cache{$b} <=> $cache{$a} } @files;
closedir DIR;

if (! exists ($files[$nfiles]))
{
    $nfiles = scalar @files;
    print "The directory `$wdir/working_rels/' contains $nfiles relations files.\n";
}
my $last_file = `ls -go --time-style=long-iso "$wdir/working_rels/$files[$nfiles-1]"  | cut -d" "  -f4,5`;
chomp $last_file;

warn "I will clean up the files whose modified date is before : $last_file.\n";
warn "Are you OK to continue? [y/N] (30s timeout)\n";
my $r;
eval {
    local $SIG{'ALRM'} = sub { die "alarm\n" }; # NB: \n required
    alarm 30;
    $r = <STDIN>;
    alarm 0;
};
if ($@) {
    die unless $@ eq "alarm\n"; # propagate unexpected errors
}
die "Aborting...\n" unless $r =~ /^(y|yes|Y)/i;

my $ntrunc = 0;
my $nfiles0 = $nfiles;
for (@files) {
    last if ($nfiles-- == 0);
    $ntrunc++ if (clean_1file ("$wdir/working_rels/$_") >= 0);
}

my $p = int ( $ntrunc * 1000 / $nfiles0 ) / 10;
print "$p% ($ntrunc/$nfiles0) of files have been recovered.\n";

# move and gzip
opendir DIR, "$wdir/results_rels/"
  or die "Cannot open directory `$wdir/results_rels/': $!\n";
my @results_files;
if ($is_gzip) {
    @results_files = grep /\.rels\.\d+-\d+\.gz$/, readdir DIR;
    print "move relations files of results_rels to output_rels...\n";
    rename "$wdir/results_rels/$_", "$wdir/output_rels/$_" for (@results_files);
} else {
    @results_files = grep /\.rels\.\d+-\d+$/, readdir DIR;
    print "compress relations files with gzip and move to output_rels...\n";
    for (@results_files) {
        system "gzip $wdir/results_rels/$_";
        rename "$wdir/results_rels/$_.gz", "$wdir/output_rels/$_.gz";
    }
}
close DIR;


