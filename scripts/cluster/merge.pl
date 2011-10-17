#!/usr/bin/env perl

# This script has been taken from the rsa768 scripts, and mildly adapted.
# Here is the readme from the cacao rsa768 fileset:
#       This is another mergeing tool. It does the combination of
#       scan_holes, merge_results, and clean_file, while together being
#       hardened against possible bugs. It recovers the stats line from
#       the files, or builds them by counting the relations found. Stats
#       lines are saved in separate files. Final files are free of
#       [cacao] lines.
#
# The changes so far are as follows:
#  - the relation format has been changed from franke-kleinjung to cado.
#  - the nextprime binary is no longer used. Instead, we rely on the
#    `roots' binary to provide all the roots for a given range of
#    special-qs with the -q0 -q1 args. Note that the format of `roots'
#    differs depending on whether the program is used for a single prime
#    or for a prime range. This reduces dramatically the number of calls
#    to external programs, which is better when the relation yield per
#    special q is small.
#  - [TODO] since las does not sort roots before sieving, we have to accomodate
#    this. For the moment we have kludge, pending a fix in las.
#  - [TODO] the [cacao] stats line are not provided in cado, so we can't
#    do much with them.
#  - [TODO] there's some duplication with the scan-holes code. We should
#    think about refactoring this code somehow.

sub usage {
    print <<EOF;
Usage:
    merge.pl <input directory> -o <output directory> -n <name> [options]
Options:
    -s <number>         split in ranges of given width (default 1000000)
    --range <q0>-<q1>   treat only files within range <q0>-<q1>
    -f                  allows no relation for one sq

Note that existing files in the output directory which match the names of
expected result files will result in the program doing a fast-forward
over the given range, thereby not re-computing the output files.

Also note that this script is *not* intended to be used with identical
input and output directories.
EOF
    exit 1;
}

use strict;
use warnings;
# use Data::Dumper;
use File::Basename;
use Cwd qw(abs_path);

my $wdir = abs_path(dirname(dirname($0)));
my $scriptsdir = $0;
# dirname:
$scriptsdir =~ s{/[^/]+$}{} || do { $scriptsdir = '.'; };

# my $nextprime = "$scriptsdir/gmp_nextprime";
my $roots = "$wdir/bin/roots";
$roots = "$wdir/scripts/roots" unless -x $roots;
$roots = "$scriptsdir/roots" unless -x $roots;

# die "Need nextprime program ; could not find $nextprime" unless -x $nextprime;
die "Need roots program ; could not find `$wdir/bin/roots'" unless -x $roots;

my @input_dirs;
my $output_dir;
my $name;
my $forced_read = 0;

my $slice = 1000000;

my $forced_range;

my @ARGV_orig=@ARGV;

if (scalar @ARGV == 4 && $ARGV[0] =~ /^\d+$/ && $ARGV[1] =~ /^\d+$/) {
    # Be compatible with merge_results.pl's command line.
    $forced_range = [ ];
    push @$forced_range, shift @ARGV;
    push @$forced_range, shift @ARGV;
    push @input_dirs, shift @ARGV;
    $output_dir = shift @ARGV;
}

while (defined(my $x = shift @ARGV)) {
    if ($x eq '-s') {
        $slice = shift @ARGV or die "$x needs an argument\n";
        next;
    } elsif ($x eq '--range') {
        my $range = shift @ARGV or die "$x needs an argument\n";
        $range =~ /^(\d+)-(\d+)$/ or die "$x takes sssssss-eeeeeee\n";
        $forced_range = [ $1, $2 ];
        next;
    } elsif ($x eq '-o') {
        $output_dir = shift @ARGV or die "$x needs an argument\n";
        next;
    } elsif ($x eq '-n') {
        $name = shift @ARGV or die "$x needs an argument\n";
        next;
    } elsif ($x eq '-f') {
        $forced_read = 1;
        next;
    } elsif (-d $x) {
        push @input_dirs, $x;
        next;
    } else {
        die "Unparsed argument $x\n";
    }
}

my $file_base = "$name.rels";
my $poly_file = "$wdir/$name.poly";
$poly_file = "$scriptsdir/$name.poly" unless -e $poly_file;

# die "Need nextprime program ; could not find $nextprime" unless -x $nextprime;
die "Need poly file ; could not find `$wdir/$name.poly'" unless -e $poly_file;

my $file_pattern = qr/^(?:.*\/)?$file_base.(\d+)-(\d+)(?:\.gz)?$/;

# Here is a trimmed-down version of scan-holes.pl. We bark whenever
# anything not-contiguous shows up.
sub check_range {
    my @ranges=();

    for my $d (@_) {
        opendir my $dh, $d or die "$d: $!";
        for my $f (readdir $dh) {
            next if $f =~ /^\./;
            if ($f =~ $file_pattern) {
                if (defined($forced_range)) {
                    next if $1 >= $forced_range->[1];
                    next if $2 <= $forced_range->[0];
                }
                push @ranges, [ $1, $2, "$d/$f" ];
            }
        }
        close $dh;
    }

    @ranges = sort { $a->[0] <=> $b->[0] } @ranges;

    my $current = shift @ranges or return undef;

    my $err=0;
    my ($cs,$ce) = @$current;
    for my $x (@ranges) {
        # $cs,$ce are always defined here.
        my ($s, $e) = @$x;
        die "Uh: $e <= $s" if $e <= $s;
        if ($s == $ce) {
            $ce = $e;
        } elsif ($s < $ce) {
            print STDERR "Overlap with $x->[2]\n";
            $err=1;
        } else {
            print STDERR "Hole before $x->[2]\n";
            $err=1;
            $ce = $e;   # Arrange in order to silence further warnings.
        }
    }
    die "Please fix errors before proceeding" if $err;

    unshift @ranges, $current;

    return ($cs, $ce, \@ranges);
}

die "Missing input dir ; ARGV=@ARGV_orig" unless scalar @input_dirs;
die "Missing output dir ; ARGV=@ARGV_orig" unless -d $output_dir;

my ($start, $end, $ranges) = check_range @input_dirs;

die "Empty file list" unless scalar @$ranges;

if (defined($forced_range)) {
    die "forced range not included within [$start..${end}["
        unless $start <= $forced_range->[0]
        && $forced_range->[0] <= $forced_range->[1]
        && $forced_range->[1] <= $end;
}

# These are global variables used within the loop.
my $data_fh;
## (stats disabled) my $stats_fh;
my $q0=$start;
my $q1;

my $q;
my $r;
my $nrels;
## (stats disabled) my $stats;

my $fh;
## (stats disabled) my $sh;

my ($rmin,$rmax);
my $ifile;
## (stats disabled) my $sfile;

my $last_done;

my @lookahead = ();

# look ahead somewhat aggressively. It's much, much cheaper than a
# fork+exec anyway, so better arrange so that we have some work to do.
sub provision_qr {
    my $q = shift;
    my @r;
    my $r;
    do {
        chomp($q);
        my $q1 = $q + 10000;
        @r = `$roots -poly $poly_file -q0 $q -q1 $q1`;
        $q = $q1;
    } until scalar @r;
    for my $r (@r) {
        chomp($r);
        my @rr = split ' ', $r;
        $q = shift @rr;
        @rr = sort { $a <=> $b} @rr;
        for my $rr (@rr) {
            push @lookahead, [$q, $rr];
        }
    }
}

sub expand_qr {
    provision_qr $lookahead[$#lookahead]->[0];
}

sub forthcoming_qr {
    my $i = shift @_ || 0;
    while (scalar @lookahead <= $i) {
        expand_qr;
    }
    return $lookahead[$i] or die;
}

sub pop_qr {
    my $h = shift @lookahead;
    if (scalar @lookahead == 0) {
        provision_qr $h->[0];
    }
}

sub qr_cmp {
    my ($qr0, $qr1) = @_;
    return $qr0->[0] <=> $qr1->[0] ||  $qr0->[1] <=> $qr1->[1];
}

sub qr_checknext {
    my ($qr0) = @_;

    # FIXME. 20101215. This is a kludge here. Since we don't sort the
    # roots in las, it's more difficult to predict the output. We resort
    # to the following. Pop from @lookahead all ideals with matching q
    # into a hash.  If one matches the given r, then pop it for real, and
    # re-push the rest. It's not very efficient, but has the advantage of
    # keeping code changes down to a minimum.
    my $rs={};
    while ((my $qr1 = forthcoming_qr)->[0] == $qr0->[0]) {
        $rs->{$qr1->[1]}=$qr1;
        pop_qr;
    }
    if ($forced_read == 0 && scalar keys %$rs == 0) {
        die "Unexpected q,r ($qr0->[0], $qr0->[1]). No more ideals above $qr0->[0]. Maybe a discarded ideal ?\n";
    }

    if (exists($rs->{$qr0->[1]})) {
        delete $rs->{$qr0->[1]};
        unshift @lookahead, values %$rs;
        return;
    } else {
        if ($forced_read) {
          return;
        } else {
          die "Unexpected q,r : $qr0->[0],$qr0->[1])\n";
        }
    }
    # (resume normal code). Never reached for the moment, since the
    # kludge takes over.
    my $qr1 = forthcoming_qr;
        if (qr_cmp($qr0, $qr1) == 0) {
        pop_qr;
        return;
    }


    print STDERR "Unexpected q,r : expected ($qr1->[0],$qr1->[1]), not ($qr0->[0],$qr0->[1])\n";
    my $qr2 = forthcoming_qr(1);
    if (qr_cmp($qr0, $qr2) == 0) {
        print STDERR "Assuming discarded ideal ($qr1->[0],$qr1->[1])\n";

## (stats disabled)         if (defined($sh)) {
## (stats disabled)             my $d;
## (stats disabled)             if (!defined($d = <$sh>)) {
## (stats disabled)                 die "Cannot read stat data from $sfile";
## (stats disabled)             }
## (stats disabled)             die "Discard line not found in $sfile ; found\n$d"
## (stats disabled)                 unless $d =~ /^# \[cacao\] Discard $qr1->[0] $qr1->[1]/;
## (stats disabled)             print $stats_fh $d;
## (stats disabled)         } else {
## (stats disabled)             print $stats_fh "# [cacao] Discard $qr1->[0] $qr1->[1]\n";
## (stats disabled)         }
        pop_qr;
        pop_qr;
    } else {
        print STDERR "Failed comparison with next-to-come ($qr2->[0],$qr2->[1])\n";
        die;
    }
}



sub open_new_files {
    $q1 = $q0 + $slice;
    $q1 -= ($q1 % $slice);
    if ($q1 > $end) {
        $q1 = $end;
    }

    my $iname = "$file_base.$q0-$q1";

    print STDERR "Saving $output_dir/$iname.gz\n";
    
    ## (stats disabled) my $sname = "stats.$q0-$q1";

    open $data_fh,  "| gzip --best > $output_dir/tmp-$iname.gz";
    ## (stats disabled) open $stats_fh, "| gzip --best > $output_dir/tmp-$sname.gz";
    # if not splitting, collide $stats_fh and $data_fh
}

sub close_current_files {
    return unless defined $data_fh;
    close $data_fh;
    ## (stats disabled) close $stats_fh;
    my $iname = "$file_base.$q0-$q1";
    ## (stats disabled) my $sname = "stats.$q0-$q1";
    # If not splitting, do not close stats_fh of course.
    rename "$output_dir/tmp-$iname.gz", "$output_dir/$iname.gz";
    ## (stats disabled) rename "$output_dir/tmp-$sname.gz", "$output_dir/$sname.gz";
}

sub may_fast_forward {
    if ($q0 >= $end) {
        return 0;
    }

    my $lq1 = $q0 + $slice;
    $lq1 -= ($lq1 % $slice);
    if ($lq1 > $end) {
        $lq1 = $end;
    }
    my $iname = "$file_base.$q0-$lq1";
    ## (stats disabled) my $sname = "stats.$q0-$lq1";

    # yields undefined on failure.
    ## (stats disabled) my $ans =  -f "$output_dir/$iname.gz" && -f "$output_dir/$sname.gz";
    my $ans =  -f "$output_dir/$iname.gz";
    return $ans;
}


sub fast_forward {
    my $oq0 = $q0;
    my $tskip=0;
    while ($q0 < $end) {
        $q1 = $q0 + $slice;
        $q1 -= ($q1 % $slice);
        if ($q1 > $end) {
            $q1 = $end;
        }
        my $iname = "$file_base.$q0-$q1";
        ## (stats disabled) my $sname = "stats.$q0-$q1";


        last unless -f "$output_dir/$iname.gz";
        ## (stats disabled) last unless -f "$output_dir/$sname.gz";

        my $skip=0;
        while (scalar @$ranges && $ranges->[0]->[1] <= $q1) {
            my $skipped = shift @$ranges;
            # print STDERR "Skipped $skipped->[2]\n";
            $skip++;
        }
        print STDERR "$q0-$q1 done ; skipped $skip files   \r";
        $tskip+=$skip;
        $q0 = $q1;
        @lookahead = ();
        provision_qr $q0;
    }
    if ($q0 > $oq0) {
        print STDERR "Skipped $oq0-$q0 (already done, total $tskip files)\n";
        if ($q0 >= $end) {
            exit 0;
        }
    }
}

sub finish_qr  {
    # We've just done (q,r).

    return unless defined $q;

    #print "Finished ($q,$r) ; $nrels rels\n";

    # There is a separate stats file, read from it.
## (stats disabled)     if (defined($sh)) {
## (stats disabled)         #print "Recovering stats from $sfile\n";
## (stats disabled)         die if defined $stats;
## (stats disabled)         if (!defined($stats = <$sh>)) {
## (stats disabled)             die "Cannot read stat data from $sfile";
## (stats disabled)         }
## (stats disabled)         my $errmsg="bad ($q,$r) stats line from $sfile";
## (stats disabled)         $stats =~ /^\#\s\[cacao\]\sStat:\s
## (stats disabled)             (\d+)\s(\d+)\s     # q,r
## (stats disabled)             \d+\s[\d\.]+\s     # rels, time
## (stats disabled)             [\/\w\.-]+\s[\w\.-]+\s # binary, node
## (stats disabled)             (?:grid5000.fr\s)? # Accomodate bug...
## (stats disabled)             \d+                # timestamp
## (stats disabled)         $/x
## (stats disabled)             or die "$errmsg;\n$stats";
## (stats disabled)         if ($1 ne $q || $2 ne $r) {
## (stats disabled)             die "$errmsg;\n$stats";
## (stats disabled)         }
## (stats disabled) 
## (stats disabled)         print $stats_fh $stats;
## (stats disabled)     }
## (stats disabled) 
## (stats disabled)     if (!defined($stats)) {
## (stats disabled)         $stats = "# [cacao] Stat: $q $r $nrels\n";
## (stats disabled)         print $stats_fh $stats;
## (stats disabled)     }

    $last_done = [$q, $r];

    $q = undef;
}

# my $starting_pattern = qr/^# (?:Start|\[cacao\] Discard) (\d+) (\d+)/;
my $starting_pattern = qr/^# Sieving q=(\d+); rho=(\d+);/;
## (stats disabled) my $cacao_regexp = qr/^# (?:\[cacao\]|Stat:)/;
## (stats disabled) my $cacao_stats = qr/^# (\[cacao\]\s+)?(Stat:\s+)(\d+\s+[\d\.]+\s+[\/\w\.-]+\s+[\w\.-]+)(?:\s+grid5000.fr)?(\s+\d+)?$/;




if (defined($forced_range)) {
    $q0 = $forced_range->[0];
    $end = $forced_range->[1];

    my $skip;
    $skip=0;
    while (scalar @$ranges && $ranges->[0]->[1] <= $q0) {
        shift @$ranges;
        $skip++;
    }
    print STDERR "Skipped $skip files at beginning of range\n";
    $skip=0;
    while (scalar @$ranges && $ranges->[$#$ranges]->[0] >= $end) {
        pop @$ranges;
        $skip++;
    }
    print STDERR "Skipped $skip files at end of range\n";
}

provision_qr $q0;

fast_forward;

if ($q0 >= $end) {
    print STDERR "$start-$end : Nothing to do\n";
    exit 0;
}

open_new_files;

FILE_LOOP:
while (defined(my $me = shift @$ranges)) {
    ($rmin,$rmax,$ifile) = @$me;

    ## (stats disabled) $sfile = $ifile;
    ## (stats disabled) $sfile =~ s/$file_base/stats/;
    ## (stats disabled) $sfile .= '.gz' unless $sfile  =~ /\.gz$/;


    $fh = undef;
    ## (stats disabled) $sh = undef;

    if ($ifile =~ /\.gz$/) {
        open $fh, "gzip -dc $ifile |" or die "$ifile: $!";
        ## (stats disabled) if (-f $sfile) {
            ## (stats disabled) open $sh, "gzip -dc $sfile |" or die "$sfile: $!";
            ## (stats disabled) }
    } else {
        open $fh, $ifile or die "$ifile: $!";
    }

    my $x = <$fh>;

    my $filetime = (stat $ifile)[9];

    print STDERR "$ifile\r";

    if ($rmin < $q0) {
        print STDERR "\nFast-forwarding until $q0\n";
        die unless $rmax > $q0;
        while (1) {
            $x = <$fh>;
            last if !defined($x);
            next unless $x =~ $starting_pattern;

            $q = $1;
            $r = $2;

            if ($q < $q0) {
                ## (stats disabled) <$sh> if defined($sh);
                next;
            }

            # Got it.
            eval { qr_checknext [$q, $r]; };
            die "$ifile: $@" if $@;
            $nrels = 0;
            ## (stats disabled) $stats = undef;
            # If not splitting, print the Discard line as well.
            print $data_fh $x unless $x =~ /Discard/;
            last;
        }
        if ($q < $q0) {
            # There's an infelicity in the choice of range ends, which
            # might lead to files which _appear_ to contain the first
            # prime above some bound, while the do not.
            # e.g. 2543649612-2543650098
            print STDERR "wrong starting file $ifile -- skipping to next\n";
            close $fh;
            ## (stats disabled) close $sh if defined $sh;
            $q = undef;
            next FILE_LOOP;
        }
        print STDERR "OK -- now at ($q,$r)\n";
    }

    # Loop through all q's in the file.
    while (1) {
        $x = <$fh>;
        
        if (!defined($x)) {
            # This file is over.
            # Mark the current special q as finished.
            finish_qr;

            my $qr1 = forthcoming_qr;
            if ( $forced_read == 0 ) {
              die "Missing ($qr1->[0], $qr1->[1]) in $ifile"
                 if $qr1->[0] < $rmax;
            }
            last;
        }

        # print "Got $x";

        if ($x =~ $starting_pattern) {
            my $newq = $1;
            my $newr = $2;

            # The previous prime is over.
            finish_qr;

            if ($newq > $rmax) {
                # Provision for handling extraneous data at the end of
                # merged files.
                my $qr1 = forthcoming_qr;
                print STDERR "Special-q ($newq,$newr) does not belong to $ifile\n";
                if (qr_cmp($qr1, [$newq,$newr]) == 0) {
                    print STDERR
                        "This seems to be a false duplicated"
                        . " as sometimes created by merge_results.pl\n";
                } else {
                    die "This is really unexpected. Please fix $ifile\n";
                }
                while (defined($x=<$fh>)) {
                    die "Too much extraneous data at the end of $ifile: $x"
                        if $x =~ $starting_pattern;
                }
                close $fh;
                ## (stats disabled) close $sh if defined $sh;
                $q = undef;
                next FILE_LOOP;
            }

            if ($newq > $end) {
                # Happens only at the very end, and even then, only when
                # handling a restricted range of the input data
                print STDERR "Argh, there is some remaining data"
                    if scalar @$ranges;
                close $fh;
                ## (stats disabled) close $sh if defined $sh;
                last FILE_LOOP;
            }

            eval { qr_checknext [$newq, $newr]; };
            die "$ifile: $@" if $@;

            if ($newq > $q1) {
                close_current_files;
                print STDERR "\n";
                $q0 = $q1;
                if (may_fast_forward) {
                    print STDERR "Interrupting $me->[2] for fast-forwarding\n";
                    unshift @$ranges, $me;
                    close $fh;
                    ## (stats disabled) close $sh if defined $sh;
                    fast_forward;
                    open_new_files;
                    next FILE_LOOP;
                }
                open_new_files;
            }

            $q = $newq;
            $r = $newr;
            $nrels = 0;
            ## (stats disabled) $stats = undef;

            # If not splitting, print the Discard line as well.
            print $data_fh $x unless $x =~ /Discard/;

            next;
        }

## (stats disabled)         if ($x =~ $cacao_regexp) {
## (stats disabled)             if ($x =~ $cacao_stats) {
## (stats disabled)                 my ($prefix, $mnem, $data, $stamp) = ($1,$2,$3,$4);
## (stats disabled)                 $prefix = "[cacao] ";
## (stats disabled)                 if (!defined($stamp)) { $stamp=" $filetime"; }
## (stats disabled)                 # If not splitting, do not add q,r
## (stats disabled)                 $data = "$q $r " . $data;
## (stats disabled)                 $stats = "# " . $prefix . $mnem . $data . $stamp . "\n";
## (stats disabled) 
## (stats disabled)                 print $stats_fh $stats;
## (stats disabled)                 next;
## (stats disabled)             } else {
## (stats disabled)                 die "Bad [cacao] line: $x\n";
## (stats disabled)             }
## (stats disabled)         }

        if ($x =~ /^#/) {
            # comment lines are passed unchanged
        } else {
            $nrels += ($x =~ /^[^#]/);
        }

        print $data_fh $x;
    }
    close $fh;
    ## (stats disabled) close $sh if defined $sh;
}

finish_qr;
print STDERR "\n";

my $qr1 = forthcoming_qr;
if ( $forced_read == 0 ) {
  die "Unfinished files (want ($qr1->[0],$qr1->[1]) after ($last_done->[0],$last_done->[1]) )"
     if $qr1->[0] < $q1;
}

close_current_files;
