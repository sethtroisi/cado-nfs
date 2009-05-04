#!/usr/bin/perl -w

# This script makes sure that the files.dist and files.nodist cover the 
# set of files in the source tree (as defined by the SCM system of course).
# Note that nothing is there to make sure that the files which are _required_
# for compiling the tree are actually there. That's admittedly unfortunate,
# but doing it the right way would require some cooperation from the build
# system. At the moment, CMake has no functionality to build source package
# (worse, it's got a very lame excuse for it called CPack).

# This script must be called from the top of the tree.

use strict;
use warnings;

die "Please call $0 from the top of the source tree" unless -f "cado.h";

my @files_scm;

system "git rev-parse HEAD > /dev/null 2>/dev/null";
if ($? == 0) {
    open SCM, "git ls-files |";
    @files_scm = <SCM>;
    close SCM;
} elsif (-d ".svn") {
    open SCM, "svn list -R |";
    @files_scm = <SCM>;
    close SCM;
}

chomp for @files_scm;

die "No SCM found ???" unless @files_scm;

my $known_lists=[];

sub add_known_list {
    my $n = shift;
    my $f = {};
    my $d = {};
    open my $fh, $n or return;
    while (defined($_=<$fh>)) {
        next unless /^[^#]/;
        chomp($_);
        m{/$} && do { $d->{$_}=0; next; };
        $f->{$_}=0;
    }
    close $fh;
    if (scalar keys %$d == 0 && scalar keys %$f == 0) {
        return;
    }
    return [ $n, $f, $d ];
}

push @$known_lists, &add_known_list('files.dist');
push @$known_lists, &add_known_list('files.nodist');
if (defined(my $x = &add_known_list('files.unknown'))) {
    warn "Warning: files.unknown is not-empty. Please arrange so that every file is covered by a pattern in either files.dist or files.nodist";
    push @$known_lists, $x;
}

my $error=0;
for my $file (@files_scm) {
    # Check whether this matches in any way.
    my @matches=();
    for my $k (@$known_lists) {
        my ($n,$f,$d) = @$k;
        if (exists($f->{$file})) {
            $f->{$file}++;
            push @matches, "$n:$file";
        }
        for my $p (keys %$d) {
            if ($file =~ /^$p/) {
                $d->{$p}++;
                push @matches, "$n:$p";
            }
        }
    }
    my $nm = scalar @matches;
    next if ($nm == 1);
    if ($nm == 0) {
        print STDERR "File $file is not covered by any pattern\n";
        $error++;
    } else {
        print STDERR "File $file is covered by $nm patterns: @matches\n";
        $error++;
    }
}

for my $k (@$known_lists) {
    my ($n,$f,$d) = @$k;
    for (keys %$f) {
        next if $f->{$_};
        print STDERR "Pattern $n:$_ matches nothing\n";
        $error++;
    }
    for (keys %$d) {
        next if $d->{$_};
        print STDERR "Pattern $n:$_ matches nothing\n";
        $error++;
    }
}
if ($error) {
    die "$error errors found. Please fix files.dist and files.nodist";
}

