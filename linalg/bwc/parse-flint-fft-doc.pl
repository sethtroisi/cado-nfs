#!/usr/bin/perl

use Data::Dumper;

my $current = undef;
my $all = {};


open F, "<flint-fft/doc/fft.txt";

while (defined($_=<F>)) {
    chomp($_);
    # print STDERR "$_\n";
    if (/^[a-z].*\(/) {
        $current = undef;
        my $key = $_;
        while ($key !~ /\)/) {
            last unless defined($_=<F>);
            chomp($_);
            # print STDERR "$_\n";
            s/^\s*/ /g;
            $key .= $_;
        }
        $all->{$key}="";
        $current = \$all->{$key};
        last unless defined($_=<F>);
        # print STDERR "$_\n";
        die unless /^\s*$/;
        next;
    }
    next unless $current;
    if (/^\*/ || /^\*/) {
        $current = undef;
        next;
    }
    next if /^\s*$/ && $$current eq "";
    $$current .= "$_\n";
}
close F;

my $tags={};
open F, "<flint-fft/tags";
while (defined($_=<F>)) {
    chomp($_);
    next if /^!/;
    /^(\S+)\011(\S+)\011(.*);".*$/ or die;
    push @{$tags->{$1}}, [$2, $3];
}
close F;

# print Dumper($tags);
# print Dumper($all);

for my $header (keys %$all) {
    $header =~ /(.*?)\s+(\S+)\s*\(/ or die;
    my $type = $1;
    my $func = $2;
    my $t = $tags->{$func} or die "No tag for $func";
    if (@$t != 1) {
        die "Multiple tags for $func";
    }
    my ($file, $position) = @{$t->[0]};

    # create a comment block.
    my $doc = $all->{$header};
    my @doclines = split(/^/m, $doc);
    my $prefix = 999;
    for (@doclines) {
        next if /^\s*$/;
        /^(\s*)/;
        my $l = length($1);
        if ($l < $prefix) { $prefix = $l; }
    }

    for (@doclines) {
        s/^\s{$prefix}//;
        s/^/ * /gm;
    }
    $doc = "/*\n" . join("",@doclines) . " */\n";

    print "Adding documentation to flint-fft/$file\n";
    open G, "| ex -s flint-fft/$file > /dev/null";
    print G "set nomagic\n";
    # print G "set noai\n";
    print G "$position\n";
    print G "i\n";
    print G $doc;
    print G ".\n";
    print G "wq\n";
}

