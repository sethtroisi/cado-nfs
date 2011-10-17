#!/usr/bin/perl -Tw

# This is an example git post-receive hook script. It can be used to
# automatically send emails listing new revisions to the repository
# introduced by the change being reported.
#
# Copyright 2011 Emmanuel Thomé -- This script is placed in the public domain.
#
# Usage: copy this to $GIT_DIR/hooks/post-receive (executable) in the
# repository receiving pushes.
#
# The file $GIT_DIR/config should contain the following keys (these are
# examples).
#
# hooks.emailprefix=
# hooks.mailinglist=Cado-NFS commits <cado-nfs-commits@lists.gforge.inria.fr>
# hooks.envelopesender=Cado-NFS commits <cado-nfs-commits@lists.gforge.inria.fr>
# hooks.gitweb=http://tomato.elliptic.org/gitweb/?p=cado-nfs.git;

use strict;
use warnings;
use UUID;       # Or any simple way of creating MIME delimiters, really.
use File::Basename;
use POSIX;

$ENV{'PATH'}="/bin:/usr/bin:/usr/sbin";

my $git_dir="$ENV{'GIT_DIR'}" or die 'no $GIT_DIR set';

my $repo = basename($git_dir, ".git");

sub config_get {
    my $key = shift;
    my $x = `git config $key 2>/dev/null`;
    if ($x =~ /^([^'`]*)$/) {
        # A positive answer with empty content still reaches here.
        $x=$1;
        chomp($x);
        return $x;
    }
    my $default = shift;
    return '' unless $default;
    die "Mandatory config key $key is missing" if $default eq '!';
    return $default;
}
    
my $gitweb=config_get('hooks.gitweb');
chomp($gitweb);
my $recipient=config_get('hooks.mailinglist', '!');
my $sender=config_get('hooks.envelopesender', '!');
my $subj_prefix=config_get('hooks.emailprefix', "[$repo-commits] ");

sub commit_link {
    my $commit = shift;
    return $commit unless $gitweb;
    my $url = $gitweb . "a=commitdiff;h=$commit";
    return "<a href=\"$url\">$commit</a>";
}

my @branches_touched=();
my $total_text_plain='';
my $total_text_html='';
my $total_nc = 0;
my $last_subject;

sub boundary { UUID::generate($_);UUID::unparse($_, my $x); $x; }

my $mainbound = boundary();

my $patches_mboxes='';

while (defined($_=<STDIN>)) {
    chomp($_);
    s/^\s*//;
    next if /^$/;
    s/^([0-9a-f]{40})\s+// or die "Bad input line: $_";
    my $old = $1;
    s/^([0-9a-f]{40})\s+// or die "Bad input line: $_";
    my $new = $1;
    next if $old =~ /^0+$/;
    next if $new =~ /^0+$/;
    my $ref = $_;
    $ref =~ s,^refs/heads/,,;
    push @branches_touched, $ref;
    my $format = '%h %ae %s';
    open F, "git log --pretty=\"$format\" --abbrev $old..$new |";
    my $text_plain='';
    my $text_html='';
    my $author='';
    my $nc=0;
    while (defined(my $line = <F>)) {
        my ($commit, $who, $subject) = split(' ', $line, 3);
        chomp($subject);
        if ($who ne $author) {
            $text_html .= "</ul>\n" if $author;
            $author = $who;
            $text_plain .= "\n$who\n";
            $text_html .= "<p>$who</p>\n<ul>\n";
        }
        $text_plain .= " $commit $subject\n";
        my $link = commit_link($commit);
        $text_html .= "<li>$link $subject</li>\n";
        $last_subject = $subject;
        $nc++;
    }
    close F;
    $text_html .= "</ul>\n" if $author;
    $text_plain .= "\n\n";

    my $x = "Branch: $ref";
    $text_plain = "$x, $nc new commits\n" . ('=' x length($x)) . "\n" . $text_plain;
    $text_html = "<h1>Branch: $ref, $nc new commits</h1>\n" . $text_html;;

    $total_text_plain .= $text_plain;
    $total_text_html .= $text_html;
    $total_nc += $nc;

    my $bound = boundary();

    open PATCH, "git format-patch --stdout $old..$new |";
    $x = eval { local $/=undef; <PATCH>; };
    $x =~ s/^(From (\w{7})\w{33})/--$bound
Content-Type: message\/rfc822; charset=us-ascii
Content-Disposition: attachment; filename="$2.patch"

$1/gm;

    my $now = POSIX::strftime( "%Y%m%d%H%M%S", localtime);

    $patches_mboxes .= <<EOF;
--$mainbound
Content-Type: multipart/mixed; boundary=$bound
Content-Disposition: attachment; filename="$repo-$ref-$now.mbox"

$x

--$bound--
EOF
    close PATCH;
}


my $refs = join(', ', @branches_touched);

return unless $total_nc;

my $subject = "$total_nc new commits";
if ($total_nc == 1) {
    $subject = $last_subject;
}


my $header = <<EOF;
To: $recipient
From: $sender
Subject: $subj_prefix($refs) $subject
EOF

$header=~s/\n$//m;

my $textbound = boundary();

my $body = <<EOF;
--$mainbound
Content-Type: multipart/alternative; boundary=$textbound

--$textbound
Content-Type: text/plain; charset=utf-8

$total_text_plain

--$textbound
Content-Type: text/html; charset=utf-8

$total_text_html

--$textbound--

$patches_mboxes

--$mainbound--
EOF

my $bodylength = length($body);
my $bodylines = eval { @_=split(/^/m, $body);scalar @_; };

if (@ARGV && $ARGV[0] eq '-d') {
    open SMTP, ">&STDOUT";
    my $now=strftime("%a %b %d %H:%M:%S %Y", localtime);
    print SMTP "From root $now\n";
} else {
    open SMTP, "| exim '$recipient'";
}

print SMTP <<EOF;
$header
MIME-Version: 1.0
Content-Type: multipart/mixed; boundary=$mainbound
Content-Transfer-Encoding: 8bit
Content-Disposition: inline
Content-Length: $bodylength
Lines: $bodylines

$body
EOF
close SMTP;

