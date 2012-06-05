#!/usr/bin/perl

# I have this in /etc/apache2/conf.d/refresh.conf:
# ScriptAlias /refresh.pl /var/www.cgi/refresh.pl
# <Directory "/var/www.cgi">
# AllowOverride None
# Options +ExecCGI -MultiViews +SymLinksIfOwnerMatch
# Order allow,deny
# Allow from all
# </Directory>

# And this in /git/cado-nfs/cado-nfs.git/hooks/post-receive on gforge:
# URL=http://fondant.loria.fr/refresh.pl
# post=""
# while read r1 r2 x ; do
#         if ! [ "$post" ] ; then
#                 post="h=$x"
#         else
#                 post="${post};h=$x"
#         fi
# done
# wget "--post-data=$post" -q $URL -O- > /dev/null 2>&1

# And finally /localdisk/pub/cado-nfs.git/ contains a bare repo with the proper
# ./hooks/post-receive script.

# For the gitweb part: I simply need $projectroot = "/localdisk/pub"; in /etc/gitweb.conf


 
use strict;
use warnings;
use Data::Dumper;

use CGI;

my $q = CGI->new;

print $q->header();

# print "Hello, world\n";

$ENV{'GIT_DIR'}="/localdisk/pub/cado-nfs.git/";
$ENV{PATH}='/bin:/usr/bin';

open STDOUT, ">&/dev/null";
open STDERR, ">&/dev/null";

sub rev_parse {
	my $h = shift;
	my $x=`git rev-parse $h 2>/dev/null`;
	chomp($x);
	if ($x=~/^([0-9a-f]{40})$/) {
		$x=$1;
	} else {
		$x='0' x 40;
	}
	return $x;
}

my @current=();
my @heads = $q->param('h');
for my $h (@heads) {
	# print "Head: $h\n";
	push @current,[$h,rev_parse($h)];
}

system "git fetch";

my $payload="";

for my $hr (@current) {
	my ($h,$r) = @$hr;
	my $nr = rev_parse($h);
	if ($r ne $nr) {
		$payload .= "$r $nr $h\n";
	}
}

if (length($payload)) {
	open F, "| $ENV{'GIT_DIR'}/hooks/post-receive";
	print F $payload;
	close F;

}
