#!/usr/bin/env perl

# Author: JIN Xiaoyang
# Date  : 2015-11-11
#
# Find specified pattern in a genome fasta file, the pattern is a perl regular
# expression, while the matched pasitions in each sequence will be returned.

use warnings;
use strict;

my $fa_file="/Users/jinxy/Database/hg38/genome.fa";
my $reverse=0;

$/='>';

sub usage()
{
    print "Usage: $0 [-r] <pattern> [fasta_file]\n";
    print "\te.g. $0 'ggaa..gaaa' genome.fa\n";
    print "\te.g. $0 -r 'ggaa..gaaa' genome.rev.fa\n";
}

if (defined $ARGV[0] and $ARGV[0] eq '-r') {
    $reverse=1;
    shift;
}

&usage and exit if not defined $ARGV[0];
$fa_file = $ARGV[1] if defined $ARGV[1];
my $pattern = qr/$ARGV[0]/;

open(FA, "<$fa_file") or die "Couldn't open $fa_file for reading: $!";

while(<FA>) {
    if ($_ =~ /\A(\S+)[^\n]*\n(.*)\z/s) {
	my $chr=$1;
	my $seq=$2;

	$seq =~ s/\s//g;
	$seq =~ s/>\z//;

	while ($seq =~ m/$pattern/ig) {
	    my $match_beg = pos($seq) - length($&) + 1;
	    my $match_len = $&;
	    my $match_end = pos($seq);

	    if ($reverse) {
		$match_beg -= (length($seq) + 1);
		$match_end -= (length($seq) + 1);
	    }

	    print "$chr\t$match_beg\t$match_end\t$match_len\n";
	}
    }
}
