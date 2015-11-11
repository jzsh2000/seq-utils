#!/usr/bin/env perl

# Author: JIN Xiaoyang
# Date  : 2015-11-11
#
# Find specified pattern in a genome fasta file, the pattern is a perl regular
# expression, while the matched pasitions in each sequence will be returned.

use warnings;
use strict;

my $fa_file="/Users/jinxy/Database/hg38/genome.fa";

$/='>';

sub usage()
{
    print "Usage: $0 <pattern> [fasta_file]\n";
    print "\te.g. $0 'ggaa..gaaa' genome.fa\n";
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
	    print $chr."\t".(pos($seq)-length($&)+1)."\t".pos($seq)."\t".$&."\n";
	}
    }
}
