#!/usr/bin/env perl

# Author: JIN Xiaoyang
# Date  : 2015-11-11
#
# Find specified pattern in a genome fasta file, the pattern is a simplified
# regular expression (with `.` and `[]` supported), while the matched sequences
# can cross two lines.

use warnings;
use strict;

my $fa_file="/Users/jinxy/Database/hg38/genome.fa";

sub usage()
{
    print "Usage: $0 <pattern> [fasta_file]\n";
    print "\te.g. $0 'ggaa..gaaa' genome.fa\n";
}

sub check_pattern()
{
    my $str=$_[0];
    # print "pattern = '$str'\n";
    if ($str =~ m/^[agctn.\[\]^]+$/i) {
	return $str;
    } else {
	return '';
    }
}

sub get_pattern_len()
{
    my $str=$_[0];
    $str =~ s/\[[^\[\]]+\]/./g;
    # print "pattern = '$str'\n";
    return length $str;
}

&usage and exit if not defined $ARGV[0];
$fa_file = $ARGV[1] if defined $ARGV[1];

my $pattern = &check_pattern($ARGV[0]);
if ($pattern eq '')
{
    print STDERR "Error: invalid pattern\n";
    exit 1;
}

my $pat_len =  &get_pattern_len($pattern);

open(FA, "<$fa_file") or die "Couldn't open $fa_file for reading: $!";

my $chr='';
my $chr_count=0;
my $tail_str='';

while(<FA>) {
    chomp;
    /^\s*$/ and next;

    if(/^>(\w+)/) {
	print "$chr\t$chr_count\n" if $chr ne '';
	$chr=$1;
	$chr_count = 0;
    } else {
	my @word=(($tail_str.$_)=~m/$pattern/ig);
	my $tail_str=substr($_, 1-$pat_len);
	$chr_count += scalar(@word);
    }
}

print "$chr\t$chr_count\n" if $chr ne '';
