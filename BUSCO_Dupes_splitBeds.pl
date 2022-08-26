#!/usr/bin/env perl
##############################################################################################
# Usage: perl SNPAnnotationsToTable_all.pl <annotated vcf file> <tab-delimited output file>  #
# This parses a VCF file annotated by snpEff and exports the values of the ANN INFO field to #
# a tab-delimited file. If one SNP has multiple annotations appended to it, this script will #
# create a row for each annotation, repeating the coordinate and allele information, and will#
# enumerate the "rank" of each annotation in order from first to last.                       #
##############################################################################################
use warnings;
use strict;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $outfile2 = $ARGV[2];

open(IN, "<", $infile) or die "I cannot open ".$infile.", sorry!";
#open(OUT, ">>", $outfile) or die "I cannot open".$outfile.", sorry!"; # to append
open(OUT1, ">", $outfile) or die "I cannot open ".$outfile.", sorry!";
open(OUT2, ">", $outfile2) or die "I cannot open ".$outfile2.", sorry!";

## We want to take the input bed file and, wherever column[3] is the same, we want the first instance to go in Bedfile 1 and the second instance to go in Bedfile 2 
#
#Sample lines:
##1       22466906        22468147        Name=19562at7147.faa
##1       22727536        22729424        Name=19562at7147.faa
##1       23310551        23424745        Name=23840at7147.faa
##1       23423117        23424327        Name=23840at7147.faa
##1       26187113        26249822        Name=49533at7147.faa
##1       26757484        26800677        Name=49533at7147.faa
##1       30738518        30750528        Name=10627at7147.faa
##1       31426272        31449070        Name=10627at7147.faa

my %bed1 = ();
my %bed2 = ();

my $current = "";
while(<IN>) {
	my $line = $_;
	chomp $line;
	#print $line."\n";
	my @Cols = split('\t', $line);
	#$Cols[0] = chromosome, $Cols[1] = start, $Cols[2] = end, $Cols[3] = name
	if (exists($bed1{$Cols[3]})) {
#		print "Found in bed1:" . $Cols[3] ." = ". $bed1{$Cols[3]} ." \n ";
		$bed2{$Cols[3]} = $line;
	} else {
#		print "Not in bed1:" . $Cols[3] ."\n";
		$bed1{$Cols[3]} = $line;
	}
}
close(IN);

#my $keysbed1 = keys %bed1;
#print scalar($keysbed1);

foreach my $entry (keys %bed1) {
	print $entry .": ". $bed1{$entry} ."\n";
        print OUT1 $bed1{$entry} ."\n";
}

close(OUT1);

foreach my $entry2 (keys %bed2) {
	print $entry2 .": ".$bed2{$entry2} ."\n";
        print OUT2 $bed2{$entry2} . "\n";
}
close(OUT2);

