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

open(IN, "<", $infile) or die "I cannot open ".$infile.", sorry!";
#open(OUT, ">>", $outfile) or die "I cannot open".$outfile.", sorry!"; # to append
open(OUT, ">", $outfile) or die "I cannot open ".$outfile.", sorry!";

print OUT "Type \t Impact \t Chromosome \t Position \t Ref \t Alt \t geneID \t NuclVar \t ProtVar \t Rank \n";
while(<IN>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^#/) {next;}
    my @Cols = split('\t', $line);
    my $chrom = $Cols[0];
    my $Pos = $Cols[1];
    my $RefNT = $Cols[3];
    my $AltNT = $Cols[4];
    my $Info = $Cols[7];   ### This contains the ENTIRE <INFO> field for the vcf.
                           ### For reference, the INFO field of vcf files can be gnarly --
                           ### It is a semicolon-delimited list of the defined INFO types (defined in the VCF header)
                           ### Which all start with "ABBREV=".
    my $annotations = $Info =~ s/^.*ANN=(.*$)/$1/r; # I can get away with just isolating the ANN field
                                                    # like this because it's the last entry in my annotated VCF, but
                                                    # if it were not, we could try splitting the field on semicolons
                                                    # and looping through it to find the entry that started with "ANN="
### The entry in the ANN field is a comma-delimited list of variant-annotations; each annotation is a
### pipe-delimited list of values produced by snpEff. (It has to be pipe-delimited; semicolons, commas, underscores, and dots are already used.)
### ~~~~ This should be a json-format file. You know it, I know it, and the people of the world know it.
### ~~~~ But we're not going to change entrenched data formats -- we just have to work with them.
###
### Here are some example annotations:
###annotations: T|missense_variant|MODERATE|jg23759|jg23759|transcript|jg23759.t1|protein_coding|4/4|n.529G>A|p.Glu177Lys|529/648|529/-1|177/-1||,
###             T|downstream_gene_variant|MODIFIER|65307|65307|transcript|65307_t|protein_coding||n.*349C>T|||||349|,
###             T|downstream_gene_variant|MODIFIER|jg23758|jg23758|transcript|jg23758.t1|protein_coding||n.*349C>T|||||349|
###
    my @annots = split(',', $annotations); # split the annotations for this line on commas
    my $count=0;                           # start counting to keep track of annotation rank order
    foreach(@annots) {
        $count++; # increment the count (the first one is "1"!)
        my $variant = $_ =~ s/\|/;/gr;  # Replace the pipe delimiter with a semicolon delimiter, because
                                        # splitting on the pipe is tricky (because pipes mean "or" in regular expressions)
        my @variant = split(';', $variant); # Now we split on the semicolon
        my $geneID = $variant[3];
        my $type = $variant[1];
        my $impact = $variant[2];
        my $n = $variant[9];
        my $p = $variant[10];
        # For right now, I'm printing to the screen what's also going to the outfile, to simplify debugging
        my $output = $type."\t".$impact."\t".$chrom."\t".$Pos."\t".$RefNT ."\t".$AltNT ."\t".$geneID ."\t".$n ."\t".$p."\t".$count."\n";
        print $output;
        print OUT $output;
    }
}



##
##

close(IN);
close(OUT);
