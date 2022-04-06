#!/usr/bin/env perl
#
use warnings;
use strict;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open(IN, "<", $infile) or die "I cannot open ".$infile.", sorry!";
#open(OUT, ">>", $outfile) or die "I cannot open".$outfile.", sorry!"; # to append
open(OUT, ">", $outfile) or die "I cannot open ".$outfile.", sorry!";


#my %subCounts = ();

while(<IN>) {
    my $line = $_;
    chomp $line;
    if ($line =~ /^#/) {next;}
    my @Cols = split('\t', $line);
    my $chrom = $Cols[0];
    my $Pos = $Cols[1];
    my $RefNT = $Cols[3];
    my $AltNT = $Cols[4];
    my $Info = $Cols[7];   ## <<< that's where all the chunky goodnes iss
    my $annotations = $Info =~ s/^.*ANN=(.*$)/$1/r;
    #print OUT "annotations: ".$annotations."\n";
    ### Here's what some of those look like:
    ###annotations: C|missense_variant|MODERATE|jg23729|jg23729|transcript|jg23729.t1|protein_coding|5/8|n.3650T>C|p.Val1217Ala|3650/4686|3650/-1|1217/-1||
###annotations: A|missense_variant|MODERATE|jg23729|jg23729|transcript|jg23729.t1|protein_coding|5/8|n.3665G>A|p.Ser1222Asn|3665/4686|3665/-1|1222/-1||A
###annotations: T|missense_variant|MODERATE|jg23759|jg23759|transcript|jg23759.t1|protein_coding|4/4|n.529G>A|p.Glu177Lys|529/648|529/-1|177/-1||,T|downstream_gene_variant|MODIFIER|65307|65307|transcript|65307_t|protein_coding||n.*349C>T|||||349|,T|downstream_gene_variant|MODIFIER|jg23758|jg23758|transcript|jg23758.t1|protein_coding||n.*349C>T|||||349|
###annotations: T|missense_variant|MODERATE|65620|65620|transcript|65620_t|protein_coding|1/4|n.2368G>A|p.Glu790Lys|2368/3417|2368/-1|790/-1||,T|missense_variant|MODERATE|jg23969|jg23969|transcript|jg23969.t1|protein_coding|2/5|n.2683G>A|p.Glu895Lys|2683/3732|2683/-1|895/-1||
    my @annots = split(',', $annotations);
    ## we probably just want the first one, but let's check them all, who knows how many there are...
    foreach(@annots) {
        my $variant = $_ =~ s/\|/;/gr;
        my @variant = split(';', $variant);
        print $variant[2] . "\n";
        my $geneID = $variant[3];
        my $type = $variant[1];
        my $impact = $variant[2];
        my $n = $variant[9];
        my $p = $variant[10];
        print $type."\t".$impact."\t".$chrom."\t".$Pos."\t".$RefNT ."\t".$AltNT ."\t".$geneID ."\t".$n ."\t".$p."\n";
        print OUT $type."\t".$impact."\t".$chrom."\t".$Pos."\t".$RefNT ."\t".$AltNT ."\t".$geneID ."\t".$n ."\t".$p."\n";
    }
}



##
##

close(IN);
close(OUT);
