#!/usr/bin/env perl 
#
#tig00001274_pilon_RagTag        GFACS   gene    604210  604896  .       -       .       ID=jg40541
#tig00001274_pilon_RagTag        GFACS   CDS     604210  604896  .       -       .       Parent=jg40541.t1
#tig00001274_pilon_RagTag        GFACS   gene    654824  676787  .       +       .       ID=jg40548
#tig00001274_pilon_RagTag        GFACS   CDS     654824  656191  .       +       .       Parent=jg40548.t1
#tig00001274_pilon_RagTag        GFACS   CDS     676157  676555  .       +       .       Parent=jg40548.t1
#tig00001274_pilon_RagTag        GFACS   CDS     676611  676787  .       +       .       Parent=jg40548.t1
#tig00001274_pilon_RagTag        GFACS   intron  656192  676156  .       +       .       .
#tig00001274_pilon_RagTag        GFACS   intron  676556  676610  .       +       .       .
#
#
##### but we actually need it to look more like this ##### 
#1       havana  gene    150956201       150958296       .       +       .       gene_id "ENSMUSG00000102628"; gene_version "2"; gene_name "Gm37671"; gene_source "havana"; gene_biotype "TEC";
#1       havana  transcript      150956201       150958296       .       +       .       gene_id "ENSMUSG00000102628"; gene_version "2"; transcript_id "ENSMUST00000193198"; transcript_version "2"; gene_name "Gm37671"; gene_source "havana"; gene_biotype "TEC"; transcript_name "Gm37671-201"; transcript_source "havana"; transcript_biotype "TEC"; tag "basic"; transcript_support_level "NA (assigned to previous version 1)";
#1       havana  exon    150956201       150958296       .       +       .       gene_id "ENSMUSG00000102628"; gene_version "2"; transcript_id "ENSMUST00000193198"; transcript_version "2"; exon_number "1"; gene_name "Gm37671"; gene_source "havana"; gene_biotype "TEC"; transcript_name "Gm37671-201"; transcript_source "havana"; transcript_biotype "TEC"; exon_id "ENSMUSE00001342333"; exon_version "2"; tag "basic"; transcript_support_level "NA (assigned to previous version 1)";
#1       havana  gene    150983666       150984611       .       +       .       gene_id "ENSMUSG00000100595"; gene_version "2"; gene_name "Gm19087"; gene_source "havana"; gene_biotype "processed_pseudogene";
#1       havana  transcript      150983666       150984611       .       +       .       gene_id "ENSMUSG00000100595"; gene_version "2"; transcript_id "ENSMUST00000191430"; transcript_version "2"; gene_name "Gm19087"; gene_source "havana"; gene_biotype "processed_pseudogene"; transcript_name "Gm19087-201"; transcript_source "havana"; transcript_biotype "processed_pseudogene"; tag "basic"; transcript_support_level "NA (assigned to previous version 1)";
use warnings;
use strict;

my $infile = $ARGV[0]; 
my $outfile = $ARGV[1];
open (GTF, "<", $infile) or die "I cannot open ".$infile.", sorry!\n";

open (OUT, ">>", $outfile) or die "I cannot open ".$outfile.", sorry!\n";

print "OKAY, Here's how it's going to go down: \n";
print "I'm taking ".$infile." and trying to convert it to a nice embl format gtf, called ".$outfile.".\n";
print "Hold your horses...\n";

while (<GTF>) {
   my $line = $_;
   chomp $line;
   #   print "The line is ".$line."\n";
   my @Cols = split("\t",$line); 
   #print "Number of fields: ". scalar(@Cols) . "\n";
   my $Chrom = $Cols[0];
   my $source = "GFACS";
   my $Feat = $Cols[2]; 
   #print "Feature type written: ".$Feat."\n";
   $Feat =~ s/CDS/exon/g; 
   #print "Feature type is now: ".$Feat."\n";
   my $Start = $Cols[3]; # 4th column is start
   my $End = $Cols[4]; # end position 
   my $Score = "."; # we don"t have scores 
   my $Strand = $Cols[6]; #7th column is strand 
   my $Frame = $Cols[7]; #we mostly don"t have this, but keep whatever is there
   my $Attr = $Cols[8];
   #   print $Attr;
   ####### Now, how are we going to fix this stuff? 
   ### if the feature type eq intron, we don"t keep this line; go to next thing in the loop. 
   #
   ## if the feature type eq gene; 
   # then $Attr probably matches "ID=\(stuff we want\)"
   # reformat $Attr so that "ID=jg40533" becomes instead, "gene_id "jg40533""; 
   ### if instead the feature type eq exon:
   # $Attr probably looks like "Parent=jg40533.t1" 
   # We want to reformat that to look like:
   # "gene_id "jg40533"; transcript_id "jg4.0533.t1"";
   # then, we need to add a little more to help snpeff out: 
   # add to $Attr "gene_biotype "protein_coding"" and, 
   # if feature type eq exon, add to $Attr "; transcript_biotype "protein_coding";"
   if ($Feat eq "intron") {
	   #print "This Feature type is ".$Feat." so I am going to drop it \n";
	   next;
  } elsif ($Feat eq "gene") {
	  #print "This Feature type is ".$Feat." so I need to edit it. \n";
	  $Attr =~ /ID=(.*)$/; 
	  my $geneID = $1;
	  $geneID =~ s/\_g//;
	  $geneID =~ s/\..*$//;
	  #print "My gene id is: ".$geneID."\n";
	  #
	  my $newAttr = 'gene_id "' . $geneID . '"; gene_biotype "protein_coding";';
	  #print "New field will be: ".$newAttr."\n";
		my $newLine = $Chrom . "\t" .
			$source . "\t" .
			$Feat . "\t" .
			$Start . "\t" .
			$End . "\t" .
			$Score . "\t" .
			$Strand . "\t" .
			$Frame . "\t" .
			$newAttr . "\n";
			#print "My new line: \n";
			#print $newLine;
                 print OUT $newLine;
  } elsif ($Feat eq "exon") {
	  #print "This Feature type is ".$Feat." so I need to edit it, too, and add somet stuff. \n";
	 #tig00001274_pilon_RagTag        GFACS   CDS     676157  676555  .       +       .       Parent=jg40548.t1
	 $Attr =~ /Parent=(.*)$/;
	 ## Well, that doesn't work because some of our transcripts look like this: 
	 #Parent=90041_t
	 my $transcriptID = $1;
	         $transcriptID =~ /(.*)\..+$/; 
		 my $geneID = $1; 
	         $transcriptID =~ /(.*)\_.+$/;
		 $geneID = $1; 
		 #### Replace the transcript id with the gene id, because the fasta only has unique transcripts, 
		 #### and the IDs are the geneIDs. 
		 $transcriptID = $geneID;
	 #print "My parent gene id is: ".$geneID." and also, my transcript id is ".$transcriptID."\n";
	 my $newAttr = 'gene_id "' . $geneID . '"; gene_biotype "protein_coding"; transcript_id "' . $transcriptID . '"; transcript_biotype "protein_coding";';
	 #print "New field will be: ".$newAttr."\n";
		my $newLine = $Chrom . "\t" .
			$source . "\t" .
			$Feat . "\t" .
			$Start . "\t" .
			$End . "\t" .
			$Score . "\t" .
			$Strand . "\t" .
			$Frame . "\t" .
			$newAttr . "\n";
			#print "My new line: \n";
			#print $newLine;
		 print OUT $newLine;

  }
}


close GTF;
close OUT;

print "I do believe we're done! \n";

