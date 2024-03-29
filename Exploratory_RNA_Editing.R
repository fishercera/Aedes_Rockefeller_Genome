### Exploring RNA-Editing
# 
# set_wd <- function() {
#   library(rstudioapi)
#   current_path <- getActiveDocumentContext()$path 
#   setwd(dirname(current_path ))
#   print( getwd() )
# }
# set_wd()
# getwd()

library(dplyr)
library(tidyr)
setwd("D:/Aedes_Illumina_Genomes/Analysis/RNA_Editing/RNA_Editing_R/src")
###### After parsing out the Annotations from the Info field: 
setwd("../../snpEff_out/")

allSNPs.df <- read.delim("ROCK-RNA.diff-from-wgs.braker.snpeff.ChromsOnly.var_table.tab", sep="\t", 
                         header=F)

head(allSNPs.df)
colnames(allSNPs.df) <- c("Type", "Impact", "Chrom", "Pos", "Ref", "Alt", "GeneID", "ntAllele", "aaAllele")

allSNPs.df <- allSNPs.df %>% mutate(mutation = paste(Ref, Alt))
table(allSNPs.df$mutation)
table(allSNPs.df$Type)
allMutation.counts <- allSNPs.df %>% group_by(mutation) %>% count(mutation)
allMutation.counts
# removing the ones with more than one alternate 
allMutation.counts <- allMutation.counts %>% filter(n > 10)

library(ggplot2)

sumAll <- sum(allMutation.counts$n)
misMutation.counts <- allSNPs.df %>% filter(Type=="missense_variant") %>%
  group_by(mutation) %>% count(mutation)
misMutation.counts <- misMutation.counts %>% filter(n > 10)
sumMis <- sum(misMutation.counts$n)

allMutation.counts <- allMutation.counts %>% mutate(prop = n / sumAll)
head(allMutation.counts)

misMutation.counts <- misMutation.counts %>% mutate(prop = n/sumMis)
head(misMutation.counts)

allMutation.counts <- allMutation.counts %>% mutate(group="zall") ## dumb hack to keep my colors the same
misMutation.counts <- misMutation.counts %>% mutate(group="mis")

plotdf <- rbind(allMutation.counts, misMutation.counts)

plot <- ggplot(data=plotdf, aes(x=mutation, y=prop, group=group)) +
  geom_col(position="dodge", aes(fill=group, color=group)) + 
  labs(title="Proportion of substitutions in RNA reads",
       subtitle="missense variants compared to all variants",
       xlab="Substitution (Ref:Alt)", 
       ylab="Proportion of All Substitutions In Group") 



pdf("Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Figures/Substitutions-RNA-Seq_missense-vs-all.pdf")
plot
dev.off()

### we want to compare synonymous with non-synonymous 

synSNPs.df <- filter(allSNPs.df, Type=="synonymous_variant")
misSNPs.df <- filter(allSNPs.df, Type=="missense_variant")

synMutation.counts <- synSNPs.df %>% group_by(mutation) %>% count(mutation) %>% mutate(group="syn")
synMutation.counts <- synMutation.counts %>% filter(n>10)
misMutation.counts <- misSNPs.df %>% group_by(mutation) %>% count(mutation) %>% mutate(group="mis")
misMutation.counts <- misMutation.counts %>% filter(n>10)
### n>10 drops all the SNPs with two alts, because they're very low in number. 

sumSyn <- sum(synMutation.counts$n)
sumMis <- sum(misMutation.counts$n)

synMutation.counts <- synMutation.counts %>% mutate(prop=n/sumSyn)
misMutation.counts <- misMutation.counts %>% mutate(prop=n/sumMis)


plotdf <- rbind(synMutation.counts, misMutation.counts)
plot2 <- ggplot(data=plotdf, aes(x=mutation, y=prop, group=group)) + 
  geom_col(position="dodge", aes(fill=group, color=group)) + 
  theme_minimal()+
  labs(title="Frequency of substitutions at potential RNA-editing sites",
       subtitle="missense variants compared to synonymous variants",
       x="Substitution (Ref:Alt)", 
       y="Proportion of All Substitutions In Group") 
  
pdf("Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Figures/Substitutions-RNA-Seq.pdf")
plot2
dev.off()

### Get table of annotations for ROCK-RNA.diff-from-wgs missense. 
### Make a table of genes and annotations that might have RNA editing sites. 
setwd("D:/Aedes_Nanopore_Genomes/5_AnnotatedAssemblies/RK/ROCKv4/gFACs_output/")
annotations <- read.delim("ROCKv4.braker.comprehensive.entap_no_contam_lvl4.tsv", sep="\t", header=T)

#### Have to fix the gene ids so that they match the changes required to make snpEff work. What a mess!!! 
annotations$Query.Sequence <- sub("ID=", "", annotations$Query.Sequence) # No ID=
annotations$Query.Sequence <- sub(perl=TRUE, pattern="(.*)\\..*",        # No .1 expansion
                                  replacement="\\1", 
                                  x = annotations$Query.Sequence)
annotations$Query.Sequence <- sub(perl=TRUE, pattern="(.*)_.*",          # And no _g expansion
                                  replacement="\\1", 
                                  x = annotations$Query.Sequence)

# Only some of the fields from EnTAP are informative.
Informative <- annotations %>% select(Query.Sequence, Subject.Sequence, Description, Species, Seed.Ortholog, Predicted.Gene, 
                                 GO.Biological, GO.Cellular, GO.Molecular, KEGG.Terms, InterPro, 
                                 Protein.Database, Protein.Description)

write.table(Informative, row.names=F, quote=F, file="D:/Aedes_Nanopore_Genomes/5_AnnotatedAssemblies/RK/ROCKv4/Entap/RK_F3_v4.Comprehensive.Uniq.braker.out/RK_F3.v4.braker.annotations.usefulColsOnly.tsv")
write.csv(Informative, row.names=F, quote=F, file="D:/Aedes_Nanopore_Genomes/5_AnnotatedAssemblies/RK/ROCKv4/Entap/RK_F3_v4.Comprehensive.Uniq.braker.out/RK_F3.v4.braker.annotations.usefulColsOnly.csv")

RNA.edited <- Informative %>% filter(Query.Sequence %in% misSNPs.df$GeneID) # 2,512 genes
write.table(RNA.edited, "Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Supplementals/RNA_editing_missense_diff-from-wgs2.txt", row.names=F, quote=F, sep="\t")

RNA.edited2 <- Informative %>% filter(Query.Sequence %in% synSNPs.df$GeneID) # 4,773 genes with synonymous subs
write.table(RNA.edited, "Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Supplementals/RNA_editing_synonymous_diff-from-wgs.txt", row.names=F, quote=F, sep="\t")

colnames(misSNPs.df)
missenseAG <- misSNPs.df %>% filter(mutation == "A G")
missenseCT <- misSNPs.df %>% filter(mutation == "C T")

RNA.editedAG <- Informative %>% filter(Query.Sequence %in% missenseAG$GeneID) ## 474 genes
RNA.editedCT <- Informative %>% filter(Query.Sequence %in% missenseCT$GeneID) ## 638 genes
write.table(RNA.editedAG, file="Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Supplementals/RNA_editing_ag2.txt", sep="\t", row.names=F, quote=F)
write.table(RNA.editedCT,  sep="\t", row.names=F, quote=F, file="Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Supplementals/RNA_editing_ct2.txt")


#######################
# Rough, quick and dirty estimation of transition-to-transversion ratio

head(plotdf)
library(dplyr)
plotdf <- plotdf %>% mutate(type = case_when(
    mutation == "A G" ~ "Ts", 
    mutation == "C T" ~ "Ts", 
    mutation == "G A" ~ "Ts", 
    mutation == "T A" ~ "Ts", 
    TRUE ~ "Tv"
  ))

misTs <- plotdf %>% filter(group=="mis", type=="Ts")
countMisTs <- sum(misTs$n)
sum(misTs$prop)
misTv <- plotdf %>% filter(group=="mis", type=="Tv")
countMisTv <- sum(misTv$n)
sum(misTv$prop)
misTsTv <- countMisTs/countMisTv
misTsTv2 <- countMisTv/countMisTs
### missense raw Ts/Tv ratio is 1:1.04

synSNPs.df <- synSNPs.df %>% mutate(type = case_when(
  mutation == "A G" ~ "Ts", 
  mutation == "C T" ~ "Ts", 
  mutation == "G A" ~ "Ts", 
  mutation == "T A" ~ "Ts", 
  TRUE ~ "Tv"
))

### The synSNPs.df is still in its unaggregated form, so I can actually just use group_by and count -- 
synRatiosdf <- synSNPs.df %>% group_by(type) %>% count(type)
synRatiosdf$type
synRatiosdf
countSynTs <- synRatiosdf$n[1]
countSynTv <- synRatiosdf$n[2]
synTsTv2 <- countSynTv/countSynTs
### synonymous raw Ts/Tv ratio is 1:0.57 

sessionInfo()
library(dplyr)
library(ggplot2)
getwd()
setwd("../snpEff_out/")
RNAnonref50nonwgs <- read.delim("ROCK-RNA-nonref50-nonwgs-annotated.tab", sep="\t", header=F)
colnames(RNAnonref50nonwgs) <- c("Type", "Impact", "Chrom", "Pos", "Ref", "Alt", "GeneID", 
                                 "nucl", "prot", "rank")
## "rank" here refers to whether this is the first, second, third etc. annotation at this position. 

RNAnonref50nonwgsnonrepet <- RNAnonref50nonwgs %>% filter(rank==1)
## leaves us with 15,453 instead of 23,318 

RNA3 <- RNAnonref50nonwgsnonrepet

RNA3 <- RNA3 %>% mutate(mutation=paste(Ref, Alt))
table(RNA3$mutation)

RNA3 <- RNA3 %>% mutate(type = case_when(
  mutation == "A G" ~ "Ts", 
  mutation == "C T" ~ "Ts", 
  mutation == "G A" ~ "Ts", 
  mutation == "T A" ~ "Ts", 
  TRUE ~ "Tv"
))

RNA3Mutation.counts <- RNA3 %>% group_by(mutation) %>% count(mutation)
RNA3Mutation.counts <- RNA3Mutation.counts %>% filter(n>10)
RNA3Mutation.counts


RNA3.mis <- RNA3 %>% filter(Type == "missense_variant")
## 1,887
table(RNA3.mis$mutation)



###### Let's look at the following case:
###### Variant is present in ROCK-RNA reads at 0.75 non-ref allele frequency. 
###### Variant is not present in RK_G1 (whole-genome sequencing data)
###### Variant does not fall in an area of "coverage outliers" for RK_G1 
######    (which means if it were present in the WGS pool-seq, it would 
######     have been detected.)
getwd()
RNA75 <- read.table("../snpEff_out/ROCK-RNA.nonref75.snpEff_annotated.tab", sep="\t", 
                    header=TRUE)

### There are 312 variants annotated here -- but there were not that many SNPs, 
### some of them have been double-annotated. 

RNA75snps <- read.table("../ROCK-RNA-nonref75.notinwgs.notInCovOutliers.bed", sep="\t", 
                        header=F)
### There are only 292 SNPs. 
### Let's be smart about how we deal with the "double-counting" problem and 
### take a look at the ones that have multiple annotations. 

max(RNA75$Rank)
## OK, none of them have more than 2 annotations, that's good 

RNA75_doubles <- filter(RNA75, Rank == 2)
## there's only 20 of these -- filtering on Rank == 2 only gives me their second annotation though 

## Can I safely filter on Position? Or do I actually have duplicate Positions on different chromosomes? 
RNA75 %>% group_by(Position) %>% count(Position) %>% filter(n>2)
## nah we're good 

RNA75_troubles <- RNA75 %>% filter(Position %in% RNA75_doubles$Position)
## as expected, that dataframe has 40 entries -- one for each of our 20 double-counted variants 
View(RNA75_troubles)
write.table(RNA75_troubles, file="out/DoubleCounted_RNASNPs.txt", sep="\t", quote=F, row.names = F)
## The results aren't so bad. It looks like the double-counting is resulting 
## almost entirely from duplicative gene models from Braker. Alas, that can't really 
## be helped right now. 

