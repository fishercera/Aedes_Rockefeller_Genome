### circlize manual https://jokergoo.github.io/circlize_book/book/

### Basic circular Aedes aegypti genome plot for SNP density
set_wd <- function() {
  library(rstudioapi)
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
  print( getwd() )
}
set_wd()

library(circlize)
library(dplyr)
library(tidyr)
#### Chunk 1: Load and modify data ###


cytoband.file <- "D:/R_Scripts/Aedes_ROCK_cytobands_saved2.txt"

circos.clear() # I always call this before I start a Circos plot. 

circos.par(start.degree=90, gap.degree=10) # Some aesthetic parameters
circos.initializeWithIdeogram(cytoband.file) # makes the outer track

#### tutorial

set.seed(999)
mat = matrix(sample(18, 18), 3, 6) 
rownames(mat) = paste0("S", 1:3)
colnames(mat) = paste0("E", 1:6)
mat

df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df

chordDiagram(mat)
mat
#######

dupes <- read.delim("ROCKv4_DuplicatedBuscos_Sorted.bed", header=F, sep="\t")

colnames(dupes) <- c("Chromosome", "start", "end", "name")

grouped_dupes <- dupes %>% group_by(name) %>% count(name)



top_dupes <- filter(dupes, name=="Name=17610at50557.faa")
top_dupes

#######
circos.clear()
set.seed(123)
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 20), ]
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 20), ]

circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
                   border = NA)

bed1
bed2

bed1 <- bed1[1,]
bed2 <- bed2[1,]

circos.clear()
circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
                   border = NA)

bed1
bed2
circos.clear()
circos.initializeWithIdeogram(cytoband.file)

bed1
bed1$chr[1] <- "1"
bed2
bed2$chr[1] <- "3"

circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
                   border = NA)

top_dupes

circos.genomicLink(top_dupes[1,], top_dupes[4,])

bed1



circos.clear()
circos.initializeWithIdeogram(cytoband.file)

circos.clear()
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[c(1:3), ]
bed1
bed_df <- read.delim(cytoband.file, sep="\t", header=F)

str(bed_df)
colnames(bed_df) <- c("Chrom", "start", "stop", "name", "value")
bed_df$Chrom <- as.character(bed_df$Chrom)
str(bed_df)
bed_df$Chrom <- as.factor(bed_df$Chrom)

circos.clear()
circos.par(start.degree=90, gap.degree=10)
circos.genomicInitialize(bed_df)

df_1 <- bed_df[3,]
df_2 <- bed_df[40,]
df_3 <- bed_df[4, ]
df_4 <- bed_df[6, ]
df_2

circos.genomicLink(df_1, df_2)
circos.genomicLink(df_3, df_4, col=rand_color(1, transparency = 0.4))

chr1 <- filter(bed_df, Chrom=="1")
str(chr1)
chr1$Chrom <- as.factor(as.character("1"))
str(chr1)


circos.clear()
circos.genomicInitialize(chr1)

########

head(dupes)

dupes <- dupes %>% group_by(name)

str(dupes)

dupes_chr1 <- filter(dupes, Chromosome=="1")
dupes_chr1 <- ungroup(dupes_chr1)

head(dupes_chr1)

bed1 <- ""
bed2 <- ""
bed3 <- ""
bed4 <- ""

############### have to parse this out with perl, but let's take a look at it first 

dupes <- read.delim("ROCK.v4.DuplicatedBuscos_diptera.sorted.bed", header=F, sep="\t")
## 529 duplicates from the Diptera SCO set 
colnames(dupes) <- c("Chrom", "Start", "End", "Name")
dupes %>% group_by(Name) %>% count(Name) %>% arrange(desc(n))

top_dupes <- dupes %>% filter()

bed1 <- read.delim("ROCKv4.DipteraDupes.bed1", header=F, sep="\t")

bed2 <- read.delim("ROCKv4.DipteraDupes.bed2", header=F, sep="\t")

rm(dupes)

dupes <- full_join(bed1, bed2, by="V4")
colnames(dupes)
colnames(dupes) <- c("Chrom_bed1", "Start_bed1", "End_bed1", "Name", "Chrom_bed2", "Start_bed2", "End_bed2")

small_test <- dupes[c(1:5), ]

circos_genomic_bed <- read.delim("ROCKv4_Fasta_ForDupes.bed", sep="\t", header=F)

colnames(circos_genomic_bed) <- c("Chrom", "Start", "End")
circos_genomic_bed$Chrom <- as.factor(circos_genomic_bed$Chrom)
str(circos_genomic_bed)
circos.clear()
circos.par(start.degree=90, gap.degree=10)
circos.genomicInitialize(circos_genomic_bed)
circos.genomicLink(small_test[,c(1:3)], small_test[,c(5:7)], 
                   col = rand_color(nrow(small_test), transparency = 0.4))

circos.clear()
circos.par(start.degree=90, gap.degree=10)
circos.genomicInitialize(circos_genomic_bed)
circos.genomicLink(dupes[,c(1:3)], dupes[,c(5:7)], 
                   col = rand_color(nrow(dupes), transparency = 0.4))

pdf("Diptera_Duplicated_Buscos-all.pdf")
circos.clear()
circos.par(start.degree=90, gap.degree=10)
circos.genomicInitialize(circos_genomic_bed)
circos.genomicLink(dupes[,c(1:3)], dupes[,c(5:7)], 
                   col = rand_color(nrow(dupes), transparency = 0.4))
dev.off()

#######

chromsbed2 <- read.delim("ROCKv4.DipteraDupes.chromsonly.bed2", sep="\t", header=F)
colnames(chromsbed2) <- c("Chrom", "Start", "End", "Name")
colnames(bed1) <- c("Chrom", "Start", "End", "Name")
chroms_dupes <- right_join(bed1, chromsbed2, by='Name')

chroms_only_bed <- circos_genomic_bed[c(1:3),]
chroms_only_bed$Chrom <- as.factor(as.character(chroms_only_bed$Chrom))
circos.clear()
circos.par(start.degree=90, gap.degree=10)
circos.genomicInitialize(chroms_only_bed)
circos.genomicLink(chroms_dupes[,c(1:3)], chroms_dupes[,c(5:7)], 
                   col = rand_color(nrow(chroms_dupes), transparency = 0.4))


pdf("ROCKv4_DipteraDupes_ChromosomesOnly.pdf")
circos.clear()
circos.par(start.degree=90, gap.degree=10)
circos.genomicInitialize(chroms_only_bed)
circos.genomicLink(chroms_dupes[,c(1:3)], chroms_dupes[,c(5:7)], 
                   col = rand_color(nrow(chroms_dupes), transparency = 0.4))
dev.off()


##### Let's zoom in on a potential problem area 

zoom_bed <- data.frame("1", 100556000, 119190000)
colnames(zoom_bed) <- c("Chrom", "Start", "End")
zoom_bed$Chrom <- as.factor(zoom_bed$Chrom)
str(zoom_bed)

zoom_dupes <- dupes_chroms_only %>% filter(Chrom_bed1=="1", Start_bed1 > 80000000, End_bed1 < 160000000) %>%
  filter(Chrom_bed2=="1", Start_bed2 > 80000000, End_bed2 < 160000000) %>%
  filter(Start_bed1 < 131445158)

circos.clear()
circos.par(start.degree=360, gap.degree=180)
circos.genomicInitialize(zoom_bed)
circos.genomicLink(zoom_dupes[,c(1:3)], zoom_dupes[,c(5:7)], 
                   col=rand_color(nrow(zoom_dupes), transparency = 0.2))

min(zoom_dupes$Start_bed1)
max(zoom_dupes$End_bed2)

pdf("ROCKDipteraDupes_ZoomChrom1.pdf")
circos.clear()
circos.par(start.degree=360, gap.degree=180, xaxis.clock.wise=TRUE)
circos.genomicInitialize(zoom_bed)
circos.genomicLink(zoom_dupes[,c(5:7)], zoom_dupes[,c(1:3)], 
                   col=rand_color(nrow(zoom_dupes), transparency = 0.2), inverse=TRUE)
dev.off()

zoom_bed$End-zoom_bed$Start


zoom_bed1 <- zoom_dupes[,c(1:4)]
zoom_bed2 <- zoom_dupes[,c(5:7,4)]

getwd()
write.table(zoom_bed1, "Chr1Dupes_set1.bed", eol="\n", sep="\t", quote=F, row.names=F)
write.table(zoom_bed2, "Chr1Dupes_set2.bed", eol="\n", sep="\t", quote=F, row.names=F)
