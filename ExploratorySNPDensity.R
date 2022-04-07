library(dplyr)
library(ggplot2)

getwd()


SNPs.df <- read.table("D:/Aedes_Illumina_Genomes/Analysis/Variants/ROCKv4_SNPs/RK_G1.redo.SNPdensity.RK_bigwindows-750-250kb.bed", sep="\t", header=FALSE)
colnames(SNPs.df) <- c("Chrom", "Start", "End", "Count")



SNPs.df <- mutate(SNPs.df, Mid = as.integer(Start + 125000))

### Were the ten extra windows just MTs? I didn't want those... 
SNPs.df <- filter(SNPs.df, Chrom != "MT_RagTag")
## No, that was only one extra one, obviously 
plot <- ggplot(data=SNPs.df) +
  geom_line(aes(x=Mid, y=Count), color = "darkolivegreen4") +
  scale_x_continuous(expand=c(0.01,0.01)) +
  facet_grid(rows = vars(Chrom)) + 
  theme_light() + 
  labs(title="SNP density across ROCK genome") + 
  xlab("Sliding windows (750kb, incremented by 250kb)") + 
  ylab("Count of SNPs in window")

plot

summary(SNPs.df$Count)
meanSNPs <- mean(SNPs.df$Count)
stdevSNPs <- sd(SNPs.df$Count)

setwd("Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Figures/")
pdf(file="RK_G1_SNPs_Density_Across_RK_Genome_BetterCovOutlierWindows_750kbWin.pdf", width=8.5)
plot
dev.off()

### Coverage plots: 
cov <- read.table("D:/Aedes_Illumina_Genomes/Analysis/Variants/ROCKv4_SNPs/ROCK_coverage_with_Chrom3End_250kb.bed", header=FALSE)
kbcov <- read.table("D:/Aedes_Illumina_Genomes/Analysis/Variants/ROCKv4_SNPs/ROCK_coverage_1kb.bed", sep="\t", header=F)
### kbcov is 140,000 windows, so we really can't plot it nicely across the genome. 
### But I do want to compare for myself kbcov vs 250kb cov. 
colnames(cov) <- c("Chr", "start", "end", "meanQ", "medianQ", "count")
colnames(kbcov) <- c("Chr", "start", "end", "meanQ", "medianQ", "count")
cov$meanQ <- as.numeric(cov$meanQ)
cov$medianQ <- as.numeric(cov$medianQ)
summary(na.omit(cov$medianQ))
summary(na.omit(cov$meanQ)) 
### These are the median and mean mapping quality, and we're not going to use that info. Drop it. 
cov <- cov %>% select(Chr, start, end, count)
kbcov <- kbcov %>% select(Chr, start, end, count)
colnames(cov) <- c("Chrom", "Start", "End", "Count")
colnames(kbcov) <- c("Chrom", "Start", "End", "Count")
table(cov$Chrom)
table(kbcov$Chrom) # oh no it includes all the non-chromosomal scaffolds 
kbcov <- kbcov %>% filter(Chrom == c("1_RagTag", "2_RagTag", "3_RagTag"))
table(kbcov$Chrom)
cov <- cov %>%
  mutate(coverage = Count*200/250000) ## I counted these in non-sliding windows -- 250000kb each 
summary(cov$coverage)                 ## to more accurately reflect coverage where SNPs were being counted

kbcov <- kbcov %>% 
  mutate(coverage = Count*200/1000)
summary(kbcov$coverage)

#### OK, they are about the same, and the mean is right around 25x. It's a little lower in the 250kb coverage, 
#### but it's also just got a smaller range. 


## Plotting SNPs overlaid on coverage 

plot <- ggplot(data=SNPs.df) +
  geom_area(data=cov, aes(x=Start+125000, y=coverage*20), alpha=0.5, fill="darkolivegreen4")+
  ### Multiplied coverage by 20 so that it's actually visible. 
  ### We dono't have a scale 
  geom_line(data=SNPs.df, aes(x=Mid, y=Count), size=0.05, color = "black") +
  scale_x_continuous(expand=c(0.01,0.01)) +
  facet_grid(rows = vars(Chrom)) + 
  theme(
    plot.title=element_text(size=48), 
    axis.title=element_text(size=36),
    axis.text=element_text(size=24), 
    strip.text.y=element_text(size=24)
  ) + 
  labs(title="SNP density across ROCK genome", 
       x="Sliding windows (750kb, incremented by 250kb)", 
       y="Count of SNPs in window") 

plot
getwd()
jpeg("SNPDensityAcrossRockGenome_OverlaidOnCoverage.jpg", width=2000, height=2000)
plot
dev.off()

png("SNPDensityAcrossRockGenome_OverlaidOnCoverage.png", width=2000, height=2000)
plot
dev.off()

table(cov$Chrom)
plot <- ggplot(data=cov) +
  geom_area(data=cov, aes(x=Start+125000, y=coverage), alpha=0.5, fill="darkolivegreen4")+
  #geom_line(data=SNPs.df, aes(x=Mid, y=Count), color = "black") +
  scale_x_continuous(expand=c(0.01,0.01)) +
  facet_grid(rows = vars(Chrom)) + 
  theme(
    plot.title=element_text(size=30), 
    axis.title=element_text(size=20),
    axis.text=element_text(size=20), 
    strip.text.y=element_text(size=20)
  ) +
  labs(title="Illumina read coverage in 250kb windows", 
       x="Genome position of window midpoint", 
       y="Coverage (read length * read count / 250000)") 
  
plot
png("JustCoverageTracks2.png", width=2000)
plot
dev.off()



SNPs.df[is.na(match(cov$Start, SNPs.df$Start)),] ## No missing blocks 

SNPs.df <- SNPs.df %>% mutate(win.id = paste(Chrom, "_", Start, sep=""))
cov <- cov %>% mutate(win.id = paste(Chrom, "_", Start, sep=""))

df <- full_join(SNPs.df, cov, by="win.id")
colnames(df)
df <- df %>% select("Chrom.x", "Start.x", "End.x", "Count.x", "win.id", "Count.y", "coverage")
colnames(df) <- c("Chrom", "Start", "End", "SNPs.Count", "win.id", "Reads.Count", "coverage")

## Normalizing SNPs count by coverage 
df <- df %>% mutate(SNPs_by_cov = SNPs.Count/coverage)

newplotdf <- df %>% select(Chrom, Start, End, SNPs.Count, Reads.Count, coverage, SNPs_by_cov)
colnames(newplotdf) <- c("Chrom", "Start", "End", "SNPsCount", "ReadCount", "coverage", "SNPs_by_cov")
summary(newplotdf$SNPs_by_cov)


newplotdf2 <- newplotdf[!is.na(newplotdf$SNPs_by_cov), ]

newplotdf2 <- newplotdf2 %>% mutate(log2SNPs = log2(SNPs_by_cov+.01))

plot <- ggplot(data=newplotdf2) +
  geom_line(data=newplotdf2, aes(x=End-125000, y=SNPs_by_cov), color="darkolivegreen4")+
  #geom_line(data=SNPs.df, aes(x=Mid, y=Count), color = "black") +
  scale_x_continuous(expand=c(0.01,0.01)) +
  scale_y_continuous(trans='log10') + 
  facet_grid(rows = vars(Chrom)) + 
  theme_minimal() + 
  theme(aspect.ratio=0.1) +
  labs(title="SNP density across ROCK genome", 
       caption="Dividing SNP counts by coverage reveals several islands of high and low polymorphism.
       Mean=56.10 (sd = 68.67), Median=44.43, Max=1791.33") + 
  xlab("Sliding windows (750kb, incremented by 250kb)") + 
  ylab("Count of SNPs per bases mapped") 

summary(newplotdf2$SNPs_by_cov)
sd(newplotdf2$SNPs_by_cov)
library(ggplot2)
plot <- ggplot(data=newplotdf2) +
  geom_line(data=newplotdf2, aes(x=End-125000, y=SNPs_by_cov), color="darkolivegreen4")+
  #geom_line(data=SNPs.df, aes(x=Mid, y=Count), color = "black") +
  scale_x_continuous(expand=c(0.01,0.01)) +
  scale_y_continuous(trans='log10') + 
  facet_grid(rows = vars(Chrom)) + 
  theme_minimal() + 
  theme(aspect.ratio=0.1) +
  labs(title="SNP density across ROCK genome", 
       caption="Dividing SNP counts by coverage reveals several islands of high and low polymorphism.
       Mean=57.56 (sd = 70.78), Median=45.37, Max=1891.45") + 
  xlab("Sliding windows (750kb, incremented by 250kb)") + 
  ylab("Count of SNPs per bases mapped") 

plot2 <- ggplot(data=newplotdf2) +
  geom_line(data=newplotdf2, aes(x=Start+125000, y=log2SNPs), color="darkolivegreen4")+
  #geom_line(data=SNPs.df, aes(x=Mid, y=Count), color = "black") +
  scale_x_continuous(expand=c(0.01,0.01)) +
  facet_grid(rows = vars(Chrom)) + 
  theme_minimal() + 
  theme(aspect.ratio=0.1) +
  labs(title="Depressed SNP density across ROCK genome", 
       caption="log2(SNPs-count / coverage depth) reveals islands of relatively depressed polymorphism.") + 
  xlab("Sliding windows (750kb, incremented by 250kb)") + 
  ylab("log2(SNPs per bases mapped)")



corr <- cor(as.numeric(newplotdf2$SNPsCount), 
            as.numeric(newplotdf2$ReadCount), 
            use="pairwise.complete.obs", method="pearson")
var(as.numeric(newplotdf2$SNPs_by_cov))
var(as.numeric(newplotdf2$SNPsCount))
sd(as.numeric(newplotdf2$SNPs_by_cov))
sd(as.numeric(newplotdf2$SNPsCount))

table(newplotdf$Chrom)
### this is good enough to show Jeff and ask him some questions 
setwd("Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Figures/")
pdf(file="RK_G1_SNPs_Density_Across_RK_Genome_with_coverage_750kbWinLog10Y.pdf", width=8.5)
plot
dev.off()
# cov3 <- cov %>% filter(Chrom=="3_RagTag")
# max(cov3$End)


summary(newplotdf2$SNPs_by_cov)
sd(newplotdf2$SNPs_by_cov)

highSNPs <- mean(newplotdf2$SNPs_by_cov) + (2*(sd(newplotdf2$SNPs_by_cov)))
### Two standard deviations above the mean = 193.5
lowSNPs <- mean(newplotdf2$SNPs_by_cov) - (2*(sd(newplotdf2$SNPs_by_cov)))
### That doesn't work nearly as well because the standard deviation is larger than the mean... 
### Let's look at bottom quartile?
summary(newplotdf2$SNPs_by_cov)
## 1st quartile is 20.61 

logten <- log10(newplotdf2$SNPs_by_cov)
table(logten<0.3)


table(newplotdf2$SNPs_by_cov>highSNPs)

### There are 74 windows whose normalized SNPdensity is > 2 standard deviations more than the mean
highestSNPs <- mean(newplotdf2$SNPs_by_cov) + (3*(sd(newplotdf2$SNPs_by_cov)))
table(newplotdf2$SNPs_by_cov>highestSNPs)
### There are 35 where it's > 3 sd. 

###Looking at the histogram... let's take a look at the extreme outlier windows -- whose normalized SNP density is 
### > 500 

table(newplotdf2$SNPs_by_cov>500)
## There's 15 of them. 

HighestDensityWindows <- newplotdf2 %>% filter(SNPs_by_cov>500)

table(HighestDensityWindows$Chrom)

## 2 windows of high polymorphism on Chr1, 8 on Chr2, and 5 on Chr3. 

write.table(HighestDensityWindows, "D:/Aedes_Illumina_Genomes/Analysis/Variants/ROCKv4_SNPs/HighestDensitySNPWindows.tab", row.names=F, quote=F, sep="\t" )


viola <- ggplot(data=newplotdf2, aes(x=Chrom, y=SNPsCount, group = Chrom)) + 
  geom_violin(stat="ydensity", position="dodge", aes(fill=Chrom, color=Chrom), alpha=0.5) + 
  geom_boxplot(stat="boxplot", position="dodge", size=0.2, width=0.05, fill="#00000000", color="black") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  labs(title="Distribution of polymorphisms across the ROCK genome", 
       subtitle="Number of SNPs counted per 750kb sliding (250kb) window (not normalized)", 
       x="Density and distribution of SNP counts by chromosome", 
       y="SNPs per window")

       # caption="This violin-wrapped boxplot shows the smoothed distribution 
       # of windowed SNP counts for each chromosome. Boxes show interquartile range, 
       # whiskers show 1.5x interquartile range, with outliers shown as dots.")
png("../Supplementals/ViolinPlot_SNPsCounts.png")
viola
dev.off()

viola2 <- ggplot(data=newplotdf2, aes(x=Chrom, y=SNPs_by_cov, group = Chrom)) + 
  geom_violin(stat="ydensity", position="dodge", aes(fill=Chrom, color=Chrom), 
              alpha=0.5) + 
  geom_boxplot(stat="boxplot", position="dodge", size=0.2, width=0.05, fill="#00000000", color="black") +
  scale_y_continuous(trans='log10') + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  labs(title="Coverage normalized SNP density on log scale", 
       subtitle="Normalized SNP density = SNP counts divided by read depth on a per window basis", 
       x="Density and distribution of SNP counts by chromosome", 
       y="Coverage normalized SNPs per window")
png("../Supplementals/ViolinPlot_LogTransformedSNPsperBaseMapped.png")
viola2
dev.off()

# SNPsStats <- summary(newplotdf2$SNPs_by_cov)
# SNPsMean=summary(newplotdf2$SNPs_by_cov)[4]
# SNPsIQR = as.numeric(SNPsStats[5] - SNPsStats[3])
# SNPs1Q = as.numeric(SNPsStats[3])
# SNPs3Q = as.numeric(SNPsStats[5])
# 
# SNPs1Q - (1.5 * SNPsIQR)



viola3 <- ggplot(data=newplotdf2, aes(x=Chrom, y=SNPs_by_cov, group = Chrom)) + 
  geom_violin(stat="ydensity", position="dodge", aes(fill=Chrom, color=Chrom), alpha=0.5) + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  geom_boxplot(stat="boxplot", position="dodge", size=0.2, width=0.05, fill="#00000000", color="black") +
  labs(title="Coverage normalized SNP density", 
       subtitle="Normalized SNP density = SNP counts divided by read depth on a per-window basis", 
       y="SNPs count normalized by coverage", 
       x="Density and distribution of SNPs counts by chromosome")
png("../Supplementals/ViolinPlot_SNPsPerBaseMapped.png")
viola3
dev.off()

pdf("Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Supplementals/SNPDensityViolinPlots.redo.pdf")
viola
viola3
viola2
dev.off()


library("MetBrewer") ## colors are out of control here, need to decide on something uniform
color_pals <- met.brewer(name="Cassatt2", n=3, type="discrete")
# New plot, add black line for highest density SNPs
plot <- ggplot(data=newplotdf2) +
  geom_line(data=newplotdf2, aes(x=End-125000, y=SNPs_by_cov), color=color_pals[3])+
  geom_hline(yintercept=500, color=color_pals[2]) + 
  geom_hline(yintercept=1, color=color_pals[2]) + 
  
  scale_x_continuous(expand=c(0.01,0.01)) +
  scale_y_continuous(trans='log10') + 
  facet_grid(rows = vars(Chrom)) + 
  theme_minimal() + 
  theme(aspect.ratio=0.1) +
  labs(title="SNP density across ROCK genome", 
       caption="Dividing SNP counts by coverage reveals several islands of high and low polymorphism.
           Mean=57.56 (sd = 70.78), Median=45.37, Max=1891.45") + 
  xlab("Sliding windows (750kb, incremented by 250kb)") + 
  ylab("Count of SNPs per bases mapped")

pdf("Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Figures/SNPdensity_by_cov.redo.pdf")
plot
dev.off()
table(newplotdf2$SNPs_by_cov<1)
### There are 107 windows where normalized SNP density is less than 1. 
table(newplotdf2$SNPs_by_cov==0)
### There are 7 where it is zero. I have a hard time believing those, and I wonder if they 
### are just windows where there just wasn't good read coverage, so SNPs didn't get counted at all

LowSNPDensityWindows <- newplotdf2 %>% filter(SNPs_by_cov != 0) %>%
  filter(SNPs_by_cov < 1)

write.table(LowSNPDensityWindows, "D:/Aedes_Illumina_Genomes/Analysis/Variants/ROCKv4_SNPs/LowestSNPDensityWindows.tab", row.names = F, quote = F)


pdf("Z:/Shared Documents/Cera Fisher/DraftPublications/ROCKGenome/Figures/SNP_density_final_plot.pdf")
plot
dev.off()

### gene density across genome 
genes.count <- read.delim("D:/Aedes_Illumina_Genomes/Analysis/Variants/ROCKv4_SNPs/ROCKv4.genes.redo.SNPdensity.RK_bigwindows-750-250kb.bed", sep="\t", header=F)
colnames(genes.count) <- c("Chrom", "Start", "End", "genes.count")

##cov <- cov %>% mutate(win.id = paste(Chrom, "_", Start, sep=""))
##df <- full_join(SNPs.df, cov, by="win.id")

genes.count <- genes.count %>% mutate(win.id = paste(Chrom, "_", Start, sep=""))
summary(genes.count$genes.count)
genes.count <- genes.count %>% filter(Chrom != "MT_RagTag")


viola4 <- ggplot(data=genes.count, aes(x=Chrom, y=genes.count)) +
   geom_violin(stat="ydensity", position="dodge", fill="#00000033") + 
  geom_boxplot(stat="boxplot", position="dodge", size=0.2, width=0.05, fill="#FFFFFF", color="black") + 
  theme_bw() + 
  labs(title="Gene density across the ROCK genome", 
       x="Chromosome", 
       y="Count of genes in 750kb windows (sliding by 250kb)") 


viola4

png("../Supplementals/Violin_GeneDensity.png")
viola4
dev.off()
HighGeneDensityWindows <- genes.count %>% filter(genes.count > 36) ##81!
HighGeneDensityWindows <- genes.count %>% filter(genes.count > 60) ##81!
LowSNPDensityWindows <- LowSNPDensityWindows %>% mutate(win.id = paste(Chrom, "_", Start, sep=""))
LowSNPDensityWindows

match(LowSNPDensityWindows$win.id, HighGeneDensityWindows$win.id)

newplotdf2 <- newplotdf2 %>% mutate(win.id = paste(Chrom, "_", Start, sep=""))

geneDensityVsSNPDensity <- full_join(SNPs.df, genes.count, by = c("Chrom", "Start"))

cor(x=geneDensityVsSNPDensity$genes.count, y=geneDensityVsSNPDensity$Count)
##slight negative correlation
par(cex=1)

plot(geneDensityVsSNPDensity$Count, geneDensityVsSNPDensity$genes.count)

scatter1 <- ggplot(data=geneDensityVsSNPDensity, aes(x=genes.count, y=Count)) + 
  geom_point(alpha=0.5) + theme_bw() + 
  labs(title="SNP density is uncorrelated to gene density", 
       x="Count of genes per 750kb window", 
       y="Count of SNPs per 750kb window")

scatter1

sd(genes.count$genes.count)
which(geneDensityVsSNPDensity$genes.count==53)
geneDensityVsSNPDensity[which(geneDensityVsSNPDensity$genes.count==118),]
geneDensityVsSNPDensity[which(geneDensityVsSNPDensity$genes.count==53),]

which(newplotdf2$End==373250000)
newplotdf2[which(newplotdf2$End==373250000),]
which(newplotdf2$End==73250000)
newplotdf2[which(geneDensityVsSNPDensity$genes.count>60),]

newplotdf3 <- full_join(newplotdf2, genes.count, by=c("Chrom", "Start"))

scatter2 <- ggplot(data=newplotdf3, aes(x=genes.count, y=SNPs_by_cov)) + 
  scale_y_continuous(trans='log10') +
  geom_point(alpha=0.5) + theme_bw() + 
  labs(title="SNP density is uncorrelated to gene density", 
       x="Count of genes per 750kb window", 
       y="Coverage normalized SNPs per window")

scatter2
png("../Supplementals/CovNormSNPs_by_gene_density_scatter.png")
scatter2
dev.off()


plot(newplotdf3$genes.count, newplotdf3$SNPs_by_cov)
cor(newplotdf3$SNPs_by_cov, newplotdf3$genes.count)

HighGenDensity <- newplotdf3[which(newplotdf3$genes.count>60),c(1,2,3,4,5,6,7,11)]
LowGeneDensity <- newplotdf3[which(newplotdf3$genes.count<5),c(1,2,3,4,5,6,7,11)]
