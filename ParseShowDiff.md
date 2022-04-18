Parsing ShowDiffs
================

Take the output from MUMmer4 Show-Diff and parse the nature of the
breakpoints according to MUMmer’s manual of how to interpret what the
breakpoints are.

``` r
dff <- read.delim("D:/Aedes_Nanopore_Genomes/6_Analysis/GenomeAlignment/RK_on_LVP/RKv4-on-LVPChroms-Show-Diff-pid85.l10000.REF.txt", sep="\t", header=TRUE)
colnames(dff) <- c("Chromosome", "Type", "Ref.Start", "Ref.End", "Ref.Length", "Extra6", "Extra7")
```

We need dplyr to break these out into different tables, because the
columns hold different kinds of information based on the type of
breakpoint.

``` r
library(dplyr)

dff.GAP <- filter(dff, Type == "GAP")
# columns names: IDR GAP gap-start gap-end gap-length-R gap-length-Q gap-diff
dff.BRK <- filter(dff, Type == "BRK")
# columns  names: IDR BRK gap-start gap-end gap-length
dff.JMP <- filter(dff, Type=="JMP")
# columns names: IDR JMP gap-start gap-end gap-length
dff.INV <- filter(dff, Type=="INV")
# columns names: IDR INV gap-start gap-end gap-length
dff.SEQ <- filter(dff, Type=="SEQ")
# columns names: IDR SEQ gap-start gap-end gap-length prev-sequence next-sequence
dff.DUP <- filter(dff, Type=="DUP")
# columns names: IDR DUP dup-start dup-end dup-length
print(paste("Are there any DUP? No: length(dff.DUP$Chromosome) = ", length(dff.DUP$Chromosome), sep=" "))
```

    ## [1] "Are there any DUP? No: length(dff.DUP$Chromosome) =  0"

``` r
table(dff$Type)
```

    ## 
    ##   BRK   GAP   INV   JMP   SEQ 
    ##     6 13220   127   573   568

##### Get these numbers broken down by chromosome

``` r
# "Other Breakpoint"
table(dff.BRK$Chromosome)  
```

    ## 
    ## 1 2 3 
    ## 2 2 2

``` r
# "Inversion with possible relocation"
table(dff.INV$Chromosome)
```

    ## 
    ##  1  2  3 
    ## 26 43 58

``` r
# "Relocation - same chromosome"
table(dff.JMP$Chromosome)
```

    ## 
    ##   1   2   3 
    ## 108 254 211

``` r
# "Relocation - different chromosome"
table(dff.SEQ$Chromosome)
```

    ## 
    ##   1   2   3 
    ## 118 246 204

``` r
# "Alignment gaps - all"
table(dff.GAP$Chromosome)
```

    ## 
    ##    1    2    3 
    ## 2895 5541 4784

This is the data that goes into **Table Breakpoints.** We’ll just need
to do a little more parsing of the “GAP” data to understand whether the
gap represents a deletion or an insertion in the ROCK genome relative to
LVP.

### Parsing Gaps

Adjust column names

``` r
colnames(dff.GAP)
```

    ## [1] "Chromosome" "Type"       "Ref.Start"  "Ref.End"    "Ref.Length"
    ## [6] "Extra6"     "Extra7"

``` r
# columns names: IDR GAP gap-start gap-end gap-length-R gap-length-Q gap-diff
colnames(dff.GAP)[c(6,7)] <- c("Qry.Length", "gap.diff")
```

From the mummer manual: <https://mummer4.github.io/manual/manual.html>

> \[GAP\] A gap between two mutually consistent ordered and oriented
> alignments. gap-length-R is the length of the alignment gap in the
> reference, gap-length-Q is the length of the alignment gap in the
> query, and gap-diff is the difference between the two gap lengths. If
> gap-diff is positive, sequence has been inserted in the reference. If
> gap-diff is negative, sequence has been deleted from the reference. If
> both gap-length-R and gap-length-Q are negative, the indel is tandem
> duplication copy difference.

``` r
## 1) Is there an insertion in the reference?
# "If gap-diff is positive, sequence has been inserted in the reference."
dff.GAP <- dff.GAP %>%
  mutate(RefInsert = 
           case_when(
             gap.diff > 0 ~ "TRUE", 
             gap.diff < 0 ~ "FALSE"
           ))

## 2) Is the indel tandem duplication copy difference? 
#"If both gap-length-R and gap-length-Q are negative, the indel is
# tandem duplication copy difference."
dff.GAP <- dff.GAP %>%
  mutate(CopyNumDiff = 
           case_when(
             Ref.Length <0 & Qry.Length <0 ~ "TRUE"
           ))


## 3) It follows that if CopyNumDiff == True and RefInsert == True, then the 
## copy number is higher in the reference than the query; if RefInsert == FALSE, 
## the query has the higher count of copy numbers
dff.GAP <- dff.GAP %>%
  mutate(RefCopyDups = 
           case_when(
             RefInsert == TRUE & CopyNumDiff == TRUE ~ "TRUE", 
             RefInsert == FALSE & CopyNumDiff == TRUE ~ "FALSE"
           ))

# Now let's count them:
table(dff.GAP$RefInsert)
```

    ## 
    ## FALSE  TRUE 
    ##  6920  6291

``` r
table(dff.GAP[which(dff.GAP$Chromosome==1),]$RefInsert)
```

    ## 
    ## FALSE  TRUE 
    ##  1517  1375

``` r
table(dff.GAP[which(dff.GAP$Chromosome==2),]$RefInsert)
```

    ## 
    ## FALSE  TRUE 
    ##  2910  2628

``` r
table(dff.GAP[which(dff.GAP$Chromosome==3),]$RefInsert)
```

    ## 
    ## FALSE  TRUE 
    ##  2493  2288

``` r
write.table(dff.GAP, file="RKv4-on-LVP.showdiffs.GAPS-in-reference.txt", row.names=FALSE, quote=FALSE)
```

#### Store these for later inspection

These are all potentially interesting for further scrutiny, but to make
sense of them later, I want to store them separately with the correct
field names (tidy data).

“BRK” - reported as “Other Breakpoint”

> \[BRK\] An insertion in the reference of unknown origin, that
> indicates no query sequence aligns to the sequence bounded by
> gap-start and gap-end. Often found around DUP elements or at the
> beginning or end of sequences.

``` r
colnames(dff.BRK) <- c("Chromosome", "Type", "Ref.Gap.Start", "Ref.Gap.End", "Ref.Gap.Length", "", "")
dff.BRK <- dff.BRK %>% select(Chromosome, Type, Ref.Gap.Start, Ref.Gap.End, Ref.Gap.Length)

write.table(dff.BRK, file="RKv4-on-LVP.showdiffs.BRKs-in-reference.txt", row.names=FALSE, quote=FALSE)
```

“JMP” - reported as “Relocation (Same Chromosome)”

> \[JMP\] A relocation event, where the consistent ordering of
> alignments is disrupted. The coordinate columns specify the
> breakpoints of the relocation in the reference, and the gap-length
> between them. A negative gap-length indicates the relocation occurred
> around a repetitive sequence, and a positive length indicates unique
> sequence between the alignments.

``` r
colnames(dff.JMP)
```

    ## [1] "Chromosome" "Type"       "Ref.Start"  "Ref.End"    "Ref.Length"
    ## [6] "Extra6"     "Extra7"

``` r
# columns names: IDR JMP gap-start gap-end gap-length
dff.JMP <- dff.JMP %>% select(Chromosome, Type, Ref.Start, Ref.End, Ref.Length)
dff.JMP <- dff.JMP %>% 
  mutate(RepetitiveRegion = 
           case_when(
             Ref.Length < 0 ~ "TRUE", 
             Ref.Length > 0 ~ "FALSE"
           ))
table(dff.JMP$RepetitiveRegion)
```

    ## 
    ## FALSE  TRUE 
    ##   435   137

``` r
write.table(dff.JMP, file="RKv4-on-LVP.showdiffs.JMPs-in-reference.txt", row.names=F, quote=F)
```

“SEQ” - reported as “Relocation (Different Chromosome)”

> \[SEQ\] A translocation event that requires jumping to a new query
> sequence in order to continue aligning to the reference. If each input
> sequence is a chromosome, these features correspond to
> inter-chromosomal translocations.

``` r
colnames(dff.SEQ)
```

    ## [1] "Chromosome" "Type"       "Ref.Start"  "Ref.End"    "Ref.Length"
    ## [6] "Extra6"     "Extra7"

``` r
# columns names: IDR SEQ gap-start gap-end gap-length prev-sequence next-sequence
colnames(dff.SEQ) <- c("Chromosome", "Type", "Ref.Start", "Ref.End", "Ref.Length", "Query.Prev.Seq", "Query.Next.Seq")

write.table(dff.SEQ, file="RKv4-on-LVP.showdiffs.SEQs-in-reference.txt", row.names=F, quote=F)
```

“INV” - reported as “Inversion (with possible relocation)”

> \[INV\] The same as a relocation event, however both the ordering and
> orientation of the alignments is disrupted. Note that for JMP and INV,
> generally two features will be output, one for the beginning of the
> inverted region, and another for the end of the inverted region.

I think that the note about the number of JMP and INV representing
double the number of actual events is probably also true for SEQ, since
the table reports a breakpoint when the alignment moves to another
scaffold and when it moves back to the original scaffold.

``` r
colnames(dff.INV)
```

    ## [1] "Chromosome" "Type"       "Ref.Start"  "Ref.End"    "Ref.Length"
    ## [6] "Extra6"     "Extra7"

``` r
# columns names: IDR INV gap-start gap-end gap-length
dff.INV <- dff.INV %>% select(Chromosome, Type, Ref.Start, Ref.End, Ref.Length)

write.table(dff.INV, file="RKv4-on-LVP.showdiffs.INVs-in-reference.txt", row.names = F, quote=F)
```

Finally, my session info:

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18363)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.0.8
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.13  knitr_1.37       magrittr_2.0.1   tidyselect_1.1.1
    ##  [5] R6_2.5.1         rlang_1.0.2      fastmap_1.1.0    fansi_1.0.0     
    ##  [9] stringr_1.4.0    tools_4.1.1      xfun_0.27        utf8_1.2.2      
    ## [13] DBI_1.1.2        cli_3.1.0        htmltools_0.5.2  ellipsis_0.3.2  
    ## [17] assertthat_0.2.1 yaml_2.2.1       digest_0.6.28    tibble_3.1.6    
    ## [21] lifecycle_1.0.1  crayon_1.4.2     purrr_0.3.4      vctrs_0.3.8     
    ## [25] glue_1.6.0       evaluate_0.14    rmarkdown_2.11   stringi_1.7.6   
    ## [29] compiler_4.1.1   pillar_1.6.4     generics_0.1.1   pkgconfig_2.0.3
