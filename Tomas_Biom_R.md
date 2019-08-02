---
title: "R initial analysis"
output: 
  html_document: 
    keep_md: yes
---



# Import data

These are the libraries used here:


```r
library(phyloseq)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(tidyr)
library(tibble)
library(ggplot2)
library(vegan)
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.5-5
```

```r
library(readr)
```

Importing biom file:


```r
#Import file
TomAR <- import_biom("otu_table_TomasAR.clean.fungi.biom", refseqfilename = "otu_table_TomasAR.clean.fungi.fasta")
```

```
## Warning in strsplit(conditionMessage(e), "\n"): input string 1 is invalid
## in this locale
```

```r
TomAR
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 723 taxa and 81 samples ]
## sample_data() Sample Data:       [ 81 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 723 taxa by 7 taxonomic ranks ]
## refseq()      DNAStringSet:      [ 723 reference sequences ]
```

```r
#Renaming taxonomy levels on tax_table
colnames(tax_table(TomAR)) <- c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")

head(tax_table(TomAR))
```

```
## Taxonomy Table:     [6 taxa by 7 taxonomic ranks]:
##        Kingdom    Phylum             Class               
## OTU_7  "k__Fungi" "p__Ascomycota"    "c__Dothideomycetes"
## OTU_10 "k__Fungi" "p__Ascomycota"    "c__Sordariomycetes"
## OTU_12 "k__Fungi" "p__Ascomycota"    "c__Sordariomycetes"
## OTU_14 "k__Fungi" "p__Ascomycota"    "c__Sordariomycetes"
## OTU_22 "k__Fungi" "p__Ascomycota"    "c__Sordariomycetes"
## OTU_24 "k__Fungi" "p__Basidiomycota" "c__Agaricomycetes" 
##        Order              Family                    Genus             
## OTU_7  "o__Pleosporales"  "f__Didymellaceae"        "g__Endophoma"    
## OTU_10 "o__Hypocreales"   "f__unidentified"         "g__unidentified" 
## OTU_12 "o__Glomerellales" "f__Plectosphaerellaceae" "g__Gibellulopsis"
## OTU_14 "o__Glomerellales" "f__Plectosphaerellaceae" "g__Verticillium" 
## OTU_22 "o__Hypocreales"   "f__unidentified"         "g__unidentified" 
## OTU_24 "o__Boletales"     "f__Suillaceae"           "g__Suillus"      
##        Species                    
## OTU_7  "s__Endophoma_elongata"    
## OTU_10 "s__unidentified"          
## OTU_12 "s__unidentified"          
## OTU_14 "s__Verticillium_dahliae"  
## OTU_22 "s__unidentified"          
## OTU_24 "s__Suillus_pseudobrevipes"
```

# Checking sequencing depth and normalizing


```r
# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(TomAR))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred") +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](Tomas_Biom_R_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

Standardizing by sequencing depth:


```r
#Standardize abundances to the median sequencing depth
(total <- median(sample_sums(TomAR)))
```

```
## [1] 1909
```

```r
standf <- function(x, t=total) round(t * (x/sum(x)))
TomAR.std <- transform_sample_counts(TomAR, standf)
```

# Taxa abudance at different rank levels


```r
#Sample
head(sample_data(TomAR.std))
```

```
##             BarcodeSequence          Description        InputFileName
## AR.P001.Tom                 AR.P001.Tom-R1.fastq AR.P001.Tom-R1.fastq
## AR.P002.Tom                 AR.P002.Tom-R1.fastq AR.P002.Tom-R1.fastq
## AR.P003.Tom                 AR.P003.Tom-R1.fastq AR.P003.Tom-R1.fastq
## AR.P004.Tom                 AR.P004.Tom-R1.fastq AR.P004.Tom-R1.fastq
## AR.P005.Tom                 AR.P005.Tom-R1.fastq AR.P005.Tom-R1.fastq
## AR.P006.Tom                 AR.P006.Tom-R1.fastq AR.P006.Tom-R1.fastq
##             LinkerPrimerSequence ReadType Samples
## AR.P001.Tom                            R1 AR.P001
## AR.P002.Tom                            R1 AR.P002
## AR.P003.Tom                            R1 AR.P003
## AR.P004.Tom                            R1 AR.P004
## AR.P005.Tom                            R1 AR.P005
## AR.P006.Tom                            R1 AR.P006
```

```r
#summarizing by tax rank
Tom.class <- tax_glom(TomAR.std, "Phylum")
plot_bar(Tom.class, fill = "Phylum", x = "Samples") + coord_flip() + theme_gray() 
```

![](Tomas_Biom_R_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

## Exporting to tsv tabel

```r
TomAR.tsv <- psmelt(TomAR.std)

#Writting tabulated OTU table
#write_csv2(TomAR.tsv, "TomAR.tsv")
```


### Saving R object

```r
#save(TomAR, TomAR.std, TomAR.tsv, file = "TomAR_files.rda")
```
