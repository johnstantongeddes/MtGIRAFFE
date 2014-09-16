Identification of haplotype clusters for GIRAFFE gene in the Medicago HapMap
================================================================================

**Author:** John Stanton-Geddes
**Date:** 16 September 2014

## Summary

The goal of this bioinformatic analysis is to determine the number of haplotype clusters that exist for a heme oxygenase gene, Medtr8g019302.1, that was annotated in the Medicago truncatula 4.0 genome. This gene, named GIRAFFE, has been studied in the Harris lab as it has effects on plant development and nodulation[^1].

[^1] http://www.uvm.edu/~plantbio/harris.screen.php

## Data

First, I downloaded the file containing Mt4.0v1 annotation SNPs for chromosome 8 from the GWAS accessions (262) from the Medicago HapMap website (http://www.medicagohapmap.org/downloads/mt40). The SNP file is in the binary variant call format (BCF) described [here](http://www.1000genomes.org/wiki/analysis/variant-call-format/bcf-binary-vcf-version-2). I also downloaded the `bcf.csi` file which is needed for `bcftools` filtering of the file.


```r
library(stringr)
library(plyr)

chr <- "8"
site.range <- "6779092-6783193"
datadir <- "data/"

file <- paste(datadir, "chr", chr, "-filtered-set-2014Apr15.bcf", sep="")

if(!file.exists(file)) {
    system('wget http://www.medicagohapmap.org/downloads/Mt40/snps_by_chr//chr8-filtered-set-2014Apr15.bcf')
    system('wget http://www.medicagohapmap.org/downloads/Mt40/snps_by_chr//chr8-filtered-set-2014Apr15.bcf.csi')
} else print("File already downloaded!")
```

```
## [1] "File already downloaded!"
```

Next, I used [bcftools](http://samtools.github.io/bcftools/) to extract SNPs from the region containing the GIRAFFE gene. This requires that bcftools are installed and available system-wide.


```r
txt.out <- paste(datadir, "chr", chr, "_", site.range, ".txt", sep="")
#system(paste('bcftools query -H -r chr', chr,':' , site.range," -f '%CHROM\t%POS[\t%TGT]\n' ", file, " > ", txt.out, sep=""))
```

Examination of this file revealed a number of SNPs that are heterozygous within accessions.


```r
haplo.data <- read.table(txt.out)
t.haplo.data <- t(haplo.data[,2:ncol(haplo.data)])
colnames(t.haplo.data) <- paste("snp", t.haplo.data[1,], sep="")
t.haplo.data <- t.haplo.data[-1,]

# identify sites heterozygous within accessions
t.haplo.data[which(t.haplo.data == "A/C")]
```

```
## [1] "A/C" "A/C" "A/C" "A/C"
```

```r
t.haplo.data[which(t.haplo.data == "A/G")]
```

```
## [1] "A/G" "A/G" "A/G"
```

```r
t.haplo.data[which(t.haplo.data == "A/T")]
```

```
## [1] "A/T" "A/T"
```

```r
t.haplo.data[which(t.haplo.data == "C/A")]
```

```
## [1] "C/A" "C/A"
```

```r
t.haplo.data[which(t.haplo.data == "C/G")]
```

```
## [1] "C/G" "C/G" "C/G" "C/G" "C/G"
```

```r
t.haplo.data[which(t.haplo.data == "C/T")]
```

```
## [1] "C/T" "C/T"
```

```r
t.haplo.data[which(t.haplo.data == "T/A")]
```

```
## [1] "T/A" "T/A"
```

```r
t.haplo.data[which(t.haplo.data == "T/C")]
```

```
## [1] "T/C"
```

```r
t.haplo.data[which(t.haplo.data == "T/G")]
```

```
## [1] "T/G"
```

Given evidence for heterozygous sites, phase data prior to haplotype clustering. 


```r
vcf.out <- paste(datadir, "chr", chr, "_", site.range, ".vcf", sep="")
system(paste('bcftools view -f PASS -r chr', chr,':' , site.range,' -O v -o ', vcf.out, " ", file, sep=""))
```

To identify haplotypes, I used the program [Beagle](http://faculty.washington.edu/browning/beagle/beagle.html). First, I had to convert the VCF file to BEAGLE format using `vcf2beagle` [utility program](http://faculty.washington.edu/browning/beagle_utilities/utilities.html#vcf2beagle).


```r
beagle.out <- paste(datadir, "chr", chr, "_", site.range, sep="")
system(paste('cat ', vcf.out, ' | java -jar scripts/vcf2beagle.jar ? ', beagle.out, sep=""))
```

I ran BEAGLE on this file for haplotype imputation.


```r
system(paste('java -Xmx1000m -jar scripts/beagle.jar unphased=', paste(beagle.out, ".bgl.gz", sep=""), ' missing=? out=imp', sep=""))
# extract file
system('mv imp* data/.')
system(paste('gunzip ', datadir, "imp.chr", chr, "_", site.range, ".bgl.gz.phased.gz", sep=""))
```

I loaded the imputed haplotypes file into R to determine how many unique haplotypes. First I had to filter out the non-variable site that were in the VCF files.


```r
giraffe.haps <- read.table(paste(datadir, "imp.chr", chr, "_", site.range, ".bgl.gz.phased", sep=""), colClasses="character")
colnames(giraffe.haps) <- giraffe.haps[1,]
rownames(giraffe.haps) <- giraffe.haps[,2]
giraffe.haps <- giraffe.haps[-1,-c(1:2)]

giraffe.haps2 <- matrix(nrow=0, ncol=ncol(giraffe.haps))
colnames(giraffe.haps2) <- colnames(giraffe.haps)

# filter out non-variable sites
rowiter <- 1
for(i in 1:nrow(giraffe.haps)) {
    tmp <- as.character(unlist(giraffe.haps[i,]))
    chars <- unique(tmp)
    if(length(chars)>1) {
      giraffe.haps2 <- rbind(giraffe.haps2, tmp)
      rownames(giraffe.haps2)[rowiter] <- rownames(giraffe.haps)[i]   
      rowiter <- rowiter + 1
      }
    }

dim(giraffe.haps2)
```

```
## [1]  98 524
```

```r
# tranpose
t.giraffe.haps <- t(giraffe.haps2)
t.giraffe.haps[1:5,1:8]
```

```
##         chr8:6779225 chr8:6779281 chr8:6779291 chr8:6779358 chr8:6779389
## HM001   "C"          "A"          "G"          "G"          "T"         
## HM001.1 "C"          "A"          "G"          "G"          "C"         
## HM002   "C"          "A"          "G"          "G"          "C"         
## HM002.1 "C"          "A"          "G"          "G"          "C"         
## HM003   "C"          "A"          "G"          "G"          "T"         
##         chr8:6779435 chr8:6779675 chr8:6779696
## HM001   "T"          "T"          "C"         
## HM001.1 "T"          "T"          "C"         
## HM002   "T"          "T"          "C"         
## HM002.1 "T"          "T"          "C"         
## HM003   "T"          "T"          "C"
```


As a simple first pass of diversity in this gene, I determined the number of unique haplotypes.


```r
haps.even <- list()
haps.odd <- list()

for(i in 1:nrow(t.giraffe.haps)) {
    tmp <- paste(t.giraffe.haps[i,2:ncol(t.giraffe.haps)], collapse="")
    if(i%%2 == 0) haps.even <- c(haps.even, tmp) else {
        haps.odd <- c(haps.odd, tmp)
    }
}

length(haps.even)
```

```
## [1] 262
```

```r
length(haps.odd)
```

```
## [1] 262
```

```r
length(unique(haps.even))
```

```
## [1] 95
```

```r
length(unique(haps.odd))
```

```
## [1] 92
```

```r
haps <- c(haps.even, haps.odd)
length(unique(haps))
```

```
## [1] 103
```

There are 103 unique haplotypes in the 262 samples (524 chromosomes).


I use the R code available from [HaploSuite](http://www.statgen.nus.edu.sg/~software/haplosuite.html) to cluster and visualize the haplotypes.
First, I format the data for use with `haplosim`, which requires filtering to bi-allelic SNPs and re-coding each SNP to 0/1.                                                                                                                                                                 

```r
# convert to HapMap format
haplo.mat <- matrix(nrow=0, ncol=ncol(giraffe.haps2))

for(i in 1:nrow(giraffe.haps2)) {
    tmp <- as.character(unlist(giraffe.haps2[i,]))
    chars <- unique(tmp)
    if(length(chars)==2) {
        tmp[which(tmp==chars[1])] <- 0
        tmp[which(tmp==chars[2])] <- 1
        haplo.mat <- rbind(haplo.mat, as.numeric(tmp))
    }
}

haplo.mat <- as.data.frame(haplo.mat)
haplo.mat[1:10, 1:5]
```

```
##    V1 V2 V3 V4 V5
## 1   0  0  0  0  0
## 2   0  0  0  0  0
## 3   0  0  0  0  0
## 4   0  0  0  0  0
## 5   0  1  1  1  0
## 6   0  0  0  0  0
## 7   0  0  0  0  0
## 8   0  0  0  0  0
## 9   0  0  0  0  0
## 10  0  0  0  0  0
```

```r
position <- giraffe.haps[-1,2]
position <- str_split_fixed(position, ":", n=2)[,2]
```

With the formatted data, I now use the `haplosim` function to cluster haplotypes.



```r
# load functions
source("R/haplosim.R")

distance.measure="physical"

haplosim.out <- haplosim(t(haplo.mat), position, miss.code=9, snp.miss=NULL, focal=FALSE,
    focal.flag=NA, focal.weight=-1, distance.measure=distance.measure, n.hap="auto",
    tolerance=1, sticky=5)
```

```
## Loading required package: lattice
```

```
## [1] "Analysing 524 chromosomes across 98 SNPs."
## [1] "No focal position defined, every SNP contributes equally to haplotype clustering."
## [1] "Input haplotypes include 0 monomorphic SNPs."
## [1] "Finished standardizing the input matrix."
## [1] "Finished calculating the weights for the SNPs."
## [1] "Finished calculating the covariance and penalty matrices"
## [1] "Penalized correlation matrix obtained, begin eigen-decomposition."
## [1] "Automatic detection found 12 possible canonical haplotypes"
## [1] "Finding 12 canonical haplotypes..."
## [1] "Canonical haplotype 1 found with 58 entries."
## [1] "Canonical haplotype 2 found with 32 entries."
## [1] "Canonical haplotype 3 found with 35 entries."
## [1] "Canonical haplotype 4 found with 24 entries."
## [1] "Canonical haplotype 5 found with 20 entries."
## [1] "Canonical haplotype 6 found with 16 entries."
## [1] "Canonical haplotype 7 found with 16 entries."
## [1] "Canonical haplotype 8 found with 14 entries."
## [1] "Canonical haplotype 9 found with 26 entries."
## [1] "Canonical haplotype 10 found with 17 entries."
## [1] "Canonical haplotype 11 found with 18 entries."
## [1] "Canonical haplotype 12 found with 12 entries."
## [1] "Performing simple sorting..."
## [1] "Performing k-means sorting..."
## [1] "Creating mosaics of unmapped haplotypes..."
## [1] "Done! Recommend the use of HAPVISUAL to visualize the results."
```

```r
names(haplosim.out)
```

```
##  [1] "unsorted"          "simple.sort"       "stepwise.sort"    
##  [4] "kmeans.unsorted"   "kmeans.sort"       "haplosim.unsorted"
##  [7] "haplosim.sort"     "hap.group"         "hap.ranking"      
## [10] "simple.ranking"    "stepwise.ranking"  "kmeans.ranking"   
## [13] "focal.flag"        "focal.pos"         "n.chr"            
## [16] "n.snp"             "eigenvalues"       "eigenvectors"     
## [19] "n.hap"
```

```r
table(haplosim.out$hap.group)
```

```
## 
##  -1   1   2   3   4   5   6   7   8   9  10  11  12 
## 236  58  32  35  24  20  16  16  14  26  17  18  12
```

Haplotype clustering reveals 12 clusters, with most chromosomes unassigned (-1) and the other groups consisting of 12 to 58 chromosomes.

Next, I check to see if any of the accessions have chromosomes assigned to different haplotype clusters.


```r
h2 <- data.frame(
  HM = rep(colnames(giraffe.haps2)[seq(1, length(colnames(giraffe.haps2)), by=2)], each = 2),
  chrom = rep(c(1,2), times = length(haplosim.out$hap.group)/2),
  hap = haplosim.out$hap.group)

head(h2)
```

```
##      HM chrom hap
## 1 HM001     1  -1
## 2 HM001     2  -1
## 3 HM002     1  -1
## 4 HM002     2  -1
## 5 HM003     1   3
## 6 HM003     2   3
```

```r
# find which accessions have chromosome assigned to different haplotype clusters
h2.mismatch <- ddply(h2, .(HM), summarize, eq = all.equal(hap[1], hap[2]))
                 
mismatch <- h2.mismatch[which(h2.mismatch$eq != "TRUE"), "HM"]

h2[which(h2$HM %in%  mismatch), ]
```

```
##        HM chrom hap
## 13  HM007     1  -1
## 14  HM007     2  10
## 193 HM108     1   2
## 194 HM108     2  -1
## 209 HM119     1   3
## 210 HM119     2  -1
## 285 HM165     1  -1
## 286 HM165     2   2
```

Only 4 accessions have chromosome with different haplotypes, and all of these consist of pairs that have one clustered chromosome paired with an unassigned chromosome. I can thus reduce the chromosomes to a single haplotype for each accession.


```r
acc.hap <- data.frame(HM = levels(h2$HM), hap = NA)
  
for(i in levels(h2$HM)) {
  hm <- h2[which(h2$HM == i), ]
  if(hm$hap[1] == hm$hap[2]) 
    acc.hap[which(acc.hap$HM == i), "hap"] <- hm$hap[1] else {
      acc.hap[which(acc.hap$HM == i), "hap"] <- max(hm$hap)
    }
  }

table(acc.hap$hap)
```

```
## 
##  -1   1   2   3   4   5   6   7   8   9  10  11  12 
## 116  29  17  18  12  10   8   8   7  13   9   9   6
```




Next, I visualize these clusters using code provided from [HaploSuite](http://www.statgen.nus.edu.sg/~software/haplosuite.html). 



```
## pdf 
##   2
```

```
## [1] "Separate chromosomes into 1 subpopulations..."
```

```
## [1] "Clustering haplotypes using method haplosim..."
```

```
## [1] "Plotting clustering for subpopulation 1 with 524 chromosomes."
```

![image](haplovisual.png)
                                                                                                                                                                                                    
