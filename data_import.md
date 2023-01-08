R Notebook
================

``` bash
#wget -i data/url_data
```

``` r
library(rmarkdown)
library(knitr)
library(phyloseq)
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(DECIPHER)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

``` r
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

``` r
library(ggplot2)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library(shiny)
library(miniUI)
library(caret)
```

    ## Loading required package: lattice

``` r
library(pls)
```

    ## 
    ## Attaching package: 'pls'

    ## The following object is masked from 'package:caret':
    ## 
    ##     R2

    ## The following object is masked from 'package:ape':
    ## 
    ##     mvr

    ## The following object is masked from 'package:stats':
    ## 
    ##     loadings

``` r
library(e1071)
library(ggplot2)
library(randomForest)
```

    ## randomForest 4.7-1.1

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library("plyr"); packageVersion("plyr")
```

    ## 
    ## Attaching package: 'plyr'

    ## The following object is masked from 'package:XVector':
    ## 
    ##     compact

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     desc

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## [1] '1.8.7'

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following object is masked from 'package:randomForest':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggrepel)
library(devtools)
```

    ## Loading required package: usethis

``` r
library(reshape2)
library(PMA)
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
library(ggnetwork)
library(intergraph)
library(scales)
library(genefilter)
library(impute)
library(phyloseqGraphTest)
library(Biostrings)
library(dada2); packageVersion("dada2")
```

    ## [1] '1.24.0'

``` r
path <- "~/CC3finalll/data/url_data" 
list.files(path)
```

    ## character(0)

``` r
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
```

``` r
#plotQualityProfile(fnFs[1:2])
```

``` r
#plotQualityProfile(fnRs[1:2])
```

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

``` r
#errF <- learnErrors(filtFs, multithread=TRUE)
```

``` r
#errR <- learnErrors(filtRs, multithread=TRUE)
```

``` r
#plotErrors(errF, nominalQ=TRUE)
```

``` r
#dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

``` r
#dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

``` r
#mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#head(mergers[[1]])
```

``` r
#seqtab <- makeSequenceTable(mergers)
#dim(seqtab)
```

``` r
#table(nchar(getSequences(seqtab)))
```

``` r
#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)
```

``` r
#sum(seqtab.nochim)/sum(seqtab)
```

``` r
#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track) <- sample.names
#head(track)
```

``` bash
#wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
```

``` r
#taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
```

``` r
#taxa.print <- taxa 
#rownames(taxa.print) <- NULL
#head(taxa.print)
```

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.40.0'

``` r
library(Biostrings); packageVersion("Biostrings")
```

    ## [1] '2.64.1'

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.3.6'

``` r
#theme_set(theme_bw())
```

``` r
#samples.out <- rownames(seqtab.nochim)
#prof <- sapply(strsplit(samples.out, "_"), `[`, 2)
#s_prof <- substr(prof,1,1)
#day <- as.character(sapply(strsplit(samples.out, "_"), `[`, 3))
#samdf <- data.frame(prof=s_prof, Day=day)
#samdf$prof <- s_prof
#samdf$Saison <- "Ete"
#samdf$Saison[samdf$Day > "10sept14"] <- "Hiver"
#rownames(samdf) <- samples.out
```

``` r
#ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
              # sample_data(samdf), 
               #tax_table(taxa))
```

``` r
#plot_richness(ps, x="Saison", measures=c("Shannon", "Simpson"), color="prof")
```

``` r
#ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
#ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
```

``` r
#plot_ordination(ps.prop, ord.nmds.bray, color="Saison", shape="prof", title="Bray NMDS")
```

``` r
#top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
#ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#ps.top20 <- prune_taxa(top20, ps.top20)
#plot_bar(ps.top20, x="Saison", fill="Class") + facet_wrap(~prof, scales="free_x")
```

``` r
#top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
#ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#ps.top20 <- prune_taxa(top20, ps.top20)
#plot_bar(ps.top20, x="Saison", fill="Family") + facet_wrap(~prof, scales="free_x")
```

``` r
#top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
#ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#ps.top20 <- prune_taxa(top20, ps.top20)
#plot_bar(ps.top20, x="Saison", fill="Genus") + facet_wrap(~prof, scales="free_x")
```

``` r
#top40 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:40]
#ps.top40 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#ps.top40 <- prune_taxa(top40, ps.top40)
#plot_bar(ps.top40, x="Saison", fill="Family") + facet_wrap(~prof, scales="free_x")
```
