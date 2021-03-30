eTax\_analysis\_ensemble\_comps
================
D Catlett
8/24/2020

# Overview

This is an analysis of ensemble taxonomic assignments computed from 4 different assignment algorithm-reference database combo's. bayes and idtax are the algorithms, pr2 and silva are the databases. Starting from initial taxonomy tables from Mar20 dada2 run.

We'll start with a large data set of 18S-V9 ASVs inferred by dada2. The data come from several studies of marine protist communities, in addition to some extensive mock community validations. The primers are 'universal', meaning we have both prokaryotes and eukaryotes in our data set. The raw MiSeq data were processed with dada2 and four different batches of taxonomic assignments were made by running the bayesian classifer (via dada2's assignTaxonomy function) and the idtaxa algorithm (via the DECIPHER package) against both the Silva nr SSU v138 database and the Protistan Ribosomal Reference database v4.12.0. No bootstrap thresholds have been designated yet - we will apply those with the ensembleTax pre-processing functions below.

First we'll load the data and the ensembleTax package, and use the ensembleTax pre-processing functions to format the taxonomy tables properly and to set our bootstrap thresholds. Included in this repo is a "asv rubric" fasta file that has ASV sequences and arbitrary identifying names.

## Data pre-processing

``` r
rm(list=ls())
library("ensembleTax")
packageVersion("ensembleTax")
```

    ## [1] '1.0.2'

``` r
library("Biostrings")
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

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

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library("stringr")

bayes.pr2 <- readRDS("initial_tax_tabs/bayes_pr2_0boot_Mar20.rds")
bayes.silva <- readRDS("initial_tax_tabs/bayes_silva_0boot_Mar20.rds")
idtax.pr2 <- readRDS("initial_tax_tabs/idtax_pr2_0boot_Mar20.rds")
idtax.silva <- readRDS("initial_tax_tabs/idtax_silva_0boot_Mar20.rds")

asv.rubric <- readDNAStringSet("asv_rubric.fasta")
bayes.pr2 <- bayestax2df(bayes.pr2, db = "pr2", boot = 60, rubric = asv.rubric, return.conf = FALSE)
bayes.silva <- bayestax2df(bayes.silva, db = "silva", boot = 60, rubric = asv.rubric, return.conf = FALSE)
idtax.pr2 <- idtax2df(idtax.pr2, db = "pr2", boot = 50, rubric = asv.rubric, return.conf = FALSE)
idtax.silva <- idtax2df(idtax.silva, db = "silva", boot = 50, rubric = asv.rubric, return.conf = FALSE)

head(bayes.pr2)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##     kingdom     supergroup       division           class             order
    ## 1 Eukaryota   Opisthokonta        Metazoa      Arthropoda         Crustacea
    ## 2 Eukaryota  Stramenopiles     Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 3 Eukaryota      Alveolata Dinoflagellata     Syndiniales              <NA>
    ## 4 Eukaryota Archaeplastida   Streptophyta   Embryophyceae   Embryophyceae_X
    ## 5 Eukaryota           <NA>           <NA>            <NA>              <NA>
    ## 6 Eukaryota           <NA>           <NA>            <NA>              <NA>
    ##             family            genus                      species
    ## 1      Maxillopoda      Paracalanus           Paracalanus_parvus
    ## 2   Raphid-pennate Pseudo-nitzschia Pseudo-nitzschia_fraudulenta
    ## 3             <NA>             <NA>                         <NA>
    ## 4 Embryophyceae_XX             <NA>                         <NA>
    ## 5             <NA>             <NA>                         <NA>
    ## 6             <NA>             <NA>                         <NA>

``` r
head(bayes.silva)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##      domain             phylum             class                order
    ## 1 Eukaryota         Arthropoda       Maxillopoda            Calanoida
    ## 2 Eukaryota           Diatomea Bacillariophyceae Bacillariophyceae_or
    ## 3 Eukaryota      Protalveolata       Syndiniales       Syndiniales_or
    ## 4 Eukaryota Phragmoplastophyta       Embryophyta                 <NA>
    ## 5      <NA>               <NA>              <NA>                 <NA>
    ## 6 Eukaryota               <NA>              <NA>                 <NA>
    ##                 family            genus
    ## 1                 <NA>             <NA>
    ## 2 Bacillariophyceae_fa Pseudo-nitzschia
    ## 3                 <NA>             <NA>
    ## 4                 <NA>             <NA>
    ## 5                 <NA>             <NA>
    ## 6                 <NA>             <NA>

``` r
head(idtax.pr2)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##     kingdom     supergroup       division           class             order
    ## 1 Eukaryota   Opisthokonta        Metazoa      Arthropoda         Crustacea
    ## 2 Eukaryota  Stramenopiles     Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 3 Eukaryota      Alveolata Dinoflagellata            <NA>              <NA>
    ## 4 Eukaryota Archaeplastida   Streptophyta   Embryophyceae   Embryophyceae_X
    ## 5 Eukaryota           <NA>           <NA>            <NA>              <NA>
    ## 6 Eukaryota      Alveolata Dinoflagellata     Syndiniales     Dino-Group-II
    ##                   family                    genus                      species
    ## 1            Maxillopoda                     <NA>                         <NA>
    ## 2         Raphid-pennate         Pseudo-nitzschia Pseudo-nitzschia_fraudulenta
    ## 3                   <NA>                     <NA>                         <NA>
    ## 4       Embryophyceae_XX                     <NA>                         <NA>
    ## 5                   <NA>                     <NA>                         <NA>
    ## 6 Dino-Group-II-Clade-14 Dino-Group-II-Clade-14_X Dino-Group-II-Clade-14_X_sp.

``` r
head(idtax.silva)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##      domain             phylum             class                order
    ## 1 Eukaryota         Arthropoda       Maxillopoda            Calanoida
    ## 2 Eukaryota           Diatomea Bacillariophyceae Bacillariophyceae_or
    ## 3 Eukaryota      Protalveolata       Syndiniales       Syndiniales_or
    ## 4 Eukaryota Phragmoplastophyta       Embryophyta                 <NA>
    ## 5      <NA>               <NA>              <NA>                 <NA>
    ## 6      <NA>               <NA>              <NA>                 <NA>
    ##                  family genus
    ## 1                  <NA>  <NA>
    ## 2  Bacillariophyceae_fa  <NA>
    ## 3 Syndiniales Group III  <NA>
    ## 4                  <NA>  <NA>
    ## 5                  <NA>  <NA>
    ## 6                  <NA>  <NA>

The data looks good and at first glance it looks like taxonomic assignments are mostly in agreement across the taxonomy tables where they are assigned, though there are differences in resolution. Here we'll subset ASVs from a particular project so they're all from the same data set:

``` r
# sanity check
identical(bayes.pr2[, 1:2], bayes.silva[, 1:2])
```

    ## [1] TRUE

``` r
identical(bayes.pr2[, 1:2], idtax.silva[, 1:2])
```

    ## [1] TRUE

``` r
identical(bayes.pr2[, 1:2], idtax.pr2[, 1:2])
```

    ## [1] TRUE

``` r
seqtab <- readRDS("initial_tax_tabs/seqtab_nochime_Mar20.rds")
samsub <- readRDS("initial_tax_tabs/sample_subsets.rds")
subber <- c(samsub$pbsurf, samsub$pbdeep)
seqtab <- seqtab[ row.names(seqtab) %in% subber , ]
rm.me <- which(colSums(seqtab) == 0)
# also remove <90 or >180
rm.me2 <- which(str_length(colnames(seqtab)) < 90 | str_length(colnames(seqtab)) > 180)
rm.me <- unique(c(rm.me, rm.me2))
asv.sub <- colnames(seqtab)[-rm.me]

bayes.pr2 <- bayes.pr2[bayes.pr2$ASV %in% asv.sub , ]
bayes.silva <- bayes.silva[bayes.silva$ASV %in% asv.sub , ]
idtax.pr2 <- idtax.pr2[idtax.pr2$ASV %in% asv.sub , ]
idtax.silva <- idtax.silva[idtax.silva$ASV %in% asv.sub , ]

identical(bayes.pr2[, 1:2], bayes.silva[, 1:2])
```

    ## [1] TRUE

``` r
identical(bayes.pr2[, 1:2], idtax.silva[, 1:2])
```

    ## [1] TRUE

``` r
identical(bayes.pr2[, 1:2], idtax.pr2[, 1:2])
```

    ## [1] TRUE

## Taxonomy mapping

Here we'll run taxmapper to map the pr2 tables onto silva, and the silva tables onto pr2. We don't care about the mapping results other than returning the mapped taxonomy tables (you can see more on that on the mapping doc), so will set streamline = TRUE. We'll be as permissive as possible to try to get the most-resolved mapped taxonomic assignments we can by using *synonym.file = "default"* and *ignore.format = TRUE*.

``` r
bayes.pr2.mapped2silva <- taxmapper(bayes.pr2, 
                                     tt.ranks = colnames(bayes.pr2)[3:ncol(bayes.pr2)],
                                     tax2map2 = "Silva",
                                     ignore.format = TRUE,
                                     synonym.file = "default",
                                     streamline = TRUE,
                                     outfilez = NULL)

bayes.silva.mapped2pr2 <- taxmapper(bayes.silva, 
                                     tt.ranks = colnames(bayes.silva)[3:ncol(bayes.silva)],
                                     tax2map2 = "pr2",
                                     exceptions = c("Bacteria", "Archaea"),
                                     ignore.format = TRUE,
                                     synonym.file = "default",
                                     streamline = TRUE,
                                     outfilez = NULL)

idtax.pr2.mapped2silva <- taxmapper(idtax.pr2, 
                                     tt.ranks = colnames(idtax.pr2)[3:ncol(idtax.pr2)],
                                     tax2map2 = "Silva",
                                     ignore.format = TRUE,
                                     synonym.file = "default",
                                     streamline = TRUE,
                                     outfilez = NULL)

idtax.silva.mapped2pr2 <- taxmapper(idtax.silva, 
                                     tt.ranks = colnames(idtax.silva)[3:ncol(idtax.silva)],
                                     tax2map2 = "pr2",
                                     exceptions = c("Bacteria", "Archaea"),
                                     ignore.format = TRUE,
                                     synonym.file = "default",
                                     streamline = TRUE,
                                     outfilez = NULL)
```

We've now mapped each taxonomy table onto a new taxonomic nomenclature -- the 2 silva tables were mapped onto pr2's nomenclature, and the 2 pr2 tables were mapped onto silva's nomenclature. Let's look at the mapped taxonomy tables:

``` r
head(bayes.pr2.mapped2silva)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##      domain        phylum             class                order
    ## 1 Eukaryota    Arthropoda       Maxillopoda                 <NA>
    ## 2 Eukaryota      Diatomea Bacillariophyceae Bacillariophyceae_or
    ## 3 Eukaryota Protalveolata       Syndiniales                 <NA>
    ## 4 Eukaryota          <NA>              <NA>                 <NA>
    ## 5 Eukaryota          <NA>              <NA>                 <NA>
    ## 6 Eukaryota          <NA>              <NA>                 <NA>
    ##                 family            genus
    ## 1                 <NA>             <NA>
    ## 2 Bacillariophyceae_fa Pseudo-nitzschia
    ## 3                 <NA>             <NA>
    ## 4                 <NA>             <NA>
    ## 5                 <NA>             <NA>
    ## 6                 <NA>             <NA>

``` r
head(bayes.silva.mapped2pr2)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##     kingdom    supergroup       division           class             order
    ## 1 Eukaryota  Opisthokonta        Metazoa      Arthropoda         Crustacea
    ## 2 Eukaryota Stramenopiles     Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 3 Eukaryota     Alveolata Dinoflagellata     Syndiniales              <NA>
    ## 4 Eukaryota          <NA>           <NA>            <NA>              <NA>
    ## 5      <NA>          <NA>           <NA>            <NA>              <NA>
    ## 6 Eukaryota          <NA>           <NA>            <NA>              <NA>
    ##           family            genus species
    ## 1    Maxillopoda             <NA>    <NA>
    ## 2 Raphid-pennate Pseudo-nitzschia    <NA>
    ## 3           <NA>             <NA>    <NA>
    ## 4           <NA>             <NA>    <NA>
    ## 5           <NA>             <NA>    <NA>
    ## 6           <NA>             <NA>    <NA>

``` r
head(idtax.pr2.mapped2silva)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##      domain         phylum             class                order
    ## 1 Eukaryota     Arthropoda       Maxillopoda                 <NA>
    ## 2 Eukaryota       Diatomea Bacillariophyceae Bacillariophyceae_or
    ## 3 Eukaryota Dinoflagellata              <NA>                 <NA>
    ## 4 Eukaryota           <NA>              <NA>                 <NA>
    ## 5 Eukaryota           <NA>              <NA>                 <NA>
    ## 6 Eukaryota  Protalveolata       Syndiniales                 <NA>
    ##                 family            genus
    ## 1                 <NA>             <NA>
    ## 2 Bacillariophyceae_fa Pseudo-nitzschia
    ## 3                 <NA>             <NA>
    ## 4                 <NA>             <NA>
    ## 5                 <NA>             <NA>
    ## 6                 <NA>             <NA>

``` r
head(idtax.silva.mapped2pr2)
```

    ##       svN
    ## 1     sv1
    ## 2    sv10
    ## 3   sv100
    ## 4  sv1000
    ## 5 sv10000
    ## 6 sv10001
    ##                                                                                                                                   ASV
    ## 1     GCTACTACCGATTGAACATTTTAGTGAGGTCCTCGGACTGTGAGCCAGGCGGGTCGCCCTGCCTGGTCTACGGGAAGACGACCAAACTGTAGTGTTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC
    ## 2     GCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGATTAGTTTCCTTTATTGGAAGGTAGTTATGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCT
    ## 3     GCTCCTACCGATTGAGTGGTCCGGTGAATAATTCGGACTGGTGCCGATTTCGGTTCTCCGAGTTCGGCGCTGGGAAGTCTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 4 GCTCCTACCGATTGAATGGTCCGGTGAAGTGTTCGGATCGTGGCGACGTGGGCGGTTCGCTGCCTGCGACGTCGCGAGAAGTCCACTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCC
    ## 5                                  GTTTCTTCCGACTAATGCTTTATGTGAGTGTCACGGATTTTAAAGAAGTGCTGTGAACATTGAGTATCGGAGGAAGAAAAAGTCGTAACAAGGTTATC
    ## 6      GCTGCTACCGATTGAGTGTCCTGGTGAATTATTTGGACCGGCAGTAATTCGAGTTTCTCGATTTACAGCTGGAAAATCTTGTAAACCCTGACACTTAGAGGAAGCAGAAGTCGTAACAAGGTTTCC
    ##     kingdom    supergroup       division           class             order
    ## 1 Eukaryota  Opisthokonta        Metazoa      Arthropoda         Crustacea
    ## 2 Eukaryota Stramenopiles     Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 3 Eukaryota     Alveolata Dinoflagellata     Syndiniales              <NA>
    ## 4 Eukaryota          <NA>           <NA>            <NA>              <NA>
    ## 5      <NA>          <NA>           <NA>            <NA>              <NA>
    ## 6      <NA>          <NA>           <NA>            <NA>              <NA>
    ##           family genus species
    ## 1    Maxillopoda  <NA>    <NA>
    ## 2 Raphid-pennate  <NA>    <NA>
    ## 3           <NA>  <NA>    <NA>
    ## 4           <NA>  <NA>    <NA>
    ## 5           <NA>  <NA>    <NA>
    ## 6           <NA>  <NA>    <NA>

Looks good. ensembleTax time:

``` r
## mapped to pr2's:
xx <- list(bayes.pr2, idtax.pr2, bayes.silva.mapped2pr2, idtax.silva.mapped2pr2)
names(xx) <- c("bayes-pr2", "idtax-pr2", "bayes-silva", "idtax-silva")
et.all4.pr2.wna <- ensembleTax(xx, count.na = TRUE, 
                               tiebreakz = NULL, assign.threshold = 0)
et.all4.pr2.nona <- ensembleTax(xx, count.na = FALSE, 
                               tiebreakz = c("idtax-pr2", "idtax-silva", "bayes-pr2"))

## mapped to silva
xx.silva <- list(bayes.pr2.mapped2silva, idtax.pr2.mapped2silva, bayes.silva, idtax.silva)
rn <- rownames(idtax.silva)
for (i in 1:length(xx.silva)) {
  rownames(xx.silva[[i]]) <- rn
}
names(xx.silva) <- c("bayes-pr2", "idtax-pr2", "bayes-silva", "idtax-silva")
et.all4.silva.wna <- ensembleTax(xx.silva, count.na = TRUE, ranknames = colnames(idtax.silva[, 3:ncol(idtax.silva)]),
                               tiebreakz = NULL, assign.threshold = 0)
et.all4.silva.nona <- ensembleTax(xx.silva, count.na = FALSE, ranknames = colnames(idtax.silva[, 3:ncol(idtax.silva)]),
                               tiebreakz = c("idtax-silva", "idtax-pr2", "bayes-silva"))
```

Ok now remove proks, metazoa, fungi, streptophyta

``` r
rm.i <- unique(c(which(et.all4.pr2.nona$kingdom %in% c("Bacteria", "Archaea")),
            which(et.all4.pr2.nona$kingdom %in% c("Bacteria", "Archaea")),
            which(et.all4.pr2.nona$division %in% c("Metazoa", "Fungi", "Streptophyta"))))
bayes.silva.mapped2pr2 <- bayes.silva.mapped2pr2[-rm.i , ]
idtax.silva.mapped2pr2 <- idtax.silva.mapped2pr2[-rm.i , ]
bayes.pr2 <- bayes.pr2[-rm.i , ]
idtax.pr2 <- idtax.pr2[-rm.i , ]
et.all4.pr2.nona <- et.all4.pr2.nona[-rm.i , ]
et.all4.pr2.wna <- et.all4.pr2.wna[-rm.i , ]

## the ones mapped to silva too:
bayes.silva <- bayes.silva[-rm.i , ]
idtax.silva <- idtax.silva[-rm.i , ]
bayes.pr2.mapped2silva <- bayes.pr2.mapped2silva[-rm.i , ]
idtax.pr2.mapped2silva <- idtax.pr2.mapped2silva[-rm.i , ]
et.all4.silva.nona <- et.all4.silva.nona[-rm.i , ]
et.all4.silva.wna <- et.all4.silva.wna[-rm.i , ]
```

Compute and plot % ASVs unassigned at each rank:

``` r
library("ggplot2")
library("reshape2")

nasum <- function(taxdf){
    notuz <- nrow(taxdf)
    x <- is.na(taxdf)
    ii <- colSums(x) / notuz
    return(ii)
}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
# with the individuals and 4-way ensembles:
xx.pr2 <- list(idtax.pr2, bayes.pr2,
               idtax.silva.mapped2pr2, bayes.silva.mapped2pr2, 
               et.all4.pr2.wna, et.all4.pr2.nona)
names(xx.pr2) <- c("idtax-pr2", "bayes-pr2",  
                 "idtax-silva","bayes-silva",
                 "ensemble-acc", "ensemble-res")
xx.pr2 <- lapply(xx.pr2, function(x) x[, -c(1,2)])
tina <- lapply(xx.pr2, nasum)
yaboi <- matrix(unlist(tina), nrow=length(xx.pr2), byrow=TRUE)
rownames(yaboi) <- names(xx.pr2)
colnames(yaboi) <- colnames(xx.pr2[[1]])
yaboi <- as.data.frame(t(yaboi))
yaboi$rankz <- rownames(yaboi)
yaboi <- melt(yaboi, id.vars = "rankz")
p.pr2.2 <- ggplot(yaboi, aes(fill = variable, x = rankz, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Taxonomic Rank", y = "Proportion of ASVs Unassigned") + 
    scale_x_discrete(limits = colnames(xx.pr2[[1]])) + 
    scale_y_continuous(limits = seq(0,1.05,0.1), expand = c(0,0)) +
    coord_cartesian(ylim = c(0, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_fill_manual(name = "Taxonomy table",
                      breaks = names(xx.pr2),
                      values = cbPalette)


# silvas:
xx.silva <- list(idtax.pr2.mapped2silva, bayes.pr2.mapped2silva,
               idtax.silva, bayes.silva, 
               et.all4.silva.wna, et.all4.silva.nona)
names(xx.silva) <- c("idtax-pr2", "bayes-pr2",  
                 "idtax-silva","bayes-silva",
                 "ensemble-acc", "ensemble-res")
xx.silva <- lapply(xx.silva, function(x) x[, -c(1,2)])
tina <- lapply(xx.silva, nasum)
yaboi <- matrix(unlist(tina), nrow=length(xx.silva), byrow=TRUE)
rownames(yaboi) <- names(xx.silva)
colnames(yaboi) <- colnames(xx.silva[[1]])
yaboi <- as.data.frame(t(yaboi))
yaboi$rankz <- rownames(yaboi)
yaboi <- melt(yaboi, id.vars = "rankz")
p.silva.2 <- ggplot(yaboi, aes(fill = variable, x = rankz, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Taxonomic Rank", y = "Proportion of ASVs Unassigned") + 
    scale_x_discrete(limits = colnames(xx.silva[[1]])) + 
    scale_y_continuous(limits = seq(0,1.05,0.1), expand = c(0,0)) +
    coord_cartesian(ylim = c(0, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_fill_manual(name = "Taxonomy",
                      breaks = names(xx.pr2),
                      values = cbPalette)
```

Nice. Let's do some comparisons of the assignments between the ensembles and the individual taxonomy tables:

``` r
library("dplyr")
# fcn that compares two taxonomy tables:
tblcomper <- function(y,x) {
  perf <- dplyr::intersect(x, y) # perfectly-matching rows (assignments)

  tmp.x <- dplyr::setdiff(x, perf)
  tmp.y <- dplyr::setdiff(y, perf)
  if (!identical(tmp.x[, 1], tmp.y[, 1])) {
    tmp.x <- sort_my_taxtab(tmp.x, ranknames = colnames(tmp.x)[3:ncol(tmp.x)])
    tmp.y <- sort_my_taxtab(tmp.y, ranknames = colnames(tmp.y)[3:ncol(tmp.x)])
  }
  xna <- is.na(tmp.x)
  yna <- is.na(tmp.y)
  # This warning happens 1 million times when there's no NA's. It doesn't matter b/c inf is always bigger so suppressing.
  ## Warning in min(which(z)): no non-missing arguments to min; returning Inf
  x.minna <- apply(xna, MARGIN = 1, function(z) suppressWarnings(min(which(z))))
  y.minna <- apply(yna, MARGIN = 1, function(z) suppressWarnings(min(which(z))))
  x.mo <- which(x.minna > y.minna)
  y.mo <- which(y.minna > x.minna)
  
  yunder.i <- c()
  yover.i <- c()
  mis.i <- c()
  
  # subset where x has more resolved assignments and only to cols where 
  # both have assignments, then see if names match
  yunder <- 0
  ymis <- 0
  if (length(x.mo) > 0){
    for (i in 1:length(x.mo)){
      if ((y.minna[x.mo[i]]-1) < 3){ # if kingdom is unassigned just add to under-classification
        yunder <- yunder+1
        # yunder.i <- append(yunder.i, tmp.x$svN)
      } else {
        tmp.tmpx <- tmp.x[x.mo[i] , 3:(y.minna[x.mo[i]]-1)]
        tmp.tmpy <- tmp.y[x.mo[i] , 3:(y.minna[x.mo[i]]-1)]
    
        if (all(tmp.tmpx == tmp.tmpy)) {
          yunder <- yunder+1
          # yunder.i <- append(yunder.i, tmp.x$svN)
        } else {
          ymis <- ymis+1
          # mis.i <- append(mis.i, tmp.x$svN)
        }
      }
    }
  }
  
  # repeat above where y is more resolved than x:
  yover <- 0
  xmis <- 0
  if (length(y.mo) > 0){
    for (i in 1:length(y.mo)){
      if ((x.minna[y.mo[i]]-1) < 3){ # if kingdom is unassigned just add to under-classification
        yover <- yover+1
        # yover.i <- append(yover.i, tmp.x$svN)
      } else {
        tmp.tmpx <- tmp.x[y.mo[i] , 3:(x.minna[y.mo[i]]-1)]
        tmp.tmpy <- tmp.y[y.mo[i] , 3:(x.minna[y.mo[i]]-1)]
    
        if (all(tmp.tmpx == tmp.tmpy)) {
          yover <- yover+1
          # yover.i <- append(yover.i, tmp.x$svN)
        } else {
          xmis <- xmis+1
          # mis.i <- append(mis.i, tmp.x$svN)
        }
      }
    }
  }
  
  if (yover+xmis != length(y.mo) || yunder+ymis != length(x.mo)) {
    stop("somethings wrong i think")
  }
  
  perf <- nrow(perf)
  # whatever's left is where both tables have the same ranks named, but different names. this is a misclassification:
  moremis <- nrow(x) - (xmis+ymis+yover+yunder+perf)
  result <- data.frame(all.match = c(perf), mis = c(ymis+xmis+moremis), over = c(yover), under = c(yunder))
  if (sum(result) == nrow(x)) {
    return(result)
    # return(list(result, mis.i, yover.i, yunder.i, perf.i))
  } else {
    stop("noooooooooooooo")
  }
}

### this is the pr2 comparisons... scroll past plot saving for silva
# make a list of all taxonomy tables:
tbl.list <- list(idtax.pr2, bayes.pr2, 
                 idtax.silva.mapped2pr2, bayes.silva.mapped2pr2,
                 et.all4.pr2.wna)
all.comp.nona <- lapply(tbl.list, FUN = tblcomper, x = et.all4.pr2.nona)
all.comp.nona <- base::do.call(base::rbind.data.frame, all.comp.nona)
row.names(all.comp.nona) <- c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-acc")
all.comp.nona <- all.comp.nona / rowSums(all.comp.nona)

tbl.list <- list(idtax.pr2, bayes.pr2, 
                 idtax.silva.mapped2pr2, bayes.silva.mapped2pr2,
                 et.all4.pr2.nona)
all.comp.wna <- lapply(tbl.list, FUN = tblcomper, x = et.all4.pr2.wna)
all.comp.wna <- base::do.call(base::rbind.data.frame, all.comp.wna)
row.names(all.comp.wna) <- c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-res")
all.comp.wna <- all.comp.wna / rowSums(all.comp.wna)

all.comp.nona$tbl <- row.names(all.comp.nona)
all.comp.wna$tbl <- row.names(all.comp.wna)

plt.all.comp.wna <- melt(all.comp.wna)
plt.all.comp.nona <- melt(all.comp.nona)

# plot the results:
p.comp1.pr2 <- ggplot(plt.all.comp.wna, aes(fill = tbl, x = variable, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Relative to ensemble-acc", y = "Proportion of ASVs") + 
    scale_x_discrete(breaks = c("all.match","mis","over","under"), 
                   labels = c("Agree", "Misclassified", "Overclassified","Underclassified")) +
    scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                   limits = c(0, 1), expand = c(0,0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_fill_manual(name = "Taxonomy",
                      breaks = c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-res"),
                      values = cbPalette[c(1:4,6)])
print(p.comp1.pr2)
```

![](eTax_analysis_ensemble_comps2_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
# plot the results:
p.comp2.pr2 <- ggplot(plt.all.comp.nona, aes(fill = tbl, x = variable, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Relative to ensemble-res", y = "Proportion of ASVs") + 
    scale_x_discrete(breaks = c("all.match","mis","over","under"), 
                   labels = c("Agree", "Misclassified", "Overclassified","Underclassified")) +
    scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                   limits = c(0, 1), expand = c(0,0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_fill_manual(name = "Taxonomy",
                      breaks = c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-acc"),
                      values = cbPalette[c(1:5)])
  
print(p.comp2.pr2)
```

![](eTax_analysis_ensemble_comps2_files/figure-markdown_github/unnamed-chunk-8-2.png)

``` r
library("cowplot")
p.all <- plot_grid(
  p.pr2.2 + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  p.comp1.pr2 + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  p.comp2.pr2 + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  align = 'h',
  labels = c("(A)", "(B)", "(C)"),
  axis = 'l',
  hjust=-1,
  nrow=1
)
```

    ## Warning: Removed 45 rows containing missing values (geom_bar).

``` r
legend_b <- get_legend(p.pr2.2 + 
                         theme(legend.position="bottom",
                              legend.text=element_text(size=14),
                              legend.title = element_text(size=14),
                              legend.justification="center",
                              legend.box.margin = unit(c(0.1, 1, 0.1, 1),"cm")))
```

    ## Warning: Removed 45 rows containing missing values (geom_bar).

``` r
p.all <- plot_grid(p.all, legend_b, nrow=2, rel_heights=c(1, 0.3))
ggsave("all_pr2_plots_protists", plot = p.all, device = "pdf", width = 18, height = 7.5)
```

Let's do silva comparisons and make a 6-panel fig:

``` r
# make a list of all taxonomy tables:
tbl.list <- list(idtax.pr2.mapped2silva, bayes.pr2.mapped2silva,
                 idtax.silva, bayes.silva,
                 et.all4.silva.wna)
all.comp.nona <- lapply(tbl.list, FUN = tblcomper, x = et.all4.silva.nona)
all.comp.nona <- base::do.call(base::rbind.data.frame, all.comp.nona)
row.names(all.comp.nona) <- c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-acc")
all.comp.nona <- all.comp.nona / rowSums(all.comp.nona)

tbl.list <- list(idtax.pr2.mapped2silva, bayes.pr2.mapped2silva,
                 idtax.silva, bayes.silva,
                 et.all4.silva.nona)
all.comp.wna <- lapply(tbl.list, FUN = tblcomper, x = et.all4.silva.wna)
all.comp.wna <- base::do.call(base::rbind.data.frame, all.comp.wna)
row.names(all.comp.wna) <- c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-res")
all.comp.wna <- all.comp.wna / rowSums(all.comp.wna)

all.comp.nona$tbl <- row.names(all.comp.nona)
all.comp.wna$tbl <- row.names(all.comp.wna)

plt.all.comp.wna <- melt(all.comp.wna)
```

    ## Using tbl as id variables

``` r
plt.all.comp.nona <- melt(all.comp.nona)
```

    ## Using tbl as id variables

``` r
# plot the results:
p.comp1.silva <- ggplot(plt.all.comp.wna, aes(fill = tbl, x = variable, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Relative to ensemble-acc", y = "Proportion of ASVs") + 
    scale_x_discrete(breaks = c("all.match","mis","over","under"), 
                   labels = c("Agree", "Misclassified", "Overclassified","Underclassified")) +
    scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                   limits = c(0, 1), expand = c(0,0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_fill_manual(name = "Taxonomy",
                      breaks = c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-res"),
                      values = cbPalette[c(1:4,6)])
print(p.comp1.silva)
```

![](eTax_analysis_ensemble_comps2_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# plot the results:
p.comp2.silva <- ggplot(plt.all.comp.nona, aes(fill = tbl, x = variable, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Relative to ensemble-res", y = "Proportion of ASVs") + 
    scale_x_discrete(breaks = c("all.match","mis","over","under"), 
                   labels = c("Agree", "Misclassified", "Overclassified","Underclassified")) +
    scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                   limits = c(0, 1), expand = c(0,0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_fill_manual(name = "Taxonomy",
                      breaks = c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-acc"),
                      values = cbPalette[c(1:5)])
  
print(p.comp2.silva)
```

![](eTax_analysis_ensemble_comps2_files/figure-markdown_github/unnamed-chunk-10-2.png)

``` r
library("cowplot")
p.silva <- plot_grid(
  p.silva.2 + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  p.comp1.silva + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  p.comp2.silva + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  align = 'h',
  labels = c("(D)", "(E)", "(F)"),
  axis = 'l',
  hjust=-1,
  nrow=1
)
```

    ## Warning: Removed 33 rows containing missing values (geom_bar).

``` r
legend_b <- get_legend(p.pr2.2 + 
                         theme(legend.position="bottom",
                              legend.text=element_text(size=14),
                              legend.title = element_text(size=14),
                              legend.justification="center",
                              legend.box.margin = unit(c(0.1, 1, 0.1, 1),"cm")))
```

    ## Warning: Removed 45 rows containing missing values (geom_bar).

``` r
p.silva <- plot_grid(p.silva, legend_b, nrow=2, rel_heights=c(1, 0.3))
ggsave("all_silva_plots_protists", plot = p.silva, device = "pdf", width = 18, height = 7.5)

p.all.2 <- plot_grid(p.all, p.silva, nrow = 2, rel_heights=c(1, 1))
ggsave("enscomps_plots_pr2ANDsilva_protists", plot = p.all.2, device = "pdf", width = 18, height = 15)
```

re-arranging panels for better manuscript fig:

``` r
p.all <- plot_grid(
  p.pr2.2 + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"), title = element_text(face = "bold")) + ggtitle("pr2 taxonomic nomenclature") + theme(plot.title = element_text(hjust = 0.5)),
  p.comp1.pr2 + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  p.comp2.pr2 + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  align = 'h',
  labels = c("(A)", "(C)", "(E)"),
  axis = 'l',
  hjust=-1,
  nrow=3
)
```

    ## Warning: Removed 45 rows containing missing values (geom_bar).

``` r
p.silva <- plot_grid(
  p.silva.2 + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm"), title = element_text(face = "bold")) + ggtitle("silva taxonomic nomenclature") + theme(plot.title = element_text(hjust = 0.5)),
  p.comp1.silva + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  p.comp2.silva + theme(legend.position = "none", plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")),
  align = 'h',
  labels = c("(B)", "(D)", "(F)"),
  axis = 'l',
  hjust=-1,
  nrow=3
)
```

    ## Warning: Removed 33 rows containing missing values (geom_bar).

``` r
legend_b <- get_legend(p.pr2.2 + 
                         theme(legend.position="right",
                              legend.text=element_text(size=14),
                              legend.title = element_text(size=14, face = "bold"),
                              legend.justification="left",
                              legend.box.margin = unit(c(0.1, 1, 0.1, 1),"cm")))
```

    ## Warning: Removed 45 rows containing missing values (geom_bar).

``` r
p.all.2 <- plot_grid(p.all, p.silva, legend_b, nrow=1)
ggsave("enscomps_plots_pr2ANDsilva_protists_3x2", plot = p.all.2, device = "pdf", width = 11, height = 12)
```
