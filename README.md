# Multiomic profiling of transcription factor binding and function in human brain

### Manuscript: Now published in [Nature Neuroscience](https://doi.org/10.1038/s41593-024-01658-8) 

### Data availability
Data can be accessed on [synapse](https://doi.org/10.7303/syn51942384.1) or through our ChIP-portal. \
\
[![ChIP-portal](https://img.shields.io/badge/ChIP--portal-blue?style=for-the-badge)]([https://chip-portal.hudsonalpha.org/experiments])


### Abstract
Transcription factors (TFs) orchestrate gene expression programs crucial for cell physiology, but our knowledge of their function in the brain is limited. Using bulk tissues and sorted nuclei from multiple human post-mortem brain regions, we generated a multi-omic resource (1121 total experiments) that includes binding maps for more than 100 TFs. We demonstrate improved measurements of TF activity, including motif recognition and gene expression modeling, upon identification and removal of regions of high TF occupancy. Further, we find that predictive TF binding models demonstrate a bias for these high occupancy sites. Neuronal TFs SATB2 and TBR1 bind unique regions depleted for such sites and promote neuronal gene expression. Several TFs, including TBR1 and PKNOX1, are enriched for risk variants associated with neuropsychiatric disorders, predominantly in neurons. These data are a powerful resource for future studies seeking to understand the role of TFs in epigenetic regulation in the human brain.

<img src="https://github.com/aanderson54/Loupe_BrainTF/blob/main/images/Figure1_Page_1.png" width="700" />

## Code Availability

#### This repository contains all code generated for Loupe et al. 2024 (click [here](https://aanderson54.github.io/Loupe_BrainTF/) to view code). Pre-processing scripts can be found in the scripts directory.


[![DOI](https://zenodo.org/badge/635343701.svg)](https://zenodo.org/badge/latestdoi/635343701)





## Data Usage
All data can be accessed at [https://doi.org/10.7303/syn51942384.1](https://doi.org/10.7303/syn51942384.1)

In addition, we have included a large portion of the data here ([data/](https://github.com/aanderson54/Loupe_BrainTF/tree/main/data)) to facilitate the use of this resource. We have included ChIP-seq peaks from experiments performed in the 4 large brain regions (CB, DLPFC, FP, and OL) and in all cell types (Bulk, NeuN+, Olig2+, Neg). Peaks have been formatted as `GenomicRanges` objects and saved as an `rds` for easy download. Each `unlist_*.rds` contains  a list of the peaks called for each TF for that brain region/cell-type. Some examples of data interaction are included below. To learn more about GenomicRanges object, see their [vignette](https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html).

#### Load data
```{r}
library(GenomicRanges)
unlist_DLPFC<-readRDS("unlist_DLPFC.rds")
head(unlist_DLPFC)

      seqnames     ranges        strand |        name     score  signalValue pValue    qValue       peak
  [1]     chr1 16644628-16644842      * |       ASH2L      1000     132.671       -1   3.50934       104
  [2]     chr6 17706875-17707084      * |       ASH2L      1000     119.128       -1   3.50934       131
  [3]    chr18 50281464-50281927      * |       ASH2L      1000     119.083       -1   3.50934       232
  [4]    chr17 39980706-39980909      * |       ASH2L      1000     117.065       -1   3.50934        87
  [5]     chr6 26596803-26597010      * |       ASH2L      1000     114.682       -1   3.50934        87
  [6]    chr17 26619998-26620117      * |       ASH2L      1000     113.692       -1   3.50934        32

```
#### Subset Peaks
To look at only CTCF peaks:
```{r}
CTCF_peaks<-unlist_DLPFC[which(unlist_DLPFC$name=="CTCF"),]
```

#### Union Peaks

```{r}
union_DLPFC<-reduce(unlist_DLPFC)
```
Large union peaks were split into peaks of up to 2kb
```{bash}
#bash
Rscript Breaking_up_merged_peaks_max_2kb_bins_parallel.R union_DLPFC.rds unlist_DLPFC.rds
```

Once large peaks have been broken up, we count the number of TFs binding at a union peak and define HOT sites.
```{r}
# Get overlaps between union and unlist
overlaps<-findOverlaps(unlist_DLPFC, union_DLPFC)
mcols(union_DLPFC)$TF<-NA #set empty var to write into

# Get TF name for each unlist peak that overlaps
tmp<-mcols(unlist_DLPFC)$name[queryHits(overlaps)]  
df<-data.frame(tf=tmp, subjectHit=subjectHits(overlaps)) 
tidy<-df %>% group_by(subjectHit) %>% summarize(name=paste(sort(unique(tf)), collapse=",")) #collapse tf names on union index

# Count the number of TFs by counting commas +1
tidy$length<-str_count(tidy$name,",")+1 
mcols(union_DLPFC)$TF[tidy$subjectHit]<-tidy$name #add list of tfs to union grange
mcols(union_DLPFC)$count[tidy$subjectHit]<-tidy$length
names(union_DLPFC)<-NULL
```
#### Define HOT sites as greater than 90th percentile
```{r}
union_DLPFC$HOT<-ifelse(union_DLPFC$count> quantile(union_DLPFC$count, probs=0.9), "HOT","NOT")
```
