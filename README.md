# Loupe_BrainTF
Transcription factors (TFs) orchestrate gene expression programs crucial for cell physiology, but our knowledge of their function in the brain is limited. Using bulk tissues and sorted nuclei from multiple human post-mortem brain regions, we generated a multi-omic resource (1030 total experiments) that includes binding maps for >100 TFs. We demonstrate improved measurements of TF activity, including motif recognition and gene expression modeling, upon identification and removal of regions of high TF occupancy. Further, we find that predictive TF binding models demonstrate a bias for these high occupancy sites. Neuronal TFs SATB2 and TBR1 bind unique regions depleted for such sites and promote neuronal gene expression. Several TFs, including TBR1 and PKNOX1, are enriched for risk variants associated with neuropsychiatric disorders, predominantly in neurons. These data are a powerful resource for future studies seeking to understand the role of TFs in epigenetic regulation in the human brain.


<img src="https://github.com/aanderson54/Loupe_BrainTF/blob/main/images/Figure1.png" width="700" />

## Code Availability

#### This repository contains all code generated for Loupe et al. 2023 (click [here](https://aanderson54.github.io/Loupe_BrainTF/) to view code). Pre-processing scripts can be found in the scripts directory.



## Data Usage
In addition, we have included a large portion of the data here ([data/](https://aanderson54.com/Loupe_BrainTF/data)) to facilitate the use of this resource. We have included ChIP-seq peaks from experiments performed in the 4 large brain regions (CB, DLPFC, FP, and OL) and in all cell types (Bulk, NeuN+, Olig2+, Neg). Peaks ahave been formatted as `GenomicRanges` objects and saved as an `rds` for easy download. Each `unlist_*.rds` contains  a list of the peaks called for each TF for that brain region/cell-type. Some examples of data interaction are included below. To learn more about GenomicRanges object, see their [vignette](https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html).

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
```{r}
CTCF_peaks<-unlist_DLPFC[which(unlist_DLPFC$name=="CTCF"),]
```
