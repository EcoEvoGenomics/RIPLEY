# RIPLEY: Reproducible Pipelines for Genomic Analyses
## Introduction
RIPLEY is a repository of self-contained, reproducible Nextflow pipelines to run the many steps involved in common analyses automatically and in one terminal command. This *enormously* simplifies the process of e.g. producing a PCA and ADMIXTURE analysis from a VCF generated with our [`genotyping_pipeline`](https://github.com/EcoEvoGenomics/genotyping_pipeline), or a phylogeny of mitochondrial haplotypes from CRAMs.

## Instructions
> Instructions for use will follow at a later time ...

## Examples
### `metadata` (.csv)
```
Sample0,SpeciesA,PopA,M
Sample1,SpeciesA,PopB,F
Sample2,SpeciesA,PopC,F
Sample3,SpeciesA,PopD,M
Sample4,SpeciesA,PopE,M
Sample5,SpeciesB,PopF,F
Sample6,SpeciesB,PopG,M
Sample7,SpeciesB,PopH,F
Sample8,SpeciesB,PopI,F
Sample9,SpeciesB,PopJ,F
```

### `ref_chrom_labels` (.csv)
```
Original1,1,"Chromosome 1"
Original2,2,"Chromosome 2"
Original3,3,"Chromosome 3"
```
