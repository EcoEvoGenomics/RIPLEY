# RIPLEY: Reproducible Pipelines for Genomic Analyses
## Introduction
RIPLEY is a repository of self-contained, reproducible Nextflow pipelines to run the many steps involved in common analyses automatically and in one terminal command. This *enormously* simplifies the process of e.g. producing a PCA and ADMIXTURE analysis from a VCF generated with our [`genotyping_pipeline`](https://github.com/EcoEvoGenomics/genotyping_pipeline), or a phylogeny of mitochondrial haplotypes from CRAMs.

## Instructions
> Instructions for use will follow at a later time ...

## Examples
### `sample_metadata` (.csv)
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

### `population_metadata` (.csv)
Rows are taken in order, so the order of populations here is the order in which they appear in plots and legends. Plots which cluster samples hierarchically always follow the clustering.
```
PopA,#1B9E77
PopB,#D95F02
PopC,#7570B3
PopD,#E7298A
PopE,#66A61E
PopF,#E6AB02
PopG,#A6761D
PopH,#666666
PopI,#1F78B4
PopJ,#B2DF8A
```

### `species_metadata` (.csv)
```
SpeciesA,#4C72B0
SpeciesB,#DD8452
```

### `ref_chrom_labels` (.csv)
```
Original1,1,"Chromosome 1"
Original2,2,"Chromosome 2"
Original3,3,"Chromosome 3"
```

### `fv_filter_flags` (.txt)
Any `vcftools` filtering flags, one per line. The filename (minus its extension) names the concatenated VCF and so must be strictly alphanumeric, e.g. `default.txt`.
```
--min-alleles 2
--max-alleles 2
--max-missing 0.8
--min-meanDP 5
--max-meanDP 30
--minDP 5
--maxDP 30
--minQ 30
--mac 1
--hwe 0
--remove-filtered-all
```
