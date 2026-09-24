# RIPLEY: Reproducible Pipelines for Genomic Analyses

*Note: This is an early release and the current documentation is a draft.*

## What is RIPLEY?
RIPLEY is a collection of reproducible analysis pipelines built for population genomics. Following and extending ["Speciation & Population Genomics: a how-to-guide"](https://speciationgenomics.github.io/), RIPLEY aims to make common analyses both accessible and reproducible - all in a single terminal command. The companion repositories [XENO](https://github.com/EcoEvoGenomics/XENO) and [MORPH](https://github.com/EcoEvoGenomics/MORPH) are designed to expedite the process from raw reads to RIPLEY-ready input.

## Installation
### Prerequesites

*In short: Nextflow 25.04.4 / 25.04.6, Apptainer, Conda, Internet access.*

The only hard dependency is [Nextflow](https://www.nextflow.io/) version 25.04 (tested with 25.04.4 and 25.04.6). You also need a container engine, like Apptainer or Podman, and for most modules you need Conda. All other software requirements are handled by RIPLEY, which means internet access is the final requirement. In principle most RIPLEY modules *can* run on a laptop, but they were all *designed* to run on a high-performance computing (HPC) cluster with Slurm. Also refer to: [site configuration](#site-configuration).

### Installing RIPLEY
There is no installation other than cloning:
```sh
git clone https://github.com/EcoEvoGenomics/RIPLEY
```
Working files will be stored in this directory, so be mindful of scarce or shared storage.


## Quickstart
### Load dependencies
This step inherently depends on your environment and [site configuration](#site-configuration). For illustration, you might have to call something alike the following:
```sh
module load singularity
module load conda
conda activate Nextflow25.04.6
```

### The RIPLEY launcher
To launch RIPLEY, you invoke the RIPLEY launcher. Every invocation requires exactly three ordered arguments:
```sh
./RIPLEY <module> <options.yaml> <site.config>
```
- `<module>` is the name of an available [module](#modules). Call `./RIPLEY` with no arguments to list the available modules.
- `<options.yaml>` is a `.yaml` file with all mandatory [common](#options-common-to-every-module) and [module-specific](#modules) options, plus any optional ones (see [here](#options-file)).
- `<site.config>` is a site configuration (see [here](#site-configuration)).

All three are mandatory and none has a default. Paths are resolved relative to the directory you launch from, so you may keep your options file and site configuration anywhere. Interrupted runs continue where they stopped, but runs can be long and it is wise to launch from a `screen` terminal.

### Options file
The options file is a collection of key-value pairs in `.yaml` format:
```yml
option_a: value
option_b: value
```
You will find an exhaustive list of [common](#common-options) and [module-specific](#modules) options below.

### Site configuration
The site config tells Nextflow where and how to run tasks. Two templates are provided in `configs/`:

| Config | Use |
|--------|-----|
| `configs/local.config` | Runs every task on the local machine with Podman. |
| `configs/saga.config` | Configured for the NRIS Saga HPC. Runs tasks as Slurm jobs with Singularity. Set `process.clusterOptions` to your own project account. |

Copy one of these and adapt it to your own system, or write a config from scratch. You should consult the documentation for your HPC of choice in order to configure this file. The Nextflow [executor](https://docs.seqera.io/nextflow/executor), [queue](https://docs.seqera.io/nextflow/reference/process/directives/queue), and [container](https://docs.seqera.io/nextflow/container) documentation is helpful reading. Do not configure other options in the site config: it is for executor, queue and container engine settings only.

## Common options

### Reference
| Option | Description | Required    |
|--------|-------------|-------------|
| `ref_genome` | Path to the reference `.fasta`. A matching `.fai` and `.gff` must sit beside it under the same base name. | Yes |
| `ref_ploidy` | Path to a [ploidy file](https://samtools.github.io/bcftools/bcftools.html#ploidy). Its sex codes must match column four of `sample_metadata` (see [Metadata](#metadata)). | Yes |
| `ref_chrom_labels` | Path to a headerless `.csv`: contig names, short plot label, long plot label. (e.g. `chr1,1,'Chromosome 1'`) | Yes |
| `ref_exclude_prefix` | List of prefixes marking unattached scaffolds, e.g. `[ NW ]`. Unattached scaffols (prefixed contigs) are never analysed. | No |
| `ref_exclude_chroms` | List of chromosome names to exclude from all inputs, e.g. `[ chrW, mtDNA ]`. | No |
| `ref_exclude_coords` | Path to a headerless `.bed`: regions to exclude from all inputs. Comment (`#`) lines are permitted, a column header row is not. | No |

Contigs that survive prefix and name filtering are the analysis set. If the input has no variants on a retained chromosome, RIPLEY stops with an error. **Note:** RIPLEY enforces *strict* alphanumericity for contig names, sample names, and metadata.

### Metadata
| Option | Description | Required    |
|--------|-------------|-------------|
| `sample_metadata` | Path to a headerless `.csv`: sample, species, population, sex. | Yes |
| `population_metadata` | Path to a headerless `.csv`: population, hex colour for plots. | Yes |
| `species_metadata` | Path to a headerless `.csv`: species, hex colour for plots. | Yes |
| `focal_populations` | Optional list of focal populations for certain analyses, e.g. `[ PopA, PopB ]`. If omitted, every population in `sample_metadata` is focal. | No |

Again: *all* contig names, sample names, species names, population names, and sex codes must be strictly alphanumeric. Every population and species named in `sample_metadata` *must* have a complete entry in `population_metadata` and `species_metadata`. Every sex code in `sample_metadata` *must* match an entry in the ploidy file. Every sample in the input VCF or CRAM file(s) *must* appear in `sample_metadata` and vice versa. Rows in `population_metadata` are read in order and that order is the order in which populations appear in plots. Note, however, that plots which cluster samples hierarchically always follow the clustering.

These files are best written in a plain text editor: remember to include a final newline and no headers or trailing spaces.

## Modules and module-specific options

All module-specific options are required (`true`/`false` flags evaluate to false if unset). Greater technical detail about each module is found at the Wiki (later).

### Module: qc_alignments
Alignment quality control: `samtools stats` metrics and binned coverage (`samtools bedcov`), plotted across all samples and within each [focal population](#common-options).

| Option | Description | Example |
|--------|-------------|---------|
| `cram` | Path to a directory of two or more `.cram` files with names `<sample>.cram`. Recurses into subdirectories. All filenames must be unique. **Not** a glob. | `/data/alignments` |
| `coverage_binsize` | Bin size in base pairs for calculation of coverage. | `100000` |

### Module: qc_variants
Variant quality control: depth, missingness, site quality, allele frequency, Hardy-Weinberg and inbreeding statistics plus SNP density, plotted across all samples and within each [`focal population`](#common-options).

| Option | Description | Example |
|--------|-------------|---------|
| `vcf` | Path to a single `<alphanumeric>.vcf.gz`, or a directory holding one `<chrom>.vcf.gz` per chromosome. **Not** a glob. | `/data/variants` |
| `qc_thinning_target` | Number of SNPs randomly sampled across the `.vcf.gz` for calculation of statistics. Lower numbers mean faster compute times. Higher numbers mean greater accuracy. Record density is not thinned. | `100000` |
| `qc_snpden_binsize` | Bin size in base pairs for calculation of record density (i.e. SNP density). | `100000` |

### Module: filter_variants
Applies `vcftools` filters chromosome-wise and concatenates the result into one whole-genome VCF. Both the per-chromosome files and the concatenated VCF are published. Optionally applies filters within-population rather than at study-wide level.

| Option | Description | Example |
|--------|-------------|---------|
| `vcf` | Path to a single `<alphanumeric>.vcf.gz`, or a directory holding one `<chrom>.vcf.gz` per chromosome. **Not** a glob. | `/data/variants` |
| `flags` | Path to a `.txt` of [`vcftools` flags](https://vcftools.github.io/man_latest.html#SITE%20FILTERING%20OPTIONS), one per line. The filename minus its extension must be strictly alphanumeric, e.g. `biallelic.txt`. | `/data/filters/biallelic.txt` |
| `popwise` | If `true`, apply the filters within each [`focal population`](#common-options) rather than across all samples. | `true` |

With `popwise: true`, frequency-dependent flags such as `--maf`, `--hwe` and `--max-missing` become per-population thresholds, the retained sites are the union across populations when merged, and samples outside the [`focal population`](#common-options)s are dropped.

### Module: population_structure
Kinship, pairwise FST, PCA and ADMIXTURE. PCA and ADMIXTURE run on an LD-pruned variant set. Parental populations identified from ADMIXTURE results used to find ancestry-informative markers (AIMs) among *unpruned* input. Heterozygosity and hybrid index calculated on AIMs.

| Option | Description | Example |
|--------|-------------|---------|
| `vcf` | Path to a single `<alphanumeric>.vcf.gz`, or a directory holding one `<chrom>.vcf.gz` per chromosome. **Not** a glob. | `/data/filtered/biallelic.vcf.gz` |
| `prune_window_kb` | LD pruning window size in kilobases. | `50` |
| `prune_step_snps` | LD pruning step in variants. | `10` |
| `prune_threshold` | LD pruning r² threshold. | `0.1` |
| `admixture_kmin` | Lowest K to run. | `2` |
| `admixture_kmax` | Highest K to run. | `8` |
| `aim_parental_threshold` | Minimum ADMIXTURE assignment for a sample to represent a parental population. Range: (0.5, 1.0]. | `0.9` |
| `aim_variance_threshold` | Minimum between-population allele frequency variance for a variant to count as ancestry-informative. Range: (0, 0.5). | `0.2` |

## Third-party software
Thank you for using RIPLEY. We kindly encourage you to cite the third-party software relevant to your use:
- *Work in progress: please see configuration files under each source/\*/nextflow.config*

______
RIPLEY v. 0.0.1a | 2026 | Erik Sandertun Røed | https://github.com/EcoEvoGenomics/RIPLEY
