include { BCFTOOLS_CALL_REGION_VARIANTS; BCFTOOLS_NORMALISE } from "../process/bcftools.nf"
include { BCFTOOLS_INDEX; BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORMALISED} from "../process/bcftools.nf"
include { BCFTOOLS_MAKE_CONSENSUS_FASTA } from "../process/bcftools.nf"
include { SAMTOOLS_EXTRACT_FASTA; SAMTOOLS_INDEX_FASTA } from "../process/samtools.nf"
include { CONCATENATE_FILES as CONCATENATE_FASTAS } from "../process/system.nf"
include { MAFFT_ALIGN } from "../process/mafft.nf"
include { IQTREE_BUILD_TREE; IQTREE_TO_PLAIN_NEWICK } from "../process/iqtree.nf"
include { METADATA_TO_SPART } from "../process/system.nf"

nextflow.preview.output = true

workflow {
    main:
    crams = Channel.fromPath("${params.mt_cramdir}/**.cram")
    BCFTOOLS_CALL_REGION_VARIANTS(crams, params.metadata, params.ref_genome, params.ref_ploidy, params.mt_genome_region)
    BCFTOOLS_INDEX(BCFTOOLS_CALL_REGION_VARIANTS.out.vcf)
    BCFTOOLS_NORMALISE(BCFTOOLS_INDEX.out.indexed_vcf, params.ref_genome)
    BCFTOOLS_INDEX_NORMALISED(BCFTOOLS_NORMALISE.out.normalised_vcf)
    BCFTOOLS_MAKE_CONSENSUS_FASTA(BCFTOOLS_INDEX_NORMALISED.out.indexed_vcf, params.mt_filt_indelgap, params.mt_filt_inclusions, params.ref_genome)
    SAMTOOLS_INDEX_FASTA(BCFTOOLS_MAKE_CONSENSUS_FASTA.out.fasta)
    SAMTOOLS_EXTRACT_FASTA(SAMTOOLS_INDEX_FASTA.out.indexed_fasta, params.mt_genome_region)
    haplotype_fastas = SAMTOOLS_EXTRACT_FASTA.out.extracted.collect()

    CONCATENATE_FASTAS(haplotype_fastas, "${params.mt_genome_region}.fasta")
    MAFFT_ALIGN(CONCATENATE_FASTAS.out.concat, params.mt_mafft_iterations)
    IQTREE_BUILD_TREE(MAFFT_ALIGN.out.aligned, params.mt_iqtree_bootstraps)
    IQTREE_TO_PLAIN_NEWICK(IQTREE_BUILD_TREE.out.contree)

    METADATA_TO_SPART(params.metadata)

    publish:
    fasta = SAMTOOLS_EXTRACT_FASTA.out
    iqtree = IQTREE_BUILD_TREE.out.all_treefiles
    alignment = MAFFT_ALIGN.out.aligned
    plain_newick = IQTREE_TO_PLAIN_NEWICK.out
    spart = METADATA_TO_SPART.out
}

output {
    fasta { path "mitochondrial_haplotypes/fasta" }
    iqtree { path "mitochondrial_haplotypes/iqtree" }
    alignment { path "mitochondrial_haplotypes/hapsolutely" }
    plain_newick { path "mitochondrial_haplotypes/hapsolutely" }
    spart { path "mitochondrial_haplotypes/hapsolutely" }
}
