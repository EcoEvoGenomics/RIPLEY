process BCFTOOLS_SAMPLE_VCF {

    label "BCFTOOLS"

    input:
    tuple path(vcf), path(csi)
    val(n_sites)

    output:
    path("${vcf.simpleName}.sample.vcf.gz")

    // Reservoir-sampling awk command courtesy of GPT UiO (GPT-5; https://gpt.uio.no/)

    script:
    """
    total_sites=\$(bcftools index -n ${vcf})
    bcftools query --format '%CHROM\\t%POS' ${vcf} > chrom_pos.txt
    
    awk '
    BEGIN {
      k = ${n_sites}; srand();
    }
    {
      if (NR <= k) {
        buf[NR] = \$0;
      } else {
        i = int(rand() * NR) + 1;
        if (i <= k) buf[i] = \$0;
      }
    }
    END {
      for (i = 1; i <= k && i in buf; i++) print buf[i];
    }
    ' chrom_pos.txt > chrom_pos_sampled.txt

    bcftools view -R chrom_pos_sampled.txt -O z -o ${vcf.simpleName}.sample.vcf.gz ${vcf}
    """
}
