process GET_GENOMICS_GENERAL {

    // See https://simonmartinlab.org/software/
    // Replaces "NaN" with "nan" for compatibility with recent Numpy versions

    label "SYSTEM"

    output:
    path("genomics_general-0.5/*")

    script:
    """
    curl -L https://github.com/simonhmartin/genomics_general/archive/refs/tags/v0.5.tar.gz > genomics_general.tar.gz
    tar -zxvf genomics_general.tar.gz
    sed -i -e 's/np.NaN/np.nan/g' genomics_general-0.5/genomics.py
    """
}

process GENOMICS_GENERAL_VCF_TO_GENO {

    label "NUMPY"

    input:
    path(genomics_general)
    path(vcf)

    output:
    path("${vcf.simpleName}.geno.gz")

    script:
    """
    python VCF_processing/parseVCF.py -i ${vcf} -o ${vcf.simpleName}.geno.gz
    """
}

process GENOMICS_GENERAL_POPGEN_WINDOWS {

    // NB! Always assumes input data is phased.

    label "NUMPY"

    input:
    path(genomics_general)
    tuple path(geno), path(target_sample_list), path(metadata)
    val(window_size)
    val(step_size)
    val(min_sites)

    output:
    path("${geno.simpleName}.csv")

    script:
    """
    awk -F, '
      FILENAME=="${target_sample_list}" {samples_in_geno[\$0]=1;next}
      samples_in_geno[\$1] {print \$1, \$3}
    ' ${target_sample_list} ${metadata} > sample.pops

    echo '-w ${window_size}' >> popgenWindows.args
    echo '-s ${step_size}' >> popgenWindows.args
    echo '-m ${min_sites}' >> popgenWindows.args
    echo '-g ${geno}' >> popgenWindows.args
    echo '-o ${geno.simpleName}.csv' >> popgenWindows.args
    echo '-T ${task.cpus}' >> popgenWindows.args
    echo '-f phased' >> popgenWindows.args
    cat sample.pops | awk '{print \$2}' | sed 's/^/-p /' | sort | uniq >> popgenWindows.args
    echo '--popsFile sample.pops' >> popgenWindows.args

    cat popgenWindows.args | xargs python popgenWindows.py
    """
}
