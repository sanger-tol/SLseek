// Modules

process EXTRACT_TRANSCRIPTS_ENDS {
    publishDir "${params.outdir}/1extracted_transcripts_ends"

    input:
    path transcripts
    val length // length of the extracted ends 

    output: 
    path "${transcripts.getBaseName()}_extracted${length}.fa", emit: extracted_transcripts

    script:
    def output_filename = "${transcripts.getBaseName()}_extracted${length}.fa"
    """
    extract_initial_nucleotides.py -f ${transcripts} -n ${length} -o ${output_filename}
    """
}

process COUNTING_KMERS_FASTK {
    publishDir "${params.outdir}/kmer_counts"

    input: 
    path extracted_transcripts
    val k
    val lower_cov_thrsld_counting

    output: 
    path "fastk_${extracted_transcripts.getBaseName()}_k${k}.ktab", emit: kmers_fastk_ktab
    path "fastk_${extracted_transcripts.getBaseName()}_k${k}.hist", emit: kmers_fastk_hist
    path "fastk_${extracted_transcripts.getBaseName()}_k${k}_LC${lower_cov_thrsld_counting}.txt", emit: kmers_fastk_dump
    path "fastk_${extracted_transcripts.getBaseName()}_k${k}_LC${lower_cov_thrsld_counting}.tsv", emit: kmers_fastk_formatted
    path "fastk_${extracted_transcripts.getBaseName()}_k${k}_LC${lower_cov_thrsld_counting}.histex", emit: kmers_fastk_histex


    script:
    """
    FastK -t -k${k} ${extracted_transcripts} -N./fastk_${extracted_transcripts.getBaseName()}_k${k}
    Histex -kA fastk_${extracted_transcripts.getBaseName()}_k${k}.hist > fastk_${extracted_transcripts.getBaseName()}_k${k}_LC${lower_cov_thrsld_counting}.histex
    Tabex -t${lower_cov_thrsld_counting} fastk_${extracted_transcripts.getBaseName()}_k${k} LIST > fastk_${extracted_transcripts.getBaseName()}_k${k}_LC${lower_cov_thrsld_counting}.txt
    awk '/:/ && /=/ { gsub(":", "", \$1); gsub("=", "", \$3); \$2 = toupper(\$2); print \$2 "\\t" \$4 }' fastk_${extracted_transcripts.getBaseName()}_k${k}_LC${lower_cov_thrsld_counting}.txt > fastk_${extracted_transcripts.getBaseName()}_k${k}_LC${lower_cov_thrsld_counting}.tsv
    """
}


process COUNTING_KMERS_JELLYFISH {
    publishDir "${params.outdir}/2kmer_counts"

    input:
    path extracted_transcripts
    val k 

    output:
    path "jf_${extracted_transcripts.getBaseName()}_k${k}.jf", emit: kmers_jellyfish_jf


    script:
    """
    jellyfish count -C -m ${k} -s 100M ${extracted_transcripts} -o jf_${extracted_transcripts.getBaseName()}_k${k}.jf
    """
}


process FORMATING_JELLYFISH_OUTPUT {
    publishDir "${params.outdir}/2kmer_counts"

    input:
    path kmers_jellyfish_jf
    val lower_cov_thrsld_counting

    output:
    path "${kmers_jellyfish_jf.getBaseName()}_LC${lower_cov_thrsld_counting}.fa", emit: kmers_jellyfish_dump
    path "${kmers_jellyfish_jf.getBaseName()}_LC${lower_cov_thrsld_counting}.histo", emit: kmers_jellyfish_histo
    path "${kmers_jellyfish_jf.getBaseName()}_LC${lower_cov_thrsld_counting}.tsv", emit: kmers_jellyfish_formatted

    script:
    """
    jellyfish dump --lower-count=${lower_cov_thrsld_counting} ${kmers_jellyfish_jf} > ${kmers_jellyfish_jf.getBaseName()}_LC${lower_cov_thrsld_counting}.fa
    jellyfish histo ${kmers_jellyfish_jf} > ${kmers_jellyfish_jf.getBaseName()}_LC${lower_cov_thrsld_counting}.histo
    awk '/^>/ {count=substr(\$0,2); next} {print \$0 "\\t" count}' ${kmers_jellyfish_jf.getBaseName()}_LC${lower_cov_thrsld_counting}.fa > ${kmers_jellyfish_jf.getBaseName()}_LC${lower_cov_thrsld_counting}.tsv
    """
}

// process HISTOGRAM_PLOTTING_JELLYFISH {}

process GETTING_PUTATIVE_SLS {
    publishDir "${params.outdir}/3draft_putative_sls"
    
    input:
    path extracted_transcripts
    path kmers_table
    val abund_thrd 
    val max_distance
    val extra_border
    val size_limit_max
    val size_limit_min
    val entropy_lim
    val nucleotide_limit

    output:
    path "${kmers_table.getBaseName()}_minCov${abund_thrd}.fasta", emit: draft_putative_sls

    script:
    def output_filename = "${kmers_table.getBaseName()}_minCov${abund_thrd}"
    """
    getting_putative_sl.py -t ${extracted_transcripts} -k ${kmers_table} -a ${abund_thrd} -m ${max_distance} -b ${extra_border} -y ${size_limit_max} -x ${size_limit_min} -e ${entropy_lim} -n ${nucleotide_limit} -o ${output_filename} 
    """
}

process CLUSTERING_PUTATIVE_SLS {
    publishDir "${params.outdir}/4final_putative_sls"    

    input:
    path draft_putative_sls
    val identity_thrsld

    output:
    path "${draft_putative_sls.getBaseName()}_clustered_nr${identity_thrsld}", emit: final_putative_sls
    path "${draft_putative_sls.getBaseName()}_clustered_nr${identity_thrsld}.clstr", emit: clustering_info


    script:
    def output_filename = "${draft_putative_sls.getBaseName()}_clustered_nr${identity_thrsld}"
    def memory_mb = task.memory.toMega().intValue()
    """
    # Getting word size
    word_size=\$(awk -v id="${identity_thrsld}" 'BEGIN{
        if(id>0.95) print 10;
        else if(id>0.90) print 8;
        else if(id>0.88) print 7;
        else if(id>0.85) print 6;
        else if(id>0.80) print 5;
        else print 4;
    }')

    # Running cd-hit-est
    cd-hit-est -i ${draft_putative_sls} -o ${output_filename} -c ${identity_thrsld} -n \${word_size} -d 0 -M ${memory_mb} -T ${task.cpus}
    """

}

process SUMMARIZING_CLUSTERS {
    publishDir "${params.outdir}/5main_results"

    input:
    path final_putative_sls
    path clustering_info

    output:
    path "${final_putative_sls.getBaseName()}_clustered_putativeSLs.tsv", emit: cluster_summary

    script:
    def output_filename = "${final_putative_sls.getBaseName()}"
    """
    cluster_analysis.py -p ${final_putative_sls} -c ${clustering_info} -o ${output_filename}
    """
}

process HEATMAP {
    publishDir "${params.outdir}/6extra_results"

    input:
    path cluster_summary

    when:
    params.heatmap

    output:
    path "${cluster_summary.getBaseName()}_heatmap.png", emit: heatmap

    script:
    def output_filename = "${cluster_summary.getBaseName()}"
    """
    heatmap.py --tsv ${cluster_summary} --output ${output_filename}
    """
}

process HISTOGRAM_PLOTTING_JELLYFISH {
    publishDir "${params.outdir}/6extra_results"

    input:
    path kmers_jellyfish_histo

    when:
    params.kmerplot

    output:
    path "${kmers_jellyfish_histo.getBaseName()}_histogram.png", emit: histogram_jelyfish

    script:
    def output_filename = "${kmers_jellyfish_histo.getBaseName()}"
    """
    histogram.py --tsv ${kmers_jellyfish_histo} --output ${output_filename} --sep ' '
    """
}

process HISTOGRAM_PLOTTING {
    publishDir "${params.outdir}/6extra_results"

    input:
    path kmers_histo

    when:
    params.kmerplot

    output:
    path "${kmers_histo.getBaseName()}_histogram.png", emit: histogram

    script:
    def output_filename = "${kmers_histo.getBaseName()}"
    """
    histogram.py --tsv ${kmers_histo} --output ${output_filename}
    """
}