nextflow.enable.dsl=2

log.info """
    FILLING FILE PIPELINE (DSL-2)
    ==============================

    ==============================
"""

// Include the modules
include { EXTRACT_TRANSCRIPTS_ENDS     } from './transsplicing_modules.nf'
include { COUNTING_KMERS_FASTK         } from './transsplicing_modules.nf'
include { COUNTING_KMERS_JELLYFISH     } from './transsplicing_modules.nf'
include { FORMATING_JELLYFISH_OUTPUT   } from './transsplicing_modules.nf'
include { GETTING_PUTATIVE_SLS         } from './transsplicing_modules.nf'
include { CLUSTERING_PUTATIVE_SLS      } from './transsplicing_modules.nf'
include { SUMMARIZING_CLUSTERS         } from './transsplicing_modules.nf'
include { HEATMAP                      } from './transsplicing_modules.nf'
include { HISTOGRAM_PLOTTING_JELLYFISH } from './transsplicing_modules.nf'
include { HISTOGRAM_PLOTTING           } from './transsplicing_modules.nf'

// Define the workflow
workflow {
    EXTRACT_TRANSCRIPTS_ENDS(params.transcripts, params.length)
    
    if (params.kmer_counting_tool == "fastk") {
        COUNTING_KMERS_FASTK(EXTRACT_TRANSCRIPTS_ENDS.out, params.k, params.lower_cov_thrsld_counting)
        HISTOGRAM_PLOTTING(COUNTING_KMERS_FASTK.out.kmers_fastk_histex)
        GETTING_PUTATIVE_SLS(EXTRACT_TRANSCRIPTS_ENDS.out
                            .combine(COUNTING_KMERS_FASTK.out.kmers_fastk_formatted)
                            .map { extracted, kmers -> tuple(params.extracted, kmers, params.lower_cov_thrsld_extracting) }, 
                         params.max_distance, 
                         params.extra_border, 
                         params.size_limit_max, 
                         params.size_limit_min, 
                         params.entropy_lim, 
                         params.nucleotide_limit)
    }

    if (params.kmer_counting_tool == "jellyfish") {
        COUNTING_KMERS_JELLYFISH(EXTRACT_TRANSCRIPTS_ENDS.out, params.k)
        FORMATING_JELLYFISH_OUTPUT(COUNTING_KMERS_JELLYFISH.out.kmers_jellyfish_jf, params.lower_cov_thrsld_counting)
        HISTOGRAM_PLOTTING(FORMATING_JELLYFISH_OUTPUT.out.kmers_jellyfish_histo)
        GETTING_PUTATIVE_SLS(EXTRACT_TRANSCRIPTS_ENDS.out,
                            FORMATING_JELLYFISH_OUTPUT.out.kmers_jellyfish_formatted,
                            params.lower_cov_thrsld_extracting, 
                            params.max_distance, 
                            params.extra_border, 
                            params.size_limit_max, 
                            params.size_limit_min, 
                            params.entropy_lim, 
                            params.nucleotide_limit)

        // HISTOGRAM_PLOTTING_JELLYFISH(FORMATING_JELLYFISH_OUTPUT.out.kmers_jellyfish_histo)
    }

    CLUSTERING_PUTATIVE_SLS(GETTING_PUTATIVE_SLS.out, params.identity_thrsld)
    SUMMARIZING_CLUSTERS(CLUSTERING_PUTATIVE_SLS.out.final_putative_sls, CLUSTERING_PUTATIVE_SLS.out.clustering_info)
    HEATMAP(SUMMARIZING_CLUSTERS.out)

}


// End of the workflow