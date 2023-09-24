nextflow.enable.dsl=2
toolParams = params.tools
process SCTK__BARCODE_CORRECTION {
     
    label 'compute_resources__sctk_barcode'

     input:
       tuple val(sampleId),
             val(technology),
             path(fastq_PE1),
             path(fastq_bc),
             path(fastq_PE2),
             path(bc_whitelist)
 
     output:
       tuple val(sampleId),
             val(technology),
             path(fastq_PE1),
             path("${sampleId}_bc_corrected.fastq.gz"),
             path(fastq_PE2),
             path("${sampleId}_bc_corrected.fastq.gz.corrected_bc_stats.tsv")

     script:
       // def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_10x_scatac_fastqs)
       // processParams = sampleParams.local
        """
          correct_barcode_from_fastq.sh \
          ${bc_whitelist} \
          false \
          $fastq_bc \
          ${sampleId}_bc_corrected.fastq.gz

        """
}
