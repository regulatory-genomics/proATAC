nextflow.enable.dsl=2
toolParams = params.tools
process add_bc2read {

     label 'compute_resources__sctk__barcode_add'
     input:
       tuple val(sampleId),
             val(technology),
             path(fastq_PE1),
             path(fastq_bc),
             path(fastq_PE2)
 
     output:
       tuple val(sampleId),
             path("${sampleId}_dex_R1.fastq.gz"),
             path("${sampleId}_dex_R2.fastq.gz")

     script:
      //  def sampleParams = params.parseConfig(sampleId, params.global, toolParams.barcode_10x_scatac_fastqs)
      //  processParams = sampleParams.local ##bu dong
        """
          barcode_10x_scatac_fastqs.sh \
          $fastq_PE1 \
          $fastq_bc \
          $fastq_PE2 \
          ${sampleId}_dex \
          false \
          true \
          ${toolParams.barcode_10x_scatac_fastqs.uncorrected_bc_tag}_${toolParams.barcode_10x_scatac_fastqs.barcode_quality_tag}

        """
}
