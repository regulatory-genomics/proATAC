nextflow.enable.dsl=2

process TRIMGALORE__TRIM {

    label 'compute_resources__trimgalore__trim'

    input:
        tuple val(sampleId),
              path(fastq_PE1),
              path(fastq_PE2)

    output:
        tuple val(sampleId),
              path("${sampleId}_dex_R1_val_1.fq.gz"),
              path("${sampleId}_dex_R2_val_2.fq.gz"),
              path("${sampleId}_dex_R1.fastq.gz_trimming_report.txt"),
              path("${sampleId}_dex_R2.fastq.gz_trimming_report.txt")

    script:
     //   def sampleParams = params.parseConfig(sampleId, params.global, toolParams.trim)
     //   processParams = sampleParams.local
        def max_threads = (task.cpus > 6) ? 6 : task.cpus
        """
        ${params.tools.trim.trim_galore_dir} \
            -j ${max_threads} \
            -o . \
            ${fastq_PE1} \
            ${fastq_PE2} \
            --paired \
            --gzip
        """
}
