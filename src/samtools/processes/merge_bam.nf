nextflow.enable.dsl=2



process SAMTOOLS__MERGE_BAM {
    label 'compute_resources__samtools__merge_bam'

    input:
        tuple val(sampleId),
              path(bams)

    output:
        tuple val(sampleId),
              path("${sampleId}.bwa.out.fixmate.possorted.merged.bam"),
              path("${sampleId}.bwa.out.fixmate.possorted.merged.bam.bai")

    script:
        //def sampleParams = params.parseConfig(sampleId, params.global)
        //processParams = sampleParams.local
        """
        set -euo pipefail

        samtools merge \
            -@ 4 \
            -O bam \
            --write-index \
            -o '${sampleId}.bwa.out.fixmate.possorted.merged.bam##idx##${sampleId}.bwa.out.fixmate.possorted.merged.bam.bai' \
            ${"'" + bams.join("' '") + "'"}
        """

}
