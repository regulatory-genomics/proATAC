nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/bwamaptools/bin/" : ""

toolParams = params.tools.bwamaptools

process BWAMAPTOOLS__MAPPING_SUMMARY {

    label 'compute_resources__default'

    input:
        tuple val(sampleId),
              path(bam),
              path(bai)

    output:
        tuple val(sampleId),
              path("${sampleId}.mapping_stats.tsv")

    script:
       // def sampleParams = params.parseConfig(sampleId, params.global, toolParams)
       // processParams = sampleParams.local
        """
        ${binDir}mapping_summary.sh \
            ${sampleId} \
            ${bam} \
        """
}

