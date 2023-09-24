nextflow.enable.dsl=2

import java.nio.file.Paths
import nextflow.config.ConfigParser
import static groovy.json.JsonOutput.*
def isParamNull(param) {
    return param == null || param == "NULL"
}
def getPublishDir = { outDir, toolName ->
    if(isParamNull(toolName))
        return "${outDir}/data"
    return "${outDir}/data/${toolName.toLowerCase()}"
}

process SIMPLE_PUBLISH {

    publishDir \
        "${getPublishDir(params.global.outdir,toolName)}", \
        mode: "${params.utils.publish?.mode ? params.utils.publish.mode: 'link'}", \
        saveAs: { filename -> "${outputFileName}" }

    label 'compute_resources__minimal'

    input:
        tuple \
            val(sampleId), \
            path(f, stageAs: 'input_file')
        val(outputFileSuffix)
        val(toolName)

    output:
        tuple \
            val(sampleId), \
            path(outputFileName)

    script:
        outputFileName = "${sampleId}${outputFileSuffix}"
        """
        if [ ! -f ${outputFileName} ]; then
            ln -s input_file "${outputFileName}"
        fi
        """
}

