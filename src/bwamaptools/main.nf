
nextflow.enable.dsl=2
// include mapping stats
include {
    BWAMAPTOOLS__MAPPING_SUMMARY as MAPPING_SUMMARY;
} from './processes/mapping_summary.nf' params(params)
include {
    BWA_MEM_PE;
} from './processes/bwa_mapping.nf'
include {
    SIMPLE_PUBLISH as PUBLISH_BAM;
    SIMPLE_PUBLISH as PUBLISH_BAM_INDEX;
    SIMPLE_PUBLISH as PUBLISH_MAPPING_SUMMARY;
} from '../../src/utils/processes/utils.nf'

include {
    SAMTOOLS__MERGE_BAM;
} from '../../src/samtools/processes/merge_bam.nf'


import java.nio.file.Paths

workflow mapping {

    take:
        fastq_dex_trim

    main:

        // map with bwa mem:
        aligned_bam = BWA_MAPPING_PE(
            fastq_dex_trim.map { it -> tuple(it[0].split("___")[0], // [val(unique_sampleId),
                                             *it[0..2] ) // val(sampleId), path(fastq_PE1), path(fastq_PE2)]
                  })

        // split by sample size:
        aligned_bam.map{ it -> tuple(it[0].split("___")[0], it[1]) } // [ sampleId, bam ]
                   .groupTuple()
                   .branch {
                       to_merge: it[1].size() > 1
                       no_merge: it[1].size() == 1
                   }
                   .set { aligned_bam_size_split }

        // merge samples with multiple files:
        bam_merged = SAMTOOLS__MERGE_BAM(aligned_bam_size_split.to_merge)

        // re-combine with single files:
        bam_merged.mix(aligned_bam_size_split.no_merge.map { it -> tuple(it[0], *it[1]) })
           .set { bam }

        // publish merged BAM files or only BAM file per sample:
        PUBLISH_BAM(bam.map{ it -> tuple(it[0..1]) }, '.bwa.out.possorted.bam', 'bam')
        PUBLISH_BAM_INDEX(bam.map{ it -> tuple(it[0],it[2]) }, '.bwa.out.possorted.bam.bai', 'bam')



    emit:
        bam

}
       


/* sub-workflows used above */  
workflow get_bwa_index {

    
    take:
        fasta_path

    main:

        bwa_fasta = Channel.fromPath(fasta_path)

        bwa_index_path = Paths.get(
                                   Paths.get(fasta_path).getParent().toString(),
                                   "*.{amb,ann,pac,0123,bwt.2bit.64}"
                                   )
        bwa_index = Channel.fromPath(bwa_index_path,
                                     glob: true,
                                     type: 'file',
                                     )
                           .ifEmpty { exit 1, "ERROR: Could not find bwa indices from: ${bwa_index_path}." }
                           .collect()
                           .toList()

        data_channel = bwa_fasta.combine(bwa_index)

    emit:
        data_channel

}

workflow BWA_MAPPING_PE {

    take:
        data // a channel of [val(unique_sampleId), val(sampleId), path(fastq_PE1), path(fastq_PE2)]
        // unique_sampleId is used to label the read group field "SM" and (part of) "LB",
        // while sampleId represents each split fastq file for a unique sample.

    main:
        /*
           1) create a channel linking bwa index files from genome.fa in params, and
           2) combine this channel with the items in the data channel
        */
        bwa_inputs = get_bwa_index(params.tools.bwamaptools.bwa_fasta).combine(data)

        aligned_bam = BWA_MEM_PE(bwa_inputs)


        // publish output:

        MAPPING_SUMMARY(aligned_bam)
        PUBLISH_MAPPING_SUMMARY(MAPPING_SUMMARY.out, '.mapping_stats.tsv', 'reports/mapping_stats')

    emit:
        aligned_bam
}

    
