nextflow.enable.dsl=2
// process imports:
include {SCTK__BARCODE_CORRECTION;} from './processes/barcode_correction.nf'
include {add_bc2read;} from './processes/add_barcode2read.nf'
include {
    SIMPLE_PUBLISH as PUBLISH_BC_STATS;
} from '../utils/processes/utils.nf'
//  Define the workflow
workflow barcode_correction {
     take:
        metadata
     main:
        data = Channel.from(metadata)
        .splitCsv(
        sep: '\t',
        header: true,
        strip: true
        ).map {
                          row -> tuple(
                              row.sample_name + "___" + file(row.fastq_PE1_path)
                                  .getSimpleName()
                                  .replaceAll(row.sample_name,""),
                              row.technology,
                              file(row.fastq_PE1_path, checkIfExists: true),
                              row.fastq_barcode_path,
                              file(row.fastq_PE2_path, checkIfExists: true)
                              )
               } 
      
 // gather barcode whitelists from params into a channel:
        wl = Channel.empty()
        wl_cnt = 0
        params.data.barcode_correction.whitelist.each {k,v ->  if(v!=''){
        wl = wl.mix( Channel.of(tuple(k, file(v)) ))
                wl_cnt = wl_cnt + 1
        }
                                                       }
        /* TO DO: fix ability to skip barcode correction */
        if(wl_cnt == 0) {
        //    if(!params.containsKey('quiet')) {
        //      println("No whitelist files were found in 'params.data.barcode_correction.whitelist'. Skipping barcode correction for standard-type samples.")
            
        // run barcode demultiplexing on each read+barcode:
            fastq_dex = add_bc2read(data)
        } else {
            // join wl to the data channel:
        data_wl = wl.cross( data.map { it -> tuple(it[1], it[0], it[2], it[3], it[4]) } ) // technology, sampleId, R1, R2, R3
                        .map { it -> tuple(it[1][1], it[1][0],           // sampleId, technology
                                           it[1][2], it[1][3], it[1][4], // R1, R2, R3
                                           it[0][1]                      // whitelist
                                           ) }

        // run barcode correction against a whitelist:
        fastq_bc_corrected = SCTK__BARCODE_CORRECTION(data_wl)
        PUBLISH_BC_STATS(fastq_bc_corrected.map { it -> tuple(it[0], it[5]) }, '.corrected_bc_stats.tsv', 'reports/barcode')

        // run barcode demultiplexing on each read+barcode:
        fastq_dex = add_bc2read(
                fastq_bc_corrected.map { it -> tuple(*it[0..4]) }
            )

}
emit:
        fastq_dex
}
