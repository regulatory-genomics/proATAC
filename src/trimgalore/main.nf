nextflow.enable.dsl=2
include {
    TRIMGALORE__TRIM;
} from './processes/trim.nf'

include {
     SIMPLE_PUBLISH as PUBLISH_FASTQS_TRIMLOG_PE1;
     SIMPLE_PUBLISH as PUBLISH_FASTQS_TRIMLOG_PE2;
} from '../utils/processes/utils.nf'

workflow adapter_trimming {

    take:
        fastq_dex

    main:

        // run adapter trimming:
        switch(params.tools.trim.software) {
            case 'Trim_Galore':
                fastq_dex_trim = TRIMGALORE__TRIM(fastq_dex);
                PUBLISH_FASTQS_TRIMLOG_PE1(fastq_dex_trim.map{ it -> tuple(it[0], it[3]) }, '.R1.trimming_report.txt', 'reports/trim');
                PUBLISH_FASTQS_TRIMLOG_PE2(fastq_dex_trim.map{ it -> tuple(it[0], it[4]) }, '.R2.trimming_report.txt', 'reports/trim');
                break;
        }

    emit:
        fastq_dex_trim

}

