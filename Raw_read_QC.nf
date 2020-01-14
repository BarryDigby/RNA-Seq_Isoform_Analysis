#!/usr/bin/env nextflow

params.reads = "/data/MA5112/Practicals/RNA-Seq/Stringtie_Practical/Data/*{1,2}.fastq.gz"
Channel
        .fromFilePairs( params.reads )
        .set { read_files_fastqc }

process fastqc {
        publishDir "fastqc/Pre-Trim/", mode: 'copy'

        input:
        set val(name), file(reads) from read_files_fastqc

        output:
        file "*_fastqc.{zip,html}" into fastqc_results

        script:
        """
        fastqc -q $reads
        """
}

process multiqc {
        publishDir "fastqc/Pre-Trim/", mode: 'copy'

        input:
        file ('*') from fastqc_results.collect().ifEmpty([])

        output:
        file "multiqc_report.html" into multiqc_report
        file "multiqc_data"

        script:
        """
        multiqc .
        """
}
         
