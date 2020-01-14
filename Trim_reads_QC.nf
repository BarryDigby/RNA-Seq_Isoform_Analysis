#!/usr/bin/env nextflow

params.reads = "/data/MA5112/Practicals/RNA-Seq/Stringtie_Practical/Data/*r{1,2}.fastq.gz"
Channel
        .fromFilePairs( params.reads )
        .set { input_fq }

process bbduk {
        publishDir "trimmed_reads/", mode: 'copy'

        input:
        set val(key), file(reads) from input_fq

        output:
        file "*trimmed.fastq.gz" into fastqc_input
        file "*.stats.txt"

        script:
        """
        bbduk \
        -Xmx4g \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out1=${key}_r1.trimmed.fastq.gz \
        out2=${key}_r2.trimmed.fastq.gz \
        literal=AGATCGGAAGAG \
        minlen=30 \
        ktrim=r \
        k=12 \
        qtrim=rl \
        trimq=20 \
        stats=${key}.stats.txt
        """
}

process fastqc {
        publishDir "fastqc/Post-Trim/", mode: 'copy'

        input:
        set val(pair_id), file(trimmed_reads) from fastqc_input

        output:
        file "*.{zip,html}" into fastqc_output

        script:
        """
        fastqc -t 8 -q ${trimmed_reads}
        """
}

process multiqc {
        publishDir "fastqc/Post-Trim/", mode:'copy'

        input:
        file("*") from fastqc_output.collect().ifEmpty([])

        output:
        file "multiqc_report.html" into multiqc_report
        file "multiqc_data"

        script:
        """
        multiqc .
        """
}
