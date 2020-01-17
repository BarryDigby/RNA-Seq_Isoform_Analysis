#!/usr/bin/env nextflow

/*
 * Pipeline for RNA-seq analysis using Stringtie
 * for genome guided transcriptome assembly and quantification.
 */

/*
 * STEP 0: Set path to variables
 * Genome file will only be used once, we can read in as "file"
 * GTF file will be used multiple times. Split over n channels 
 * Fastq files are PE -- read in as pairs. 
 */ 

params.genome = "/data/MA5112/Practicals/RNA-Seq/Stringtie_Practical/Reference/chr22.fa"
genome_fasta = files( params.genome )

params.annot = "/data/MA5112/Practicals/RNA-Seq/Stringtie_Practical/Annotation/chr22.gtf"
Channel
	.fromPath( params.annot )
	.into { gtf1; gtf2; gtf3; gtf4 }

params.reads = "/data/MA5112/Practicals/RNA-Seq/Stringtie_Practical/trimmed_reads/*_r{1,2}.trimmed.fastq.gz"
Channel
	.fromFilePairs( params.reads )
	.set { read_ch }

/* 
 * STEP 1: Extract reference genome Exon/Splice sites
 * INPUTS:
 * Call GTF file from one of GTF channels defined in STEP 0. 
 * OUTPUTS:
 * Two files, place directly into channel.
 * chr22.ss will be used during 2 processes, place into two channels. 
 */

process Extract_Splice_Sites {
	publishDir "Reference/", mode:'copy'
	
	input:
	file(gtf) from gtf1

	output: 
	file "chr22.ss" into splice_sites, splice_sites1
	file "chr22.exon" into exon_sites

	script:
	"""
	hisat2_extract_splice_sites.py ${gtf} > chr22.ss
	hisat2_extract_exons.py ${gtf} > chr22.exon
	"""
}

/*
 * STEP 2: Index reference genome
 * INPUTS: 
 * Call outputs of STEP 1 & reference genome
 * OUTPUTS:
 * Collects Hisat2 index files (8 files ending in "ht2")
 * by using the "*" wildcard. 
 */

process Genome_Index { 
	publishDir "Reference/", mode:'copy'

	input:
	file ss from splice_sites
	file exon from exon_sites
	file fasta from genome_fasta

	output:
	file "${fasta.baseName}.*.ht2" into hs2_indices

	script:
	"""
	hisat2-build -p 8 --ss ${ss} --exon ${exon} ${fasta} ${fasta.baseName}.hisat2_index
	"""
}	

/* 
 * STEP 3: Align the trimmed reads to the reference genome
 * INPUTS:
 * Trimmed fastq files, chr22.ss from STEP 1,
 * genome index files from STEP 2 (multiple files put into channel, 
 * use collect() to capture all files vs. FIFO). 
 * OUTPUTS:
 * Summary txt file, bam files. Pipe output to samtools (Sam to Bam)
 * BASH:
 * To define a bash variable in nextflow you must export 
 * the variable. When calling it, escape the $ using \$. 
 */


process Hisat2_Alignment {
	publishDir "BAMS/", mode:'copy'

	input:
	set val(key), file(reads) from read_ch
	file hisat_index from hs2_indices.collect()
	file ss from splice_sites1.collect()
	
	output:
	file "${key}.summary.txt" into alignment_stats
	set val(key), file("${key}.bam") into hisat_bams
        
	script:
	index_base = hisat_index[0].toString() - ~/.\d.ht2l?/
        """
	export SM=`cut -d'_' -f2 <<< ${key}`
	export LB=`cut -d'_' -f1,2,3 <<< ${key}`
	export ID=`cut -d'_' -f2,3 <<< ${key}`
	export PU=`zcat < ${reads[0]} | head -n1 | awk 'BEGIN{FS=":"; OFS"."} {print \$3"_"\$4"_"\$10}'`
	export PL="ILLUMINA"	
	
	echo ${index_base}
	
	hisat2 -p 8 \
	--rg-id=\$ID \
	--rg SM:\$SM \
	--rg LB:\$LB \
	--rg PL:\$PL \
	--rg PU:\$PU \
	--known-splicesite-infile $ss \
	--no-mixed \
	--no-discordant \
	-x $index_base \
	--dta \
	--rna-strandness RF \
	-1 ${reads[0]} \
	-2 ${reads[1]} \
	--summary-file ${key}.summary.txt | samtools view -Sbh - > ${key}.bam
	"""
}

/* 
 * STEP 4: Sort and Index Bam files
 * INPUT: 
 * bams from STEP 3. 
 * OUTPUT: 
 * sorted bams into channel 
 */

process Sort_Index_Bams {
	publishDir "BAMS/", mode:'copy'

	input:
	set val(key), file(bam) from hisat_bams

	output:
	set val(key), file("${key}.bam") into hisat_bams1
	file "${key}.bam.bai" into indexed

	script:
	def avail_mem = task.memory == null ? '' : "-m ${task.memory.toBytes() / task.cpus}"
	"""
	samtools sort \\
	$bam \\
	-@ ${task.cpus} $avail_mem \\
	-o ${key}.bam
	samtools index ${key}.bam
	"""
}

// This command places the sorted bams into two channels

hisat_bams1.into { hisat_bams2; hisat_bams3 }

/* 
 * STEP 5: Assemble Transcripts using reference GTF as guide. 
 * INPUT: 
 * bams from STEP 4 and GTF from STEP 0. 
 * must use .combine() here so the GTF file
 * will be used with each incoming bam file. 
 * (if you do not it will run the script on
 * a random bam file once and end)
 * OUTPUT:
 * collect the GTF files
 */

process Assemble_Transcripts{
	publishDir "Assembly/", mode:'copy'

	input:
	set key,bam,gtf from hisat_bams2.combine(gtf2)

	output:
	file("${key}.gtf") into stringtie_transcripts

	script:
	"""
	stringtie \
	${bam} \
	-G ${gtf} \
	-l ${key} \
	-o ${key}.gtf \
	-p ${task.cpus} \
	"""
}

/* 
 * STEP 6: Merge all of the transcripts
 * INPUT:
 * GTF from STEP 0, GTF files from STEP 5. 
 * OUTPUT:
 * 2 files. gtf_list.txt is not intuitively an 'output',
 * but is being created and used in the process
 * and must be defined. 
 */

process Merge_Transcripts{
	publishDir "Assembly/", mode:'copy'

	input:
	file(gtf) from gtf3
	file(stringtie_gtf) from stringtie_transcripts
	
	output:
	file "merged_transcriptome.gtf" into merged_ch
	file "gtf_list.txt"

	script:
	"""
	echo -e \"${stringtie_gtf.join("\\n")}\" >> gtf_list.txt ;
	stringtie --merge \
	-p 4 \
	-G ${gtf} \
	-o merged_transcriptome.gtf \
	gtf_list.txt
	"""
}

// We will need the merged transcript file twice, place in 2 channels. 

merged_ch.into {merged1; merged2}

/*
 * STEP 7: GFFcompare to check quality of 
 *         merged transcriptome vs. reference. 
 * INPUT: 
 * GTF file from STEP 0, merged transcriptome STEP 6.
 * OUTPUT:
 * 4 files with the chosen prefix.
 */ 

process Assess_Model_Transcriptome{
        publishDir "gffcompare/", mode:'copy'

        input:
        file(gtf) from gtf4
        file(merged_trans) from merged1

        output:
        file "*stringtie_merged.*" into output_ch

        script:
        """
        gffcompare \
        -r ${gtf} \
        -o stringtie_merged \
        ${merged_trans}
	"""
}

/*
 * STEP 8: Estimate transcript abundance in bams 
 *         vs. merged transcriptome
 * INPUT:
 * bam and merged GTF, combine like in STEP 5.
 * OUTPUT:
 * A directory for each sample
 */

process Estimate_Abundance{
        publishDir "Ballgown/", mode:'copy'

        input:
        set key,bam,merged from hisat_bams3.combine(merged2)

        output:
        file("${key}") into ballgown_data

        script:
        """
        stringtie \
        -p ${task.cpus} \
        -G ${merged} \
        -o ${key}/${key}.gtf \
        ${bam} \
        -e
        """
}
