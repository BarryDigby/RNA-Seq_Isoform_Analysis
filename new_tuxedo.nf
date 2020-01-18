#!/usr/bin/env nextflow

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

process Extract_Splice_Sites {
	publishDir "Reference/", mode:'copy'
	
	input:
	file(gtf) from gtf1

	output: 
	file "chr22.ss" into splice_sites, splice_sites_
	file "chr22.exon" into exon_sites

	script:
	"""
	hisat2_extract_splice_sites.py ${gtf} > chr22.ss
	hisat2_extract_exons.py ${gtf} > chr22.exon
	"""
}

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

process Hisat_Alignment {
	publishDir "BAMS/", mode:'copy'

	input:
	set val(key), file(reads) from read_ch
	file hisat_index from hs2_indices.collect()
	file ss from splice_sites_.collect()
	
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

hisat_bams1.into { hisat_bams2; hisat_bams3 }

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

process Merge_Transcripts{
	publishDir "Assembly/", mode:'copy'

	input:
	file(gtf) from gtf3
	file(gtf_to_merge) from stringtie_transcripts.collect()
	
	output:
	file "merged_transcriptome.gtf" into merged_ch
	file "gtf_list.txt"

	script:
	"""
	echo -e \"${gtf_to_merge.join("\\n")}\" >> gtf_list.txt ;
	stringtie --merge \
	-p 4 \
	-G ${gtf} \
	-o merged_transcriptome.gtf \
	gtf_list.txt
	"""
}

merged_ch.into {merged1; merged2}

process Assess_Model{
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
