#!/usr/bin/env nextflow
nextflow.enable.dsl=1
/*
================================================================================
                           UMI-TOOLS DEDUPLICATE
================================================================================
"Wrapper" pipeline for umi_tools dedup -I mapped.bam --paired -S deduplicated.bam
--------------------------------------------------------------------------------
 @Repository
 https://github.com/savytskanatalia/
 @Author
 Savytska Natalia, 2024. Genome Biology of Neurodegenerative Diseases, DZNE TU
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run quanTEr.nf -c quanTEr.config --reads "/absolute/path/to/files/*_R{1,2}.fastq.gz" ... -with-docker savytskanatalia/quantefication2 
    Mandatory arguments:
      --bam_files [str]              Full path to mapped files. Example: "input/*.bam" 
      --outdir [str]                 Output directory path. Default: output/      
      --paired [str]		     If sequencing was paired-end (2) or single-end (1)
    """.stripIndent()
}


paired=params.paired
bam=Channel.fromPath(params.bam_files)
outdir=params.outdir




process umi_dedup {
        container='savytskanatalia/umitools'
        echo true
        publishDir "${params.outdir}/dedup_bam", mode: 'move', pattern: "*.bam"
        publishDir "${params.outdir}/raw_bam", mode: 'move', pattern: "*.bai"
        input:
        path(bam)
        val s from umiopt
        output:
        file "*"
        script:
        if( paired == 1)
            """
            echo Processing sample ${bam.baseName}
            samtools index $bam
            umi_tools dedup -I $bam -S ${bam.baseName}_deduplicated.bam


            """
        else if(paired == 2)
            """
            echo Processing sample ${bam.baseName}
            samtools index $bam
            umi_tools dedup -I $bam --paired -S ${bam.baseName}_deduplicated.bam


            """
        else
            error "Invalid protocol indicated: ${paired}. Valid options: 1 - single-end sequencing; 2 - paired-end sequencing. "
}
