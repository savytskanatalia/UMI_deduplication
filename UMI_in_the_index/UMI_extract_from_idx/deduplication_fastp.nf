//  reads three files with similar IDs - R1,R2 and index >>> to extract UMI from index and put it in read IDs
// "/path/*{index,R1,R2}_001.fastq.gz" 
// from path or from pairs? channel.fromFilePairs
// as fastp tries to work in pairs to utilize it I process R1+index and R2+index as separate instances
// need to use a single core only
// implement pecheck?! just in case to make sure the pairs are in synch
// make it two runs - first will extract UMIs, second take care of adapters in proper pairs with extracted UMIs; we do NOT apply control filtering
// --disable_adapter_trimming --disable_quality_filtering 
// --trim_poly_g because we had nova_seq 

params.my_files = "/path/*{index,R1,R2}_001.fastq.gz" 



process FASTP_UMI {
    publishDir 'results'
    input:
    tuple val(sampleId), file(reads)
    output:
    path("*_R{1,2}_PROC.fastq.gz")

    script:
    """
    echo HI THERE ${reads[0]} ${reads[1]} ${reads[2]} $sampleId

    fastp -i ${reads[0]} -I ${reads[2]} -o ${sampleId}_R1_UMI.fastq.gz -O ${sampleId}_trash1.out.fq --umi --umi_loc=read2 --umi_len=12 -Q -A -L -w 1 -p -j ${sampleId}_R1_fastp.json -h  ${sampleId}_R1_fastp.html --disable_adapter_trimming --disable_quality_filtering 
    fastp -i ${reads[1]} -I ${reads[2]} -o ${sampleId}_R2_UMI.fastq.gz -O ${sampleId}_trash2.out.fq --umi --umi_loc=read2 --umi_len=12 -Q -A -L -w 1 -p -j ${sampleId}_R2_fastp.json -h  ${sampleId}_R2_fastp.html --disable_adapter_trimming --disable_quality_filtering 
    fastp -i ${sampleId}_R1_UMI.fastq.gz -I ${sampleId}_R2_UMI.fastq.gz -o ${sampleId}_R1_PROC.fastq.gz -O ${sampleId}_R2_PROC.fastq.gz --detect_adapter_for_pe -w 1 -p -j ${sampleId}_fastp.json -h  ${sampleId}_fastp.html  --failed_out ${sampleId}_failed.fastq.gz --trim_poly_g --disable_quality_filtering 

    """
}





workflow {
  fastq_files = Channel.fromFilePairs(params.my_files, size: -1)
  FASTP_UMI(fastq_files)
}

