//  reads three files with similar IDs - R1,R2 and index >>> to extract UMI from index and put it in read IDs
//   pecheck -i ${sampleId}_R1_PROC.fastq.gz -I ${sampleId}_R2_PROC.fastq.gz -j ${sampleId}.json


params.my_files = "/path/*{R1,R2}.fastq.gz" 



process PECHECK {
    publishDir 'results'
    input:
    tuple val(sampleId), file(reads)
    output:
    path("*_pechecked.json")

    script:
    """
    echo HI THERE ${reads[0]} ${reads[1]} $sampleId
    pecheck -i ${reads[0]} -I ${reads[1]} -j ${sampleId}_pechecked.json

    """
}





workflow {
  fastq_files = Channel.fromFilePairs(params.my_files)
  PECHECK(fastq_files)
}

