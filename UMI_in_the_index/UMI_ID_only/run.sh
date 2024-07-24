

# Call to add UMI IDs (UMIs have been already removed from the reads as well as indices by the qubic (sequencing facility/company/collaborator...); and UMIs were save as sequences into separate index files)

nextflow UMI_ID.nf -c deduplication_fastp.config --my_files "/mnt/md0/natalia/EPO/raw_data/polyA_RNA/*{index,R1,R2}_001.fastq.gz" -with-docker savytskanatalia/quantefication2 --outdir "umi_extracted"

# merge lane fastqs

# Mapping command (example for sequencing data, that was chemically modified and may have become sensitive to qiality filtering):
nextflow run /mnt/md0/natalia/scripts/star_pe_dsl1.nf -c /mnt/md0/natalia/scripts/star_pe_dsl1.config --reads "polyA_RNA/*_R{1,2}.fastq.gz" --sorted 2 --read_length 105 --staropt "--alignEndsType EndToEnd --quantTranscriptomeBan IndelSoftclipSingleend" -with-docker savytskanatalia/quantefication2 --outdir "results" --star_index "/mnt/md0/natalia/mm39_105bp/" --threads 5

# deduplication

