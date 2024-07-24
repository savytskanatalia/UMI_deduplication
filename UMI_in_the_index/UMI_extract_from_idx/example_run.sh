# example scripts to use if UMIs sequences from *index.fastq.gz were NOT removed from the reads themselves (*{1,2}.fastq.gz)
# and one needs to extract those from reads and place the IDX into read names
# case if UMI is present in both reads from a pair
# this strategy relies on fastp for extraction of UMI and adding its sequence to read names

# 3 raw sequencing files present: Read1, Read2, Index
nextflow run /mnt/md0/natalia/EPO/raw_data/trial.nf -c /mnt/md0/natalia/EPO/raw_data/deduplication_fastp.config -with-docker savytskanatalia/quantefication2 --my_files "raw_data/polyA_RNA/*_{R1,R2,index}.fastq.gz"


# check to make sure pairs became unpaired (and they did, as fastp processing leads to loss of initial read order)

nextflow run pecheck.nf -c /mnt/md0/natalia/EPO/raw_data/deduplication_fastp.config --my_files "raw_data/processed/*_R{1,2}_UMI.fastq.gz" -with-docker savytskanatalia/quantefication2
for f in results/*.json; do echo "" >> synch.stat  && echo "$f" | tr "\n" "\t" >> synch.stat && grep -Po '"result":.*?[^\\]"|"message":.*?[^\\]"' $f | tr "\n" "\t" >> synch.stat ; done
# they did


# re-pairing using bbtools 
while read i; do echo $i; repair.sh in1=${i}_R1_PROC.fastq.gz in2=${i}_R2_PROC.fastq.gz out1=results/${i}_R1_RP.fastq.gz out2=results/${i}_R2_RP.fastq.gz outs=results/${i}_singletons.fq repair; done < sample_id.txt

# now can QC as pairs and remove any leftover adapters or polyG artifacts
# per sample example (X)
fastp -i X_R1_UMI.fastq.gz -I X_R2_UMI.fastq.gz -o X_R1_PROC.fastq.gz -O X_R2_PROC.fastq.gz --detect_adapter_for_pe -w 1 -p -j X_fastp.json -h X_fastp.html --failed_out X_failed.fastq.gz --trim_poly_g 

# map in any convenient way

# Now deduplicate - preserve a single read per molecular identifier and drop actual duplicates
nextflow run umi_dedup_dsl1.nf -c umi_dedup_dsl1.config -with-docker savytskanatalia/umitools --bam_files "dir/*.bam" --outdir "dir/results1" --paired 2 --umi_opt "--umi-separator=':'"
