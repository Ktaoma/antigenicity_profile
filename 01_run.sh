########################################################################
## REFERENCE-BASED GENOME ASSEMBLY
########################################################################

# Make result directories
mkdir -p result/fastq
mkdir -p result/stat
mkdir -p result/tmp
mkdir -p result/bam
mkdir -p result/consensus

# Summary statistics of raw fastq (data/03_tmp/b49.fastq.gz)
seqkit stat $1 > result/stat/summary_statistcs.tsv

#Select file based on pattern
seqkit locate -p AGTAGAAACAAGG $1 | gzip -9 > result/tmp/reverse.tsv.gz
seqkit locate -p AGCAAAAGCAGG  $1 | gzip -9 > result/tmp/forward.tsv.gz

# Select the read id with correct orientation
Rscript ./script/02_parse_primer_pattern.R result/tmp both
seqkit grep -f result/tmp/QC_id.tsv $1 | gzip -9 > result/fastq/filtered.fq.gz
seqkit stat result/fastq/filtered.fq.gz | tail -n+2 >> result/stat/summary_statistcs.tsv
rm -rf result/fastq/*.tsv.gz

# Simply Trim primer by removing first and last 50 bases at the sequence 
seqkit seq -m 500 -Q 10 result/fastq/filtered.fq.gz |\
  seqkit subseq -r 50:-50 | gzip -9 > result/fastq/filtered_02.fq.gz
seqkit stat result/fastq/filtered_02.fq.gz | tail -n+2 >> result/stat/summary_statistcs.tsv
  
# Run makeblastdb just for single time to make database for blastn 
makeblastdb -in db/sequences_DNA.fasta -dbtype nucl

# This will give you reads mapped on influenza virus sequence
seqkit fq2fa result/fastq/filtered_02.fq.gz |\
  blastn -query - -db db/sequences_DNA.fasta -outfmt '6 qseqid sseqid pident length qlen' > result/tmp/blast_hits.tsv

# Get only read with 80 percent similarity mapping
cat result/tmp/blast_hits.tsv | awk '{prop=100*($4/$5);print $0 "\t" prop}'| awk '$3 >= 80 && $6 >= 80' > result/tmp/blast_hits_filtered.tsv

# 4th index is HA segment and 6th index is NA segment
# we only use HA segment for calculating antigencity score
# for the full pipeline we might need to report the whole genome of Flu (8 segment)
# This for testing only for antigenicity prediction
for k in 4 6;
do

   # Get read of k segment
   cat db/sequences_DNA.fasta | grep ">" | grep "segment $k"  |\
        awk '{print $1}' | sed 's/>//g' | grep -f - result/tmp/blast_hits_filtered.tsv | cut -f1 | sort -u |\
            seqkit grep -f - result/fastq/filtered_02.fq.gz | gzip -9 > result/fastq/filtered_${k}.fq.gz
            
  
   # Get the name of best reference of k segment
   name=$(cat db/sequences_DNA.fasta | grep ">" |\
            grep "segment $k" | awk '{print $1}' |\
              sed 's/>//g' | grep -f - result/tmp/blast_hits_filtered.tsv |\
                cut -f2 | sort | uniq -c | sort -k1,1nr | awk 'NR==1{print $2}')

    # Generating consensus
    seqkit fx2tab db/sequences_DNA.fasta | grep $name | seqkit tab2fx > result/tmp/best_ref.fa
    seqkit fx2tab db/sequences_DNA.fasta | grep $name | seqkit tab2fx |\
        minimap2 -a - result/fastq/filtered_${k}.fq.gz |\
            samtools sort -o result/tmp/filtered_${k}.bam && samtools index result/tmp/filtered_${k}.bam

    samtools sort result/tmp/filtered_${k}.bam -o result/bam/filtered_${k}_sorted.bam
    samtools index result/bam/filtered_${k}_sorted.bam
    
    samtools consensus -m simple -d 1 -A result/bam/filtered_${k}_sorted.bam |\
           seqkit fx2tab | awk -v sample="sample" -v seg=$k 'BEGIN{OFS="\t"}{print sample "_" seg,$2}' |\
                    seqkit tab2fx > result/consensus/draft_segment_${k}.fasta
    
    # Calculate depth
    samtools depth result/bam/filtered_${k}_sorted.bam |\
           awk -v sample="sample" -v seg=$k '{print sample "_" seg,$2,$3}' > result/stat/depth_segment_${k}.tsv
  
    # Summary statistics of genome assembly
    seqkit stat result/fastq/filtered_${k}.fq.gz | tail -n+2 >> result/stat/summary_statistcs.tsv
  
done
