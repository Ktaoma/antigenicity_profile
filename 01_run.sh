#makeblastdb -in ../../../db/influenza/phase2/sequences_DNA.fasta -dbtype nucl


i=$1
type=$2
leng=$3

#sn=$(echo $i | cut -d '/' -f7 | sed -e 's/.fastq.gz//g')
mkdir result
seqkit stat $i > result/sample.tsv



function h {
  path=result/$sn/${type}
mkdir -p $path

echo 'Get summary statistic of uncurated file'
seqkit stat $i -Ta | sed -e 's/,//g' >> ${path}/read_count.tsv

echo 'Get primer pattern'
pp=$path/pattern/$type
mkdir -p $pp
seqkit locate -m 2 -p AGTAGAAACAAGG $i | awk -v ss=$sn '{print $0,ss}' | gzip -9 > ${pp%/type}/reverse.tsv.gz
seqkit locate -m 2 -p AGCAAAAGCAGG  $i | awk -v ss=$sn '{print $0,ss}' | gzip -9 > ${pp%/type}/forward.tsv.gz

echo 'Extract primer'
mkdir -p ${pp}
Rscript ./02_parse_primer_pattern.R $pp $type

echo 'Export as Fastq'
seqkit grep -f $pp/QC_id.tsv $i | gzip -9 > $path/filtered.fq.gz

echo 'Trim primer'
cutadapt --revcomp -a AGCAAAAGCAGG...CCTTGTTTCTACT -e 0.25 -m $leng -q 10 -o $path/sample_final_qc.fq.gz $path/filtered.fq.gz
seqkit fq2fa $path/sample_final_qc.fq.gz | gzip -9 > $path/sample_final.fa.gz

echo 'Sort reads into bucket'
zcat $path/sample_final.fa.gz | blastn -query - -db ../../../db/influenza/phase2/sequences_DNA.fasta -outfmt '6 qseqid sseqid pident length qlen' -word_size 28 > $path/blast_hits.tsv
cat ${path}/blast_hits.tsv | awk '{prop=100*($4/$5);print $0 "\t" prop}'| awk '$3 >= 80 && $6 >= 80' > ${path}/blast_hits_filtered.tsv

echo 'Get summary statistics of QC fasq'
seqkit stat $path/sample_final.fa.gz -Ta | sed -e 's/,//g' | awk 'NR!=1 {print $0}' >> ${path}/read_count.tsv

for k in 4;
do

    echo 'Get read of k segment'
    cat ../../../db/influenza/phase2/sequences_DNA.fasta | grep ">" | grep "segment $k"  |\
        awk '{print $1}' | sed 's/>//g' | grep -f - $path/blast_hits_filtered.tsv | cut -f1 | sort -u |\
            seqkit grep -f - $path/sample_final.fa.gz | gzip -9 > $path/sample_final_segment_${k}.fa.gz

    echo 'Get summary statistics of k segment'
    seqkit stat $path/sample_final_segment_${k}.fa.gz -Ta | sed -e 's/,//g' | awk 'NR!=1 {print $0}' >> ${path}/read_count.tsv


    echo 'Get the name of best reference of k segment'
    name=$(cat ../../../db/influenza/phase2/sequences_DNA.fasta | grep ">" | grep "segment $k"  | awk '{print $1}' |\
           sed 's/>//g' | grep -f - ${path}/blast_hits_filtered.tsv | cut -f2 | sort | uniq -c | sort -k1,1nr | awk 'NR==1{print $2}')


    echo 'Remove supplmentary reads'
    seqkit fx2tab ../../../db/influenza/phase2/sequences_DNA.fasta | grep $name | seqkit tab2fx |\
        minimap2 -a - $path/sample_final_segment_${k}.fa.gz |\
            samtools sort -o $path/sample_final_segment_${k}_tmp.bam && samtools index $path/sample_final_segment_${k}_tmp.bam

    samtools view -f 2048  $path/sample_final_segment_${k}_tmp.bam | cut -f1 | sort -u > $path/supplement_${k}.txt
    samtools view -h $path/sample_final_segment_${k}_tmp.bam | grep -vf $path/supplement_${k}.txt | samtools view -b -o $path/sample_final_segment_${k}.bam -


    echo 'Remove read with excessive INDEL and softclip'
    samtools view -h $path/sample_final_segment_${k}.bam | awk '$5==60 {print $1,$6}' > ${path}/ReadCIGAR.tsv

    Rscript ./03_parse_CIGAR.R $path
    mv ${path}/ReadCIGAR.tsv ${path}/ReadCIGAR_${k}.tsv
    mv ${path}/CIGAR_Tab.tsv ${path}/CIGAR_Tab_${k}.tsv

    cat ${path}/CIGAR_Tab_${k}.tsv  | awk '$10 != "keep" {print $1}' > $path/sample_final_segment_${k}_ID_filtered.txt
    samtools view -h $path/sample_final_segment_${k}.bam | grep -vf $path/sample_final_segment_${k}_ID_filtered.txt | samtools view -b -o $path/sample_final_segment_${k}_filtered.bam -
    samtools sort $path/sample_final_segment_${k}_filtered.bam -o $path/sample_final_segment_${k}_filtered_sorted.bam
    samtools index $path/sample_final_segment_${k}_filtered_sorted.bam

    #echo 'Reference based assembly'
    samtools consensus -m simple -d 30 -A $path/sample_final_segment_${k}_filtered_sorted.bam |\
           seqkit fx2tab | awk -v sample=$sn -v seg=$k 'BEGIN{OFS="\t"}{print sample "_" seg,$2}' |\
                    seqkit tab2fx > $path/draft_segment_${k}.fasta

    #echo 'Calculate depth'
    samtools depth $path/sample_final_segment_${k}_filtered_sorted.bam |\
           awk -v sample=$sn -v seg=$k '{print sample "_" seg,$2,$3}' > $path/depth_segment_${k}.tsv

    #echo 'Self alignment'
    samtools fastq $path/sample_final_segment_${k}_filtered_sorted.bam | gzip -9 > $path/sample_final_segment_${k}_filtered_sorted.fq.gz
    minimap2 -a $path/draft_segment_${k}.fasta $path/sample_final_segment_${k}_filtered_sorted.fq.gz > $path/self_aln_${k}.sam
    samtools sort $path/self_aln_${k}.sam -o $path/self_aln_${k}.bam
done  
}
