
## mapping using hisat2
```bash
## example hisat2 command
hisat2 -p 4 --summary-file hisat.summary  -x /home/lakhanp/database/Human/GRCh38p12.gencode30/hisat2_index/GRCh38.primary_assembly.genome.chr.fa -1 /home/lakhanp/analysis/39_Ben_RNAseq2/raw_data/CL_RNA_seq_mix05_batch1/25_DMSO_12hr_1_S54_L006_R1_001.fastq.gz,/home/lakhanp/analysis/39_Ben_RNAseq2/raw_data/CL_RNA_seq_mix05_batch23/25_DMSO_12hr_1_S1_L001_R1_001.fastq.gz,/home/lakhanp/analysis/39_Ben_RNAseq2/raw_data/CL_RNA_seq_mix05_batch23/25_DMSO_12hr_1_S1_L002_R1_001.fastq.gz -2 /home/lakhanp/analysis/39_Ben_RNAseq2/raw_data/CL_RNA_seq_mix05_batch1/25_DMSO_12hr_1_S54_L006_R2_001.fastq.gz,/home/lakhanp/analysis/39_Ben_RNAseq2/raw_data/CL_RNA_seq_mix05_batch23/25_DMSO_12hr_1_S1_L001_R2_001.fastq.gz,/home/lakhanp/analysis/39_Ben_RNAseq2/raw_data/CL_RNA_seq_mix05_batch23/25_DMSO_12hr_1_S1_L002_R2_001.fastq.gz | samtools view -bS - | samtools sort  -O bam -o DMSO_12hr_1_hisat2.bam

```


## index BAM files
```bash
for i in `cat sample_name.list`
do
cd $i
printf "samtools index %s_hisat2.bam
samtools flagstat %s_hisat2.bam > alignment.stats\n\n" $i $i >> generalJob.sh
cd ..
done
```


## run stringTie
```bash
for i in `cat sample_name.list`
do
cd $i
printf "##Run stringtie: just counting the transcripts and no assembly
stringtie %s_hisat2.bam -p 8 -e -B -G /home/lakhanp/database/Human/GRCh38p12.gencode30/annotation/GRCh38p12.gencode.v30.basic.annotation.sorted.gtf -o stringTie_%s/%s.gtf
error_exit \$?\n\n" $i $i $i >> generalJob.sh
cd ..
done
```


## rRNA counts
```bash
for i in `cat sample_name.list`
do
cd $i
printf "$i\t"
samtools view -@ 4 -c -L /home/lakhanp/database/Human/GRCh38p12.gencode30/annotation/GRCh38p12.rRNA2.bed $```i```_hisat2.bam
cd ..
done > rRNA_counts.tab
```

















