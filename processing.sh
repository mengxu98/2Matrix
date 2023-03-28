
#download from: http://download.cncb.ac.cn/gsa/CRA001963/
echo "Configure Cellranger environment"
#cd ~/cellranger-6.1.2
#source sourceme.bash
#cd /data/mengxu/data/PRJNA773987/CRR0730_results
#cd /data/mengxu/data/PRJNA773987/SRR166680_f_results
#cellranger count --id=SRR16668029_f \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/PRJNA773987/SRR16668029 \--sample=SRR16668029_f
#
#cellranger count --id=CRR073031 \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001963/CRR073031 \--sample=CRR073031

cd /data/mengxu/data/CRA001963

for i in 22 23 24 25 26 27 28 29 30 31 ;do
  {
  # Download fastq file
  wget http://download.cncb.ac.cn/gsa/CRA001963/CRR0730${i}/CRR0730${i}_f1.fastq.gz
  wget http://download.cncb.ac.cn/gsa/CRA001963/CRR0730${i}/CRR0730${i}_f2.fastq.gz
  
  if [ ! -d CRR0730${i} ]
  then
  	mkdir CRR0730${i}
  fi
  
  echo "CRR0730${i}_S1_L001_R1_001.fastq.gz is not available, start move file!"
  mv CRR0730${i}_f1.fastq.gz CRR0730${i}_S1_L001_R1_001.fastq.gz
  mv CRR0730${i}_S1_L001_R1_001.fastq.gz /data/mengxu/data/CRA001963/CRR0730${i}
  
  mv CRR0730${i}_r2.fastq.gz CRR0730${i}_S1_L001_R2_001.fastq.gz
  mv CRR0730${i}_S1_L001_R2_001.fastq.gz /data/mengxu/data/CRA001963/CRR0730${i}
  
  echo "CRR0730${i} files are available!"
  
  cd /data/mengxu/data/CRA001963/CRR0730${i}
  
  echo "CRR0730${i} cellranger count start!"
  cellranger count --id=CRR0730${i} \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001963/CRR0730${i} \--sample=CRR0730${i} \--localcores=16 \--localmem=64
  
  if [ ! -f CRR0730${i}_f_S1_L001_R2_001.fastq.gz ]
    then
  	echo "CRR0730${i}_f fastp QC start!"
  fastp -i CRR0730${i}_S1_L001_R1_001.fastq.gz -I CRR0730${i}_S1_L001_R2_001.fastq.gz -o CRR0730${i}_f_S1_L001_R1_001.fastq.gz -O CRR0730${i}_f_S1_L001_R2_001.fastq.gz -w 16 --html CRR0730${i}_QC.html --json CRR0730${i}_QC_R1.json
  fi
  
  echo "CRR0730${i}_f cellranger count start!"
  cellranger count --id=CRR0730${i}_f \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001963/CRR0730${i} \--sample=CRR0730${i}_f \--localcores=16 \--localmem=64
  
  #cd /data/mengxu/usefulaf
  #singularity exec --cleanenv \--bind $AF_SAMPLE_DIR:/workdir \--pwd /usefulaf/bash usefulaf.sif \./simpleaf index \-f /workdir/human_CR_3.0/fasta/genome.fa \-g /workdir/human_CR_3.0/genes/genes.gtf \-l 91 -t 16 -o /workdir/human_CR_3.0_splici
  
  #singularity exec --cleanenv \--bind $AF_SAMPLE_DIR:/workdir \--pwd /usefulaf/bash usefulaf.sif \./simpleaf quant \-1 /workdir/CRR073030_S1_L001_R1_001.fastq.gz \-2 /workdir/CRR073030_S1_L001_R2_001.fastq.gz \-i /workdir/human_CR_3.0_splici/index \-o /workdir/quants/pbmc1k_v3 \-f u -c v3 -r cr-like \-m /workdir/human_CR_3.0_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \-t 16
  
  # Del bam file
  rm  -f /data/mengxu/data/CRA001963/CRR0730${i}/CRR0730${i}/outs/possorted_genome_bam.bam
  rm  -f /data/mengxu/data/CRA001963/CRR0730${i}/CRR0730${i}_f/outs/possorted_genome_bam.bam
  
  cd /data/mengxu/data/CRA001963
  
  } 
done


