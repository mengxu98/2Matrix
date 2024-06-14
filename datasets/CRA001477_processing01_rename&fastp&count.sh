

#download from: http://download.cncb.ac.cn/gsa/CRA001477/
echo "Configure Cellranger environment"
#cd ~/cellranger-6.1.2
#source sourceme.bash
#cd /data/mengxu/data/PRJNA773987/CRR0492_results
#cd /data/mengxu/data/PRJNA773987/SRR166680_f_results
#cellranger count --id=SRR16668029_f \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/PRJNA773987/SRR16668029 \--sample=SRR16668029_f
#
#cellranger count --id=CRR049231 \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001477/CRR049231 \--sample=CRR049231

cd /data/mengxu/data/CRA001477

for i in 27 28 29 30 ;do
  {
  
  if [ ! -d CRR0492${i} ]
  then
  	mkdir CRR0492${i}
  fi
  
  echo "CRR0492${i}_S1_L001_R1_001.fastq.gz is not available, start move file!"
  mv CRR0492${i}_f1.fastq.gz CRR0492${i}_S1_L001_R1_001.fastq.gz
  mv CRR0492${i}_S1_L001_R1_001.fastq.gz /data/mengxu/data/CRA001477/CRR0492${i}
  
  mv CRR0492${i}_r2.fastq.gz CRR0492${i}_S1_L001_R2_001.fastq.gz
  mv CRR0492${i}_S1_L001_R2_001.fastq.gz /data/mengxu/data/CRA001477/CRR0492${i}
  
  echo "CRR0492${i} files are available!"
  
  cd /data/mengxu/data/CRA001477/CRR0492${i}
  
  echo "CRR0492${i} cellranger count start!"
  cellranger count --id=CRR0492${i} \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001477/CRR0492${i} \--sample=CRR0492${i} \--localcores=16 \--localmem=64
  
  if [ ! -f CRR0492${i}_f_S1_L001_R2_001.fastq.gz ]
    then
  	echo "CRR0492${i}_f fastp QC start!"
  fastp -i CRR0492${i}_S1_L001_R1_001.fastq.gz -I CRR0492${i}_S1_L001_R2_001.fastq.gz -o CRR0492${i}_f_S1_L001_R1_001.fastq.gz -O CRR0492${i}_f_S1_L001_R2_001.fastq.gz -w 16 --html CRR0492${i}_QC.html --json CRR0492${i}_QC_R1.json
  fi

  echo "CRR0492${i}_f cellranger count start!"
  cellranger count --id=CRR0492${i}_f \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001477/CRR0492${i} \--sample=CRR0492${i}_f \--localcores=16 \--localmem=64
  
  
  #cd 
  #这行命令只需要运行一次，在 /data/mengxu/usefulaf 路径下运行
  #singularity exec --cleanenv \--bind $AF_SAMPLE_DIR:/workdir \--pwd /usefulaf/bash usefulaf.sif \./simpleaf index \-f /workdir/human_CR_3.0/fasta/genome.fa \-g /workdir/human_CR_3.0/genes/genes.gtf \-l 91 -t 16 -o /workdir/human_CR_3.0_splici
  
  
  #singularity exec --cleanenv \--bind $AF_SAMPLE_DIR:/workdir \--pwd /usefulaf/bash usefulaf.sif \./simpleaf quant \-1 /workdir/CRR049230_S1_L001_R1_001.fastq.gz \-2 /workdir/CRR049230_S1_L001_R2_001.fastq.gz \-i /workdir/human_CR_3.0_splici/index \-o /workdir/quants/pbmc1k_v3 \-f u -c v3 -r cr-like \-m /workdir/human_CR_3.0_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \-t 16
  
  #Del bam file
  rm  -f /data/mengxu/data/CRA001477/CRR0492${i}/CRR0492${i}/outs/possorted_genome_bam.bam
  rm  -f /data/mengxu/data/CRA001477/CRR0492${i}/CRR0492${i}_f/outs/possorted_genome_bam.bam
  
  cd /data/mengxu/data/CRA001477
  
  } 
done

