
#It is only used to process fastq files!

cd /data/mengxu/data/E-MTAB-6149

for i in 1247 BT1249 BT1A BT1B BT1C BT2A BT2B ;do
  {
  
  if [ ! -d ${i} ]
  then
  	mkdir ${i}
  	
  fi
  mv ${i}.R1.fastq.gz /data/mengxu/data/E-MTAB-6149/${i}
  mv ${i}.R2.fastq.gz /data/mengxu/data/E-MTAB-6149/${i}
  mv ${i}.R3.fastq.gz /data/mengxu/data/E-MTAB-6149/${i}
  	
  cd /data/mengxu/data/E-MTAB-6149/${i}
  
  mv ${i}.R1.fastq.gz ${i}_S1_L001_R1_001.fastq.gz
  mv ${i}.R2.fastq.gz ${i}_S1_L001_I1_001.fastq.gz
  mv ${i}.R3.fastq.gz ${i}_S1_L001_R2_001.fastq.gz
  
  cellranger count --id=${i} \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6149/${i} \--sample=${i} \--localcores=16 \--localmem=64 #\--nosecondary
  
  if [ ! -f ${i}_f_S1_L001_R2_001.fastq.gz ]
    then
  	echo "${i} fastp QC start!"
  fastp -i ${i}_S1_L001_R1_001.fastq.gz -I ${i}_S1_L001_R2_001.fastq.gz -o ${i}_f_S1_L001_R1_001.fastq.gz -O ${i}_f_S1_L001_R2_001.fastq.gz -w 16 --html ${i}_QC.html --json ${i}_QC_R1.json
  fi

  echo "${i}_f cellranger count start!"
  cellranger count --id=${i}_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6149/${i} \--sample=${i}_f \--localcores=16 \--localmem=64
  
  
  #cellranger count --id=BT1249 \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6149/BT1249 \--sample=BT1249 \--localcores=16 \--localmem=64 #\--nosecondary
  
  #Del bam file
  rm  -f /data/mengxu/data/E-MTAB-6149/${i}/${i}/outs/possorted_genome_bam.bam
  rm  -f /data/mengxu/data/E-MTAB-6149/${i}/${i}_f/outs/possorted_genome_bam.bam
  
  cd /data/mengxu/data/E-MTAB-6149
  
  } 
done

