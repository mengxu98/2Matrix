
#It is only used to process *.sra files to fastq files!

cd /data/mengxu/data/PRJNA773987

for i in 29 30 31 32 33 34 35 36 37 38 39 ;do
  {
  
  if [ -f SRR166680${i} ]
  then
  	mv SRR166680${i} SRR166680${i}.sra
  fi
  
  if [ ! -d SRR166680${i} ]
  then
  	mkdir SRR166680${i}
  	mv SRR166680${i}.sra /data/mengxu/data/PRJNA773987/SRR166680${i}
  fi
  
  cd /data/mengxu/data/PRJNA773987/SRR166680${i}
  
  fastq-dump --split-files --gzip SRR166680${i}.sra
  
  mv SRR166680${i}_1.fastq.gz SRR166680${i}_S1_L001_I1_001.fastq.gz
  mv SRR166680${i}_2.fastq.gz SRR166680${i}_S1_L001_R1_001.fastq.gz
  mv SRR166680${i}_3.fastq.gz SRR166680${i}_S1_L001_R2_001.fastq.gz
  
  cd /data/mengxu/data/PRJNA773987
  
  } 
done

