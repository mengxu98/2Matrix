
# It is only used for data quality control!

cd /data/mengxu/data/PRJNA773987

for i in 29 30 31 32 33 34 35 36 37 38 39 ;do
  {
  
  cd /data/mengxu/data/PRJNA773987/SRR166680${i}
  
  if [ ! -f SRR166680${i}_f_S1_L001_R2_001.fastq.gz ]
  then
      echo "SRR166680${i}_f fastp QC start!"
      
      fastp -i SRR166680${i}_S1_L001_R1_001.fastq.gz -I SRR166680${i}_S1_L001_R2_001.fastq.gz -o SRR166680${i}_f_S1_L001_R1_001.fastq.gz -O SRR166680${i}_f_S1_L001_R2_001.fastq.gz -w 16 --html SRR166680${i}_QC.html --json SRR166680${i}_QC.json

  fi
  
  echo SRR166680${i} QC had completed!
  
  cd /data/mengxu/data/PRJNA773987
  
  } 
done

