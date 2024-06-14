echo "start"
for i in 29 33 34 35;do
{
  
  cd /data/mengxu/data/PRJNA773987/SRR166680${i}_results
  if [ -f "SRR166680${i}_S1_L001_I1_001.fastq.gz" -a -f "SRR166680${i}_S1_L001_R1_001.fastq.gz" -a -f "SRR166680${i}_S1_L001_R2_001.fastq.gz" ]
  then
  	echo "SRR166680${i} files are available, start cellranger count!"
  else
  	echo "SRR166680${i} not find all files"
  	if [ -f "SRR166680${i}_1.fastq.gz" -a -f "SRR166680${i}_2.fastq.gz" -a -f "SRR166680${i}_3.fastq.gz" ]
  		then
  			#mv SRR166680${i}_1.fastq.gz SRR166680${i}_S1_L001_I1_001.fastq.gz
  			#mv SRR166680${i}_2.fastq.gz SRR166680${i}_S1_L001_R1_001.fastq.gz
  			#mv SRR166680${i}_3.fastq.gz SRR166680${i}_S1_L001_R2_001.fastq.gz
  			echo "SRR166680${i} copy files"
  			printf "SRR166680${i} copy"
  		else
  			#fastq-dump --split-files --gzip SRR166680${i}.sra
  			#mv SRR166680${i}_1.fastq.gz SRR166680${i}_S1_L001_I1_001.fastq.gz
  			#mv SRR166680${i}_2.fastq.gz SRR166680${i}_S1_L001_R1_001.fastq.gz
  			#mv SRR166680${i}_3.fastq.gz SRR166680${i}_S1_L001_R2_001.fastq.gz
  			echo "SRR166680${i} start fastq-dump"
  			
  	fi 
  fi 
  
  if [ -f "/data/mengxu/data/PRJNA773987/SRR166680${i}_results/SRR166680${i}.sra" ]
  then
	echo "111"
else
echo "222"
  fi 
} 
done


echo "Configure Cellranger environment"
cd ~/cellranger-6.1.2
source sourceme.bash
#cd /data/mengxu/data/PRJNA773987/SRR16668029_results
#cd /data/mengxu/data/PRJNA773987/SRR16668029_f_results
#cellranger count --id=SRR16668029_f \--transcriptome=/data/mengxu/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/PRJNA773987/SRR16668029 \--sample=SRR16668029_f
#
#cellranger count --id=CRR073031 \--transcriptome=/data/mengxu/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001963/CRR073031 \--sample=CRR073031

#chang shi dan ge fastq wen jian neng bu neng jin xing cellranger count!!!

cd /data/mengxu/data/PRJNA773987

for i in 29 30 31 32 33 34 35 36 37 38 39 ;do
  {
  
  if [ ! -d SRR166680${i} ]
  then
  	mkdir SRR166680${i}
  fi
  
  echo "SRR166680${i}_S1_L001_R1_001.fastq.gz is not available, start move file!"
  mv SRR166680${i}.fastq.gz SRR166680${i}_S1_L001_R1_001.fastq.gz
  mv SRR166680${i}_S1_L001_R1_001.fastq.gz /data/mengxu/data/PRJNA773987/SRR166680${i}
  echo "SRR166680${i}_S1_L001_R1_001.fastq.gz is available!"
  
  cd /data/mengxu/data/PRJNA773987/SRR166680${i}
  
  echo "SRR166680${i} cellranger count start!"
  cellranger count --id=SRR166680${i} \--transcriptome=/data/mengxu/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/PRJNA773987/SRR166680${i} \--sample=SRR166680${i}
  
  if [ ! -f SRR166680${i}_f_S1_L001_R2_001.fastq.gz ]
  then
      echo "SRR166680${i}_f fastp QC start!"
      fastp -i SRR166680${i}_S1_L001_R1_001.fastq.gz -I SRR166680${i}_S1_L001_R2_001.fastq.gz -o SRR166680${i}_f_S1_L001_R1_001.fastq.gz -O SRR166680${i}_f_S1_L001_R2_001.fastq.gz -w 16 --html SRR166680${i}_QC.html --json SRR166680${i}_QC.json
      #fastp -i SRR16668029_S1_L001_I1_001.fastq.gz -o SRR16668029_f_S1_L001_I1_001.fastq.gz -w 16 --html SRR166680${i}_QC.html --json SRR166680${i}_QC.json
      #fastp -i SRR16668029_S1_L001_R1_001.fastq.gz -I SRR16668029_S1_L001_R2_001.fastq.gz -o SRR16668029_f_S1_L001_R1_001.fastq.gz -O SRR16668029_f_S1_L001_R2_001.fastq.gz -w 16 --html SRR16668029_f_QC.html --json SRR16668029_f_QC.json
  fi
  
  echo "SRR166680${i}_f cellranger count start!"
  cellranger count --id=SRR166680${i}_f \--transcriptome=/data/mengxu/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/PRJNA773987/SRR166680${i} \--sample=SRR166680${i}_f
  
  #cellranger count --id=SRR16668029_f \--transcriptome=/data/mengxu/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/PRJNA773987/SRR16668029 \--sample=SRR16668029_f
  
  cd /data/mengxu/data/PRJNA773987
  
  } 
done


cellranger mkfastq --id=SRR16668038 \--run=/data/mengxu/data/PRJNA773987/SRR16668038
