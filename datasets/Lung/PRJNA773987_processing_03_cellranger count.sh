
# It is only used for cellranger count!
#cd ~
#sudo vim .bashrc
#export PATH=/opt/cellranger-6.1.2:$PATH

#echo "Configure Cellranger environment"
#cd ~/cellranger-6.1.2
#source sourceme.bash

cd /data/mengxu/data/PRJNA773987

for i in 29 30 31 32 33 34 35 36 37 38 39 ;do
  {
  
  cd /data/mengxu/data/PRJNA773987/SRR166680${i}
  
  echo "SRR166680${i} cellranger count start!"
  cellranger count --id=SRR166680${i} \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/PRJNA773987/SRR166680${i} \--sample=SRR166680${i} \--localcores=16 \--localmem=64 #\--nosecondary
  # nosecondary 只获得表达矩阵，不进行后续的降维、聚类和可视化分析
  echo "SRR166680${i}_f cellranger count start!"
  cellranger count --id=SRR166680${i}_f \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/PRJNA773987/SRR166680${i} \--sample=SRR166680${i}_f \--localcores=16 \--localmem=64 #\--nosecondary
  
  #Del bam file
  rm  -f /data/mengxu/data/PRJNA773987/SRR166680${i}/SRR166680${i}/outs/possorted_genome_bam.bam
  rm  -f /data/mengxu/data/PRJNA773987/SRR166680${i}/SRR166680${i}_f/outs/possorted_genome_bam.bam
  
  cd /data/mengxu/data/PRJNA773987
  
  } 
done

