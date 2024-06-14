

#
echo "Configure Cellranger environment"
#cd ~/cellranger-6.1.2
#source sourceme.bash

cd /data/mengxu/data/E-MTAB-6653/BT1375

cellranger count --id=BT1375 \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/BT1375 \--sample=BT1375 \--localcores=16 \--localmem=64

cellranger count --id=BT1375_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/BT1375 \--sample=BT1375_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/BT1375/BT1375/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/BT1375/BT1375_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/BT1376

cellranger count --id=BT1376 \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/BT1376 \--sample=BT1376 \--localcores=16 \--localmem=64

cellranger count --id=BT1376_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/BT1376 \--sample=BT1376_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/BT1376/BT1376/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/BT1376/BT1376_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/BT1377

cellranger count --id=BT1377 \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/BT1377 \--sample=BT1377 \--localcores=16 \--localmem=64

cellranger count --id=BT1377_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/BT1377 \--sample=BT1377_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/BT1377/BT1377/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/BT1377/BT1377_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/BT1378

cellranger count --id=BT1378 \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/BT1378 \--sample=BT1378 \--localcores=16 \--localmem=64

cellranger count --id=BT1378_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/BT1378 \--sample=BT1378_f \--localcores=16 \--localmem=64
#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/BT1378/BT1378/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/BT1378/BT1378_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/scrBT1425_hg19

cellranger count --id=scrBT1425_hg19 \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1425_hg19 \--sample=scrBT1425_hg19 \--localcores=16 \--localmem=64

cellranger count --id=scrBT1425_hg19_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1425_hg19 \--sample=scrBT1425_hg19_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1425_hg19/scrBT1425_hg19/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1425_hg19/scrBT1425_hg19_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/scrBT1426_hg19

cellranger count --id=scrBT1426_hg19 \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1426_hg19 \--sample=scrBT1426_hg19 \--localcores=16 \--localmem=64

cellranger count --id=scrBT1426_hg19_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1426_hg19 \--sample=scrBT1426_hg19_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1426_hg19/scrBT1426_hg19/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1426_hg19/scrBT1426_hg19_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/scrBT1427_hg19

cellranger count --id=scrBT1427_hg19 \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1427_hg19 \--sample=scrBT1427_hg19 \--localcores=16 \--localmem=64

cellranger count --id=scrBT1427_hg19_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1427_hg19 \--sample=scrBT1427_hg19_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1427_hg19/scrBT1427_hg19/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1427_hg19/scrBT1427_hg19_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/scrBT1428m

cellranger count --id=scrBT1428m \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1428m \--sample=scrBT1428m \--localcores=16 \--localmem=64

cellranger count --id=scrBT1428m_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1428m \--sample=scrBT1428m_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1428_hg19/scrBT1428_hg19/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1428_hg19/scrBT1428_hg19_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/scrBT1429m

cellranger count --id=scrBT1429m \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1429m \--sample=scrBT1429m \--localcores=16 \--localmem=64

cellranger count --id=scrBT1429m_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1429m \--sample=scrBT1429m_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1429_hg19/scrBT1429_hg19/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1429_hg19/scrBT1429_hg19_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/scrBT1430m

cellranger count --id=scrBT1430m \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1430m \--sample=scrBT1430m \--localcores=16 \--localmem=64

cellranger count --id=scrBT1430m_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1430m \--sample=scrBT1430m_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1430m/scrBT1430m/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1430m/scrBT1430m_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/scrBT1431m

cellranger count --id=scrBT1431m \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1431m \--sample=scrBT1431m \--localcores=16 \--localmem=64

cellranger count --id=scrBT1431m_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1431m \--sample=scrBT1431m_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1431m/scrBT1431m/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1431m/scrBT1431m_f/outs/possorted_genome_bam.bam

#------------------------------------------------------------------------------
cd /data/mengxu/data/E-MTAB-6653/scrBT1432m

cellranger count --id=scrBT1432m \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1432m \--sample=scrBT1432m \--localcores=16 \--localmem=64

cellranger count --id=scrBT1432m_f \--transcriptome=/data/mengxu/Human_reference/hg19 \--fastqs=/data/mengxu/data/E-MTAB-6653/scrBT1432m \--sample=scrBT1432m_f \--localcores=16 \--localmem=64

#Del bam file
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1432m/scrBT1432m/outs/possorted_genome_bam.bam
rm  -f /data/mengxu/data/E-MTAB-6653/scrBT1432m/scrBT1432m_f/outs/possorted_genome_bam.bam

