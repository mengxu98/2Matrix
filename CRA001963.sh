# Example: CRA001963
# Make sure the Cellranger and fastp softwares are installed on your device
# Cellranger and References download: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
echo "Configure Cellranger environment......"
# cd ~/cellranger-6.1.2
# source sourceme.bash

cd ~ # Or set by yourself
mkdir CRA001963
cd CRA001963

for i in 22 23 24 25 26 27 28 29 30 31; do
  {
    # Download fastq files from: http://download.cncb.ac.cn/gsa/CRA001963/
    wget http://download.cncb.ac.cn/gsa/CRA001963/CRR0730${i}/CRR0730${i}_f1.fastq.gz
    wget http://download.cncb.ac.cn/gsa/CRA001963/CRR0730${i}/CRR0730${i}_f2.fastq.gz

    if [ ! -d CRR0730${i} ]; then
      mkdir CRR0730${i}
    fi

    echo "CRR0730${i}_S1_L001_R1_001.fastq.gz is not available, start move file......"
    mv CRR0730${i}_f1.fastq.gz CRR0730${i}_S1_L001_R1_001.fastq.gz
    mv CRR0730${i}_S1_L001_R1_001.fastq.gz CRR0730${i}

    mv CRR0730${i}_r2.fastq.gz CRR0730${i}_S1_L001_R2_001.fastq.gz
    mv CRR0730${i}_S1_L001_R2_001.fastq.gz CRR0730${i}

    echo "CRR0730${i} files are available......"

    cd CRR0730${i}

    echo "Running cellranger count for CRR0730${i}......"
    # Please note that the parameter value of 'transcriptome'(References) must be specified by yourself
    cellranger count --id=CRR0730${i} \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001963/CRR0730${i} \--sample=CRR0730${i} \--localcores=16 \--localmem=64

    # QC by fastp
    if [ ! -f CRR0730${i}_f_S1_L001_R2_001.fastq.gz ]; then
      echo "Running fastp QC for CRR0730${i}......"
      fastp -i CRR0730${i}_S1_L001_R1_001.fastq.gz -I CRR0730${i}_S1_L001_R2_001.fastq.gz -o CRR0730${i}_fastp_S1_L001_R1_001.fastq.gz -O CRR0730${i}_fastp_S1_L001_R2_001.fastq.gz -w 16 --html CRR0730${i}_QC.html --json CRR0730${i}_QC_R1.json
    fi

    echo "Running cellranger count for CRR0730${i}_fastp......"
    # Please note that the parameter value of 'transcriptome' must be specified by yourself
    cellranger count --id=CRR0730${i}_fastp \--transcriptome=/data/mengxu/Human_reference/refdata-gex-GRCh38-2020-A \--fastqs=/data/mengxu/data/CRA001963/CRR0730${i} \--sample=CRR0730${i}_fastp \--localcores=16 \--localmem=64

    # Del bam file
    rm -f CRR0730${i}/CRR0730${i}/outs/possorted_genome_bam.bam
    rm -f CRR0730${i}/CRR0730${i}_fastp/outs/possorted_genome_bam.bam

    cd -
  }
done
