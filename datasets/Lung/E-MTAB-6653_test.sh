#!/bin/bash

path=/data/mengxu/data/E-MTAB-6653
function getdir(){
  for element in `ls $fd`
  do  
    dir_or_file=$fd"/"$element
    if [ -d $path$dir_or_file ]
    then 
      #getdir $dir_or_file
      cd $path$dir_or_file
      echo $path$dir_or_file
      
      #path=$1
      files=$(ls $path$dir_or_file)
      for filename in $files
      do
        echo $filename fastp start!
        echo $files
        #fastp -i BT1375_S5_L004_R1_001.fastq.gz -I BT1375_S5_L004_R2_001.fastq.gz -o BT1375_f_S5_L004_R1_001.fastq.gz -O BT1375_f_S5_L004_R2_001.fastq.gz -w 16 --html BT1375_S5_L004_QC.html --json BT1375_S5_L004_QC.json
        
        echo $filename fastp done!
      done
        
      cd $path
    else
      
      #echo $dir_or_file
      echo This a file!
        
    fi  
  done
}
root_dir="/opt/datas"
getdir $root_dir
