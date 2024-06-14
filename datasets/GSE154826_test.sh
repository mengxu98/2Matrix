
#!/bin/sh
Folder_A="/data/mengxu/data/GSE154826"   #/GSE154826_amp_batch_ID_1
for file_a in ${Folder_A}/*
do  

    temp_folder=`basename $file_a  .tar.gz`  
    
    mkdir /data/mengxu/data/GSE154826/$temp_folder
    
    if [ ! -d /data/mengxu/data/GSE154826/$temp_folder ]
    then
      tar -zxvf $temp_file -C /data/mengxu/data/GSE154826/$temp_folder
    fi
    
    Folder_B="/data/mengxu/data/GSE154826/$temp_folder"  
    for file_b in ${Folder_B}/*
    do  
      cd /data/mengxu/data/GSE154826/$temp_folder
      temp_file=`basename $file_b "_"`
      
      echo ${temp_file}
      
      #echo ${temp_file#*_*_*-*_*}
      
      echo ${temp_file#*_*_*-*_*_*}
      
      mv ${temp_file#} ${temp_file#*_*_*-*_*_*_*}
      mv ${temp_file#} ${temp_file#*_*_*-*_*_*}
      mv ${temp_file#} ${temp_file#*_*_*-*_*}
      
    done  
done

