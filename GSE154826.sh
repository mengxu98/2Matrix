#!/bin/sh
cd ~
Folder_User="GSE154826" # Please set by yourself
mkdir $Folder_User
cd $Folder_User
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/GSE154826_sample_annots.csv.gz
gunzip GSE154826_sample_annots.csv.gz

id=$(cut -d ',' -f 3 GSE154826_sample_annots.csv)
for i in ${id}; do
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/GSE154826_amp_batch_ID_$i.tar.gz
done

for file_a in ${Folder_User}/*; do
    temp_folder=$(basename $file_a .tar.gz)

    mkdir $temp_folder

    if [ ! -d $temp_folder ]; then
        tar -zxvf $file_a -C $temp_folder
    fi

    Folder_B="$temp_folder"
    for file_b in ${Folder_B}/*; do
        cd $temp_folder
        temp_file=$(basename $file_b "_")

        echo ${temp_file}
        #echo ${temp_file#*_*_*-*_*}
        echo ${temp_file#*_*_*-*_*_*}
        mv ${temp_file#} ${temp_file#*_*_*-*_*_*_*}
        mv ${temp_file#} ${temp_file#*_*_*-*_*_*}
        mv ${temp_file#} ${temp_file#*_*_*-*_*}

    done
done
