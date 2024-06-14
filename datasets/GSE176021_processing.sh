
#!/bin/sh
path="/data/mengxu/data/GSE176021/"
for file_a in ${path}*
do  

    temp_folder=`basename $file_a  .tar.gz`
    
    echo 1 $temp_folder
    
    if [ ${temp_folder#*_*-*_*_*.*} == vdj ]
    then
      new=`basename $temp_folder  .vdj`
      echo 4 $new
      echo 3 $file_a
      echo 2 ${temp_folder#*_*-*_*_*}
      mv $file_a ${new}_vdj.tar.gz
      gunzip ${new}_vdj.tar.gz
      tar -zxvf ${new}_vdj.tar
      
    else
      tar -zxvf $file_a
      
    fi
    
done

