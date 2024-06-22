
#First step: build aggr.csv
#sample_id	molecule_h5
#CRR049227	/data/mengxu/data/CRA001477/CRR049227/CRR049227/outs/molecule_info.h5
#CRR049228	/data/mengxu/data/CRA001477/CRR049228/CRR049228/outs/molecule_info.h5

echo -e "sample_id,molecule_h5 \nCRR049227,/data/mengxu/data/CRA001477/CRR049227/CRR049227/outs/molecule_info.h5 \nCRR049228,/data/mengxu/data/CRA001477/CRR049228/CRR049228/outs/molecule_info.h5" >CRR049227_28.csv

echo -e "sample_id,molecule_h5 \nCRR049227_f,/data/mengxu/data/CRA001477/CRR049227/CRR049227_f/outs/molecule_info.h5 \nCRR049228_f,/data/mengxu/data/CRA001477/CRR049228/CRR049228_f/outs/molecule_info.h5" >CRR049227_28_f.csv

echo -e "sample_id,molecule_h5 \nCRR049229,/data/mengxu/data/CRA001477/CRR049229/CRR049229/outs/molecule_info.h5 \nCRR049230,/data/mengxu/data/CRA001477/CRR049230/CRR049230/outs/molecule_info.h5" >CRR049229_30.csv

echo -e "sample_id,molecule_h5 \nCRR049229_f,/data/mengxu/data/CRA001477/CRR049229/CRR049229_f/outs/molecule_info.h5 \nCRR049230_f,/data/mengxu/data/CRA001477/CRR049230/CRR049230_f/outs/molecule_info.h5" >CRR049229_30_f.csv

#Second step
for i in 27_28 27_28_f 29_30 29_30_f ;do
  {
  
    mkdir /data/mengxu/data/CRA001477/CRR0492${i}
    cd /data/mengxu/data/CRA001477
    cellranger aggr --id=CRR0492${i} \--csv=CRR0492${i}.csv \--normalize=mapped \--localcores=16 \--localmem=64

  }
done

#mkdir /data/mengxu/data/CRA001477/CRR049227_28
#cellranger aggr --id=CRR049227_28 \--csv=CRR049227_28.csv \--normalize=mapped \--localcores=16 \--localmem=64

#mkdir /data/mengxu/data/CRA001477/CRR049227_28_f
#cellranger aggr --id=CRR049227_28_f \--csv=CRR049227_28_f.csv \--normalize=mapped \--localcores=16 \--localmem=64

#mkdir /data/mengxu/data/CRA001477/CRR049229_30
#cellranger aggr --id=CRR049229_30 \--csv=CRR049229_30.csv \--normalize=mapped \--localcores=16 \--localmem=64

#mkdir /data/mengxu/data/CRA001477/CRR049229_30_f
#cellranger aggr --id=CRR049229_30_f \--csv=CRR049229_30_f.csv \--normalize=mapped \--localcores=16 \--localmem=64

