
# Unnecessary

cd /data/mengxu/data

for i in 29 30 31 32 33 34 35 36 37 38 39 ;do
  {

  mkdir /data/mengxu/data/PRJNA773987_count
  mkdir /data/mengxu/data/PRJNA773987_count/SRR166680${i}
  mkdir /data/mengxu/data/PRJNA773987_count/SRR166680${i}/outs
  mkdir /data/mengxu/data/PRJNA773987_count/SRR166680${i}/outs/filtered_feature_bc_matrix
  
  cd /data/mengxu/data/PRJNA773987/SRR166680${i}/SRR166680${i}/outs/filtered_feature_bc_matrix
  cp barcodes.tsv.gz /data/mengxu/data/PRJNA773987_count/SRR166680${i}/outs/filtered_feature_bc_matrix
  cp features.tsv.gz /data/mengxu/data/PRJNA773987_count/SRR166680${i}/outs/filtered_feature_bc_matrix
  cp matrix.mtx.gz /data/mengxu/data/PRJNA773987_count/SRR166680${i}/outs/filtered_feature_bc_matrix
  
  } 
done

