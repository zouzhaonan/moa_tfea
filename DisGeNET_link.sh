##########################################
### DisGeNet につなげる CTD で精度評価 ###
########################################


# DisGeNet
cp ~/ChIP-Atlas/200615_chem_disease/disgenet.tsv ~/ChIP-Atlas/200615_chem_disease/200615_curated_chemical_disease.tsv ~/ChIP-Atlas/200615_chem_disease/toMESHID.tsv ./ # score=$10 #'

genedisease="disgenet.tsv"
cat chem_tf_es_997.tsv | awktt '
BEGIN {
  while (getline < "'$genedisease'") disease[$2]=disease[$2]"@@@@"$5
} length(disease[$2]) > 0 {
  print $1,$2,disease[$2],$3
}' | awktt '{
  count = split ($3, x, "@@@@")
  for (i = 2; i <= count; i++) {
    print $1,$2,x[i],$4
  }
}' |awktt '
BEGIN {
  while (getline < "toMESHID.tsv") mesh[$1]=$4
} mesh[$3]>0 {
  print $1,$2,mesh[$3],$4
}' |sort -k4 -rn -t$'	' | awk '!a[$1,$3]++' > chem_tf_diz_es_997.tsv



