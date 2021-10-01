##################################
### TFEA の精度評価 (KEGG DRUG) ####
##################################

# chem-TF matrix

mkdir ~/ChIP-Atlas/200805/
cd ~/ChIP-Atlas/200805/

# ChemicalID_TF_ES #

cd ~/ChIP-Atlas/CTD/200515/Results/hg19/
awk -F "	" '{ print FILENAME"	"$3"	"$9 }' *.tsv | awk -F "	" '
BEGIN {
  OFS="	"
}{ 
  gsub ("-", "	", $1)
  gsub (".tsv", "", $1)
  gsub ("-", "", $3)
  print $0
}' | awk -F "	" '{ print $1"	"$4"	"$5 }' | sort -r -n -k3 -t$'	' > ~/ChIP-Atlas/200805/chemID_TF_ES_whole.tsv


cd ~/ChIP-Atlas/200805/
# awkct 1 chemID_TF_ES_whole.tsv     434 chemicals
# awkct 2 chemID_TF_ES_whole.tsv      997 chemicals
# wl chemID_TF_ES_whole.tsv              7122092

cat chemID_TF_ES_whole.tsv|awktt '!a[$1,$2]++' >chemID_TF_ES.tsv
# awkct 1 chemID_TF_ES.tsv           434 chemicals
# awkct 2 chemID_TF_ES.tsv           997 chemicals
# wl chemID_TF_ES.tsv                309744

# TF list
awku 2 chemID_TF_ES.tsv | awktt '{print NR, $0}' >TF_list.tsv
awku 1 chemID_TF_ES.tsv | awktt '{print NR, $0}'>chem_list.tsv

cat chemID_TF_ES.tsv|awktt '
BEGIN {
  while (getline < "TF_list.tsv") no[$2]=$1
}{
  print $1,no[$2],$3
}' >chemID_TFno_ES.tsv

# # local
# mkdir /mnt/c/Users/Zou/"Google Drive"/okilab/ChIP-Atlas/200805/
# cd /mnt/c/Users/Zou/"Google Drive"/okilab/ChIP-Atlas/200805/
# sget
# get /home/zou/ChIP-Atlas/200805/chemID_TF_ES.tsv

split -l 2 --additional-suffix=".txt" -d -a 3 chem_list.tsv chem_temp_

for i in `awk 'BEGIN {for (i=0; i<=216; i++) printf "%03d
", i}'`
do
  qsub -e /dev/null -o /dev/null ~/ChIP-Atlas/200805/chem_tf.sh $i
done

for i in `awk 'BEGIN {for (i=0; i<=216; i++) printf "%03d
", i}'`
do
  cat "chem_tf_matrix_temp_"$i".tsv">>chem_tf_matrix.tsv
done

# 997 TF でそろえる

mkdir ~/ChIP-Atlas/200807/
cd ~/ChIP-Atlas/200807/

cat ~/ChIP-Atlas/200805/chem_tf_matrix.tsv | awktt '{
  for (i=1;i<=997;i++) {
    print $998,i,$i
  }
}' | awktt '
BEGIN {
  while ("cat ~/ChIP-Atlas/200805/TF_list.tsv" | getline) tf[$1]=$2
} {
  print $1,tf[$2],$3
}' | sort -k3 -rn -t$'	' >chem_tf_es_997.tsv

# chem_tf_997 で正誤評価
# iwata_answer を KEGG DRUG に限定する #
cat ~/ChIP-Atlas/200805/iwata_answer.tsv | awktt '$5 ~ /KEGG DRUG/' >kegg_drug_iwata_answer.tsv

cat chem_tf_es_997.tsv | awktt '
BEGIN {
  while (getline < "kegg_drug_iwata_answer.tsv") a[$2,$4]++
} {
  if (length(a[$1,$2]) > 0) print $0,"1"
  else print $0,"0"
}' >mapping_mushi_kegg_drug_997.tsv

===============
R
library(ROCR)

  file.name <- "mapping_mushi_kegg_drug_997.tsv"
  rocdata <- read.table (file.name, sep="	")
  pred <- prediction(rocdata[,3], rocdata[,4])
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  auc
  # 0.6642023
  
  # PR curve
  file.name="mapping_mushi_kegg_drug_997.tsv"
  rocdata <- read.table (file.name, sep="	")
  pred <- prediction(rocdata[,3], rocdata[,4])
  perf <- performance(pred, "prec", "rec")
  aucpr.tmp <- performance(pred,"aucpr")
  aucpr <- as.numeric(aucpr.tmp@y.values)
  aucpr
  # 0.009236002
===============
cat mapping_mushi_kegg_drug_997.tsv | awktt '$4==1' | wl   # 56
wl mapping_mushi_kegg_drug_997.tsv       # 432698
# no skill, 56/432698 = 1.3e-4

# bychem (chem-tf) #
mkdir ~/ChIP-Atlas/200807/bychem
for chem in $(awku 1 mapping_mushi_kegg_drug_997.tsv)
do
  cat mapping_mushi_kegg_drug_997.tsv | awk -F"	" -v OFS="	" -v chem="$chem" '$1==chem' >~/ChIP-Atlas/200807/bychem/$chem"_tf.tsv"
done

cd ~/ChIP-Atlas/200807/bychem

# all"1" / all"0" を取り除く #
mkdir ../allposinega_chem
for filename in `ls`
do
  cat $filename | awk -F "	" '!a[$4]++ {print $4}' | echo -e "$filename	`wc -l`"
done | awk -F "	" '$2=="1" {print $1}' >../allposinega_chem/filelist.txt

for filename in `cat ../allposinega_chem/filelist.txt`
do
  mv $filename ../allposinega_chem/
done

=====================
R
library(ROCR)
files <- list.files(pattern=".tsv")
for (file.name in files) {
  rocdata <- read.table (file.name, sep="	")
  pred <- prediction(rocdata[,3], rocdata[,4])
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  str <- paste(auc, file.name,sep="	")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200807/auc_temp.txt", append=TRUE, quote=FALSE, col.names=F)
  
  perf <- performance(pred, "prec", "rec")
  aucpr.tmp <- performance(pred,"aucpr")
  aucpr <- as.numeric(aucpr.tmp@y.values)
  str <- paste(auc, file.name,sep="	")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200807/aupr_temp_chem_tf.txt", append=TRUE, quote=FALSE, col.names=F)
}

cd ../
cat auc_temp.txt | awktt '!a[$0]++{gsub("1 ","",$1);gsub ("_tf.tsv",""); print $2,$1}' >200807_auc_chem_tf_kegg_drug_997.tsv

for chem in $(awku 1 200807_auc_chem_tf_kegg_drug_997.tsv)
do
  es=$(cat bychem/$chem"_tf.tsv"|awktt '$4==1 {printf "%s	", $2" ("$3")"}' )
  echo -e "$chem	$es"
done >chem_top_es_997.tsv

cat 200807_auc_chem_tf_kegg_drug_997.tsv | awktt '
BEGIN {
  while ("cat ~/ChIP-Atlas/CTD/ChemicalName_ChemicalID.txt" | getline) name[$2]=$1
  while (getline < "chem_top_es_997.tsv") es[$1]=$0
}{
  print $1,name[$1],$2,es[$1]
}' >chem_tf_bychem_997.tsv
