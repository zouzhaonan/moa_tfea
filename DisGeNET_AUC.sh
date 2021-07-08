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


# 正誤判定
cat chem_tf_diz_es_997.tsv | awktt '
BEGIN {
  while (getline < "200615_curated_chemical_disease.tsv") x[$2,$4]++
} {
  if (x[$1,$3] > 0) print $0,1
  else print $0,0
}'>chem_tf_es_answer.tsv


===============
R
library(ROCR)

  file.name <- "chem_tf_es_answer.tsv"
  rocdata <- read.table (file.name, sep="	")
  pred <- prediction(rocdata[,4], rocdata[,5])
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  auc
  # 0.683486
  
  # PR curve
  file.name="chem_tf_es_answer.tsv"
  rocdata <- read.table (file.name, sep="	")
  pred <- prediction(rocdata[,4], rocdata[,5])
  perf <- performance(pred, "prec", "rec")
  aucpr.tmp <- performance(pred,"aucpr")
  aucpr <- as.numeric(aucpr.tmp@y.values)
  aucpr
  # 0.05741123
===============
cat chem_tf_es_answer.tsv | awktt '$5==1' | wl   # 12069
wl chem_tf_es_answer.tsv       # 553350
# no skill, 12069/553350 = 0.022

cat chem_tf_es_answer.tsv | awktt '
BEGIN {
  while (getline < "200615_curated_chemical_disease.tsv") p[$2,$4]=$6"	"$5
  while (getline < "toMESHID.tsv") d[$4]=$5
  while ("cat ~/ChIP-Atlas/CTD/ChemicalName_ChemicalID.txt" | getline) c[$2]=$1
} {
  print $0,c[$1],d[$3],p[$1,$3]
}' >200807_chem_tf_diz_result.tsv

# bychem (chem-diz) #
mkdir ~/ChIP-Atlas/200807/bychem_diz
for chem in $(awku 1 chem_tf_es_answer.tsv)
do
  cat chem_tf_es_answer.tsv | awk -F"	" -v OFS="	" -v chem="$chem" '$1==chem' >~/ChIP-Atlas/200807/bychem_diz/$chem"_tf_diz.tsv"
done

cd ~/ChIP-Atlas/200807/bychem_diz
# all"1" / all"0" を取り除く #
mkdir ../allposinega_chem_diz
for filename in `ls`
do
  cat $filename | awk -F "	" '!a[$5]++ {print $5}' | echo -e "$filename	`wc -l`"
done | awk -F "	" '$2=="1" {print $1}' >../allposinega_chem_diz/filelist.txt

for filename in `cat ../allposinega_chem_diz/filelist.txt`
do
  mv $filename ../allposinega_chem_diz/
done

rm ~/ChIP-Atlas/200807/auc_temp.txt
=====================
R
library(ROCR)
files <- list.files(pattern=".tsv")
for (file.name in files) {
  rocdata <- read.table (file.name, sep="	")
  pred <- prediction(rocdata[,4], rocdata[,5])
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  str <- paste(auc, file.name,sep="	")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200807/auc_temp.txt", append=TRUE, quote=FALSE, col.names=F)
}

# top_5_tf (es)
cd ../
cat auc_temp.txt | awktt '!a[$0]++{gsub("1 ","",$1);gsub ("_tf_diz.tsv",""); print $2,$1}' >200807_auc_chem_tf_diz_997.tsv

for chem in $(awku 1 200807_auc_chem_tf_diz_997.tsv)
do
  cat 200807_chem_tf_diz_result.tsv | awk -F"	" -v OFS="	" -v chem="$chem" '$1==chem {print $2"("$4")"}' | awktt '!a[$0]++' | head -n 5 | awktt '
  BEGIN {
    ORS="	"
  } {$1=$1;print $0}' | awk -F"	" -v OFS="	" -v chem="$chem" '{print chem,$0}'
done >chem_top_tf5.tsv

# chem_kegglist #
awku 1 200807_auc_chem_tf_diz_997.tsv | awktt '
BEGIN {
  while ("cat ~/ChIP-Atlas/CTD/ChemicalID_KEGGID.txt" | getline) k[$1]=$2
} k[$1]>0 {
  print k[$1]
}' | awku 1 >chem_kegglist.tsv

# KEGGID ごとに ATC を取得
for i in `cat chem_kegglist.tsv`
do
  curl -s "http://rest.kegg.jp/get/"$i | grep -A 1 "Anatomical Therapeutic Chemical (ATC) classification" | awk 'NR == 2' | awk -v i="$i" '{sub(/^[ 	]+/, "")}1 { print i"	"$1 }' # ATC
  sleep 1
done >chem_kegg_atc.tsv

# no of TRUE
for file in `ls ~/ChIP-Atlas/200807/bychem_diz/`
do
  cat ~/ChIP-Atlas/200807/bychem_diz/$file | awktt '$5==1' |echo -e "$file	`wl`"|awktt '{gsub("_tf_diz.tsv",""); print $0}'
done > chem_no_of_true.tsv

# MESH class
cp ~/ChIP-Atlas/200701/chem_meshclass.tsv ./

cat 200807_auc_chem_tf_diz_997.tsv | awktt '
BEGIN {
  while ("cat ~/ChIP-Atlas/CTD/ChemicalName_ChemicalID.txt" | getline) name[$2]=$1
  while ("cat ~/ChIP-Atlas/CTD/ChemicalID_KEGGID.txt" | getline) k[$1]=$2
  while (getline < "chem_kegg_atc.tsv") atc[$1]=$2
  while (getline < "chem_top_tf5.tsv") es[$1]=$0
  while (getline < "chem_no_of_true.tsv") true[$1]=$2
  while (getline < "chem_meshclass.tsv") mesh[$1]=$2
}{
  if (length(k[$1])>0) {
    if (length(atc[k[$1]])>0) print $0,name[$1],es[$1],true[$1],mesh[$1],atc[k[$1]],"https://www.genome.jp/dbget-bin/www_bget?dr_ja:"k[$1]
    else print $0,name[$1],es[$1],true[$1],mesh[$1],"NaN","https://www.genome.jp/dbget-bin/www_bget?dr_ja:"k[$1]
  }
  else print $0,name[$1],es[$1],true[$1],mesh[$1],"NaN","http://ctdbase.org/detail.go?type=chem&acc="$1
}' >chem_tf_diz_bychem_997.tsv
