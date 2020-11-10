#############
### TFEA ####
#############

# Chemical_gene_interaction_pattern から #
# AAA results in increased/decreased expression of BBB mRNA #
# AAA results in increased/decreased expression of BBB protein #
# AAA analog results in increased/decreased expression of BBB mRNA #
# という記述のみを抽出 #

# STEP. 1/7
========================
DATE="200515"

# スパコン 
qsub -l s_vmem=16G -l mem_req=16G -e /dev/null -o /dev/null ~/ChIP-Atlas/CTD/Chemical_Gene_list.sh $DATE

## nano ~/ChIP-Atlas/CTD/Chemical_Gene_list.sh
#============================
##!/bin/sh
##$ -S /bin/sh
#
#DATE=$1
#
#mkdir ~/ChIP-Atlas/CTD/$DATE
#cd ~/ChIP-Atlas/CTD/$DATE
## CTD_chem_gene_ixns.tsv
#curl http://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz | gunzip >CTD_chem_gene_ixns.tsv
## ChemicalName  ChemicalID  CasRN  GeneSymbol  GeneID  GeneForms  Organism  OrganismID  Interaction  InteractionActions  PubMedIDs
#
## increased/decreased expression && 欠損値がない record を抽出#
## ChemicalName, ChemicalID, PubMedID, OrganismID, GeneSymbol, Up_Down #
#cat CTD_chem_gene_ixns.tsv | grep -v ^# | awk -F "\t" '{
#  ChemicalName = $1
#  ChemicalID = $2
#  PubMedID = $11
#  OrganismID = $8
#  GeneSymbol = $4
#  if ($9 == $1" results in increased expression of "$4" mRNA") Up_Down = "Up"
#  else if ($9 == $1" results in decreased expression of "$4" mRNA") Up_Down = "Down"
#  else if ($9 == $1" results in increased expression of "$4" protein") Up_Down = "Up"
#  else if ($9 == $1" results in decreased expression of "$4" protein") Up_Down = "Down"
#  else if ($9 == $1" analog results in increased expression of "$4" mRNA") Up_Down = "Up"
#  else if ($9 == $1" analog results in decreased expression of "$4" mRNA") Up_Down = "Down"
#  else Up_Down = NULL
#  printf "%s\t%s\t%s\t%s\t%s\t%s\n", ChemicalName, ChemicalID, PubMedID, OrganismID, GeneSymbol, Up_Down
#}' | awk -F "\t" '$1 && $2 && $3 && $4 && $5 && $6' | sort | uniq >chem_gene_change_temp1.tsv
#
## PubMedID "|" 処理 #
#cat chem_gene_change_temp1.tsv | awk -F "\t" '{
#  count = split ($3, PubMedID, /[|]/)
#  for (i = 1; i <= count; i++) print $1"\t"$2"\t"PubMedID[i]"\t"$4"\t"$5"\t"$6
#}' | sort | uniq >chem_gene_change.tsv
#
#mkdir ~/ChIP-Atlas/CTD/$DATE/Up
#mkdir ~/ChIP-Atlas/CTD/$DATE/Down
#
## Up or Down #
#cd ~/ChIP-Atlas/CTD/$DATE
#cat chem_gene_change.tsv | awk -F "\t" '{
#  if ($6 == "Up") print $0 >>"Up/chem_gene_Up.tsv"
#  else print $0 >>"Down/chem_gene_Down.tsv"
#}'
#
## Organism ごとに #
#for CHANGE in Up Down
#do
#  # フォルダ #
#  for ORGANISM in hg19 mm9 other
#  do
#    mkdir ~/ChIP-Atlas/CTD/$DATE/$CHANGE/$ORGANISM
#  done
#  # genelist #
#  cd ~/ChIP-Atlas/CTD/$DATE/$CHANGE/
#  cat chem_gene_$CHANGE.tsv | awk -F "\t" '{
#    ChemicalID = $2
#    PubMedID = $3
#    OrganismID = $4
#    GeneSymbol = $5
#    filename = ChemicalID"-"PubMedID"-"OrganismID".txt"
#    if (OrganismID == "9606") print GeneSymbol >>"hg19/"filename
#    else if (OrganismID == "10090") print GeneSymbol >>"mm9/"filename
#    else print GeneSymbol >>"other/"filename
#  }'
#done
#===============



# insilicoChIP に使える genelist を準備する #
# 1. gene 数 10 以上 #
# 2. Up/Down 両方のデータある #

# STEP. 2/7
========================
DATE="200515"
mkdir -p ~/tmp/    # insilicoChIP temporary #
for ORGANISM in hg19 mm9
do
  for CHANGE in Up Down
  do
    # スパコン #
    # gene 数が 10 未満のやつを $CHANGE/$ORGANISM/genes_under_ten/ に移動 #
    qsub -l s_vmem=16G -l mem_req=16G -e /dev/null -o /dev/null ~/ChIP-Atlas/CTD/Genelist_process_for_TFEA.sh $DATE $CHANGE $ORGANISM
  done
done

##nano ~/ChIP-Atlas/CTD/Genelist_process_for_TFEA.sh
#====================
##!/bin/sh
##$ -S /bin/sh
#
#DATE=$1 
#CHANGE=$2
#ORGANISM=$3
## gene 数が 10 未満のやつを $CHANGE/$ORGANISM/genes_under_ten/ に移動 #
#mkdir ~/ChIP-Atlas/CTD/$DATE/$CHANGE/$ORGANISM/genes_under_ten/
#cd ~/ChIP-Atlas/CTD/$DATE/$CHANGE/$ORGANISM/
#for filename in `ls | grep "\.txt$"`
#do
#  wc -l $filename
#done >~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_"$CHANGE"_temp_wc.txt"
#cd ~/ChIP-Atlas/CTD/$DATE/
#cat "filelist_"$ORGANISM"_"$CHANGE"_temp_wc.txt" | awk -F " " ' $1 < 10 { print $2 }' >"filelist_"$ORGANISM"_"$CHANGE"_temp_under_ten.txt"
#for filename in `cat "filelist_"$ORGANISM"_"$CHANGE"_temp_under_ten.txt"`
#do
#  mv $CHANGE/$ORGANISM/$filename $CHANGE/$ORGANISM/genes_under_ten/
#done
## gene 数が 10 以上のやつをリストアップ #
#ls $CHANGE/$ORGANISM/ | grep "\.txt$" >~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_"$CHANGE"_temp_above_ten.txt"
#==========================

# STEP. 3/7
=====================
for ORGANISM in hg19 mm9
do
  # Up/Down 両方のデータがある record を抽出 #
  cd ~/ChIP-Atlas/CTD/$DATE/
  cat "filelist_"$ORGANISM"_Up_temp_above_ten.txt" | awk -F "\t" -v ORGANISM="$ORGANISM" '    # 外部変数の渡し方に注意！！ #
    BEGIN {
      filename="filelist_"ORGANISM"_Down_temp_above_ten.txt"    # ここでは "$" いらない #
      while ((getline < filename) > 0)
        x[$1] = $1
    }
    length(x[$1]) > 0
  ' >"filelist_"$ORGANISM"_common_temp_above_ten.txt"
  # filename.txt -> filename #
  cat ~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_common_temp_above_ten.txt" | awk '{
    gsub (/\.txt$/, "", $1)
    print $0
  }' >~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_common_above_ten.txt"
  # TFEA 出力フォルダ #
  mkdir -p ~/ChIP-Atlas/CTD/$DATE/Results/$ORGANISM/
done

# STEP. 4/7
==========================
# TFEA #
# 出力は pval でソート済み #
cd ~/ChIP-Atlas/CTD/$DATE/
# job 数調整 #
for ORGANISM in hg19 mm9
do
  wc -l "filelist_"$ORGANISM"_common_above_ten.txt"
done
# 900 hg19
# 607 mm9

# STEP. 5/7
==================
# hg19_450*2_groups #
split -l 450 --additional-suffix=".txt" -d ~/ChIP-Atlas/CTD/$DATE/filelist_hg19_common_above_ten.txt filelist_hg19_common_above_ten_
# mm9_305*2_groups #
split -l 305 --additional-suffix=".txt" -d ~/ChIP-Atlas/CTD/$DATE/filelist_mm9_common_above_ten.txt filelist_mm9_common_above_ten_
# 合計 4 groups に分けて qsub #

# STEP. 6/7
==================
# GROUP. 1
ORGANISM="hg19"
for filename in `cat ~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_common_above_ten_00.txt"`    # index = "00"
do
  qsub -l s_vmem=16G -l mem_req=16G -e /dev/null -o /dev/null ~/bin/insilicoChIP -R -a ~/ChIP-Atlas/CTD/$DATE/Up/$ORGANISM/$filename.txt -b ~/ChIP-Atlas/CTD/$DATE/Down/$ORGANISM/$filename.txt gene $ORGANISM ~/ChIP-Atlas/CTD/$DATE/Results/$ORGANISM/$filename
done

# GROUP. 2
ORGANISM="hg19"
for filename in `cat ~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_common_above_ten_01.txt"`    # index = "01"
do
  qsub -l s_vmem=16G -l mem_req=16G -e /dev/null -o /dev/null ~/bin/insilicoChIP -R -a ~/ChIP-Atlas/CTD/$DATE/Up/$ORGANISM/$filename.txt -b ~/ChIP-Atlas/CTD/$DATE/Down/$ORGANISM/$filename.txt gene $ORGANISM ~/ChIP-Atlas/CTD/$DATE/Results/$ORGANISM/$filename
done

# GROUP. 3
ORGANISM="mm9"
for filename in `cat ~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_common_above_ten_00.txt"`    # index = "00"
do
  qsub -l s_vmem=16G -l mem_req=16G -e /dev/null -o /dev/null ~/bin/insilicoChIP -R -a ~/ChIP-Atlas/CTD/$DATE/Up/$ORGANISM/$filename.txt -b ~/ChIP-Atlas/CTD/$DATE/Down/$ORGANISM/$filename.txt gene $ORGANISM ~/ChIP-Atlas/CTD/$DATE/Results/$ORGANISM/$filename
done

# GROUP. 4
ORGANISM="mm9"
for filename in `cat ~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_common_above_ten_01.txt"`    # index = "01"
do
  qsub -l s_vmem=16G -l mem_req=16G -e /dev/null -o /dev/null ~/bin/insilicoChIP -R -a ~/ChIP-Atlas/CTD/$DATE/Up/$ORGANISM/$filename.txt -b ~/ChIP-Atlas/CTD/$DATE/Down/$ORGANISM/$filename.txt gene $ORGANISM ~/ChIP-Atlas/CTD/$DATE/Results/$ORGANISM/$filename
done

# 以上で TFEA 終了 #

# STEP. 7/7
======================
# summary 作成 #

# ChemicalName ChemicalID 辞書 #
cd ~/ChIP-Atlas/CTD/$DATE/
cat CTD_chem_gene_ixns.tsv | grep -v ^# | awk -F "\t" '{ print $1"\t"$2 }' | sort | uniq >ChemicalName_ID_dictionary.tsv

for ORGANISM in hg19 mm9
do
  # ChIP-Atlas は microRNA に現時点で対応していないので，出力エラーを取り除く． #
  # 成功したもののリスト #
  cd ~/ChIP-Atlas/CTD/$DATE/Results/$ORGANISM
  rm *.bed*
  rm *.tmpForinsilicoChIP*
  for filename in `ls | grep "\.tsv$"`
  do
    wc -l $filename 
  done > ../$ORGANISM"_wc_temp.txt"
  cd ~/ChIP-Atlas/CTD/$DATE/Results
  cat $ORGANISM"_wc_temp.txt" | awk -F " " ' $1 > 0 {gsub (/[ \.]/, "\t"); print $2}' >$ORGANISM"_TEFA_done_list_temp.txt"
  
  # 失敗したもののリスト #
  cat ../"filelist_"$ORGANISM"_common_above_ten.txt" | awk -F "\t" -v ORGANISM="$ORGANISM" '
  BEGIN {
    Donelist=ORGANISM"_TEFA_done_list_temp.txt"
    while ((getline < Donelist) > 0)
      x[$1] = $1
  }
  length(x[$1]) == 0
  ' >$ORGANISM"_TEFA_failed_list_temp.txt"
  
  # insilicoChIP の top_of_pval #
  # FE > 1 TRUE/FALSE でそれぞれまとめる #
  cd ~/ChIP-Atlas/CTD/$DATE/Results/$ORGANISM
  for filename in `cat ../$ORGANISM"_TEFA_done_list_temp.txt"`
  do
    cat $filename.tsv | awk -F "\t" -v ORGANISM="$ORGANISM" -v filename="$filename" '{
      if ($11 > 1) {
        True="../"ORGANISM"_summary_TRUE_temp.tsv"
        print filename"\t"$0 >>True
        exit
      }
    }'
    cat $filename.tsv | awk -F "\t" -v ORGANISM="$ORGANISM" -v filename="$filename" '{
      if ($11 <= 1) {
        False="../"ORGANISM"_summary_FALSE_temp.tsv"
        print filename"\t"$0 >>False
        exit
      }
    }'
  done
    
  # ChemicalName を追加 #
  cd ~/ChIP-Atlas/CTD/$DATE/Results/
  for FE in TRUE FALSE
  do
    cat $ORGANISM"_summary_"$FE"_temp.tsv" | awk -F "\t" '
    BEGIN{
      OFS="\t"
    }{
      gsub ("-", "\t", $1)
      print $0
    }' | awk -F "\t" '
    BEGIN {
      OFS="\t"
      while ((getline < "../ChemicalName_ID_dictionary.tsv") > 0)
      Name_of[$2] = $1
    }{
      print Name_of[$1]"\t"$0
    }' > $ORGANISM"_summary_"$FE".tsv"
  done
done



##################################
### TFEA の精度評価 (KEGG DRUG) ####
##################################

# chem-TF matrix

mkdir ~/ChIP-Atlas/200805/
cd ~/ChIP-Atlas/200805/

# ChemicalID_TF_ES #

cd ~/ChIP-Atlas/CTD/200515/Results/hg19/
awk -F "\t" '{ print FILENAME"\t"$3"\t"$9 }' *.tsv | awk -F "\t" '
BEGIN {
  OFS="\t"
}{ 
  gsub ("-", "\t", $1)
  gsub (".tsv", "", $1)
  gsub ("-", "", $3)
  print $0
}' | awk -F "\t" '{ print $1"\t"$4"\t"$5 }' | sort -r -n -k3 -t$'\t' > ~/ChIP-Atlas/200805/chemID_TF_ES_whole.tsv


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

for i in `awk 'BEGIN {for (i=0; i<=216; i++) printf "%03d\n", i}'`
do
  qsub -e /dev/null -o /dev/null ~/ChIP-Atlas/200805/chem_tf.sh $i
done
 
  


# nano ~/ChIP-Atlas/200805/chem_tf.sh
# ==============================================
# #!/bin/sh
# #$ -S /bin/sh
# label=$1
# cd ~/ChIP-Atlas/200805/
# for chem in `awku 2 "chem_temp_"$label".txt"`
# do
#   cat ~/ChIP-Atlas/200805/chemID_TFno_ES.tsv | awk -F "\t" -v OFS="\t" -v chem="$chem" '$1==chem' >$chem"_temp.tsv"
#   for ((i=1;i<=997;i++))
#   do
#     cat $chem"_temp.tsv" | awk -F "\t" -v OFS="\t" -v i="$i" -v chem="$chem" '$1==chem {
#       if ($2==i) printf "%s\t",$3
#     }'
#     a=`cat $chem"_temp.tsv" | awk -F "\t" -v OFS="\t" -v i="$i" -v chem="$chem" '$1==chem {
#       if ($2==i) printf "%s\t",$3
#     }'`
#     if [[ $a == "" ]]; then
#       echo -n -e "0\t"
#     fi
#     if [[ $i == 997 ]]; then
#       echo $chem
#     fi
#   done
# done >"chem_tf_matrix_temp_"$label".tsv"
# ===============================================

for i in `awk 'BEGIN {for (i=0; i<=216; i++) printf "%03d\n", i}'`
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
}' | sort -k3 -rn -t$'\t' >chem_tf_es_997.tsv

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
  rocdata <- read.table (file.name, sep="\t")
  pred <- prediction(rocdata[,3], rocdata[,4])
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  auc
  # 0.6642023
  
  # PR curve
  file.name="mapping_mushi_kegg_drug_997.tsv"
  rocdata <- read.table (file.name, sep="\t")
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
  cat mapping_mushi_kegg_drug_997.tsv | awk -F"\t" -v OFS="\t" -v chem="$chem" '$1==chem' >~/ChIP-Atlas/200807/bychem/$chem"_tf.tsv"
done

cd ~/ChIP-Atlas/200807/bychem

# all"1" / all"0" を取り除く #
mkdir ../allposinega_chem
for filename in `ls`
do
  cat $filename | awk -F "\t" '!a[$4]++ {print $4}' | echo -e "$filename\t`wc -l`"
done | awk -F "\t" '$2=="1" {print $1}' >../allposinega_chem/filelist.txt

for filename in `cat ../allposinega_chem/filelist.txt`
do
  mv $filename ../allposinega_chem/
done

=====================
R
library(ROCR)
files <- list.files(pattern=".tsv")
for (file.name in files) {
  rocdata <- read.table (file.name, sep="\t")
  pred <- prediction(rocdata[,3], rocdata[,4])
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  str <- paste(auc, file.name,sep="\t")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200807/auc_temp.txt", append=TRUE, quote=FALSE, col.names=F)
}

cd ../
cat auc_temp.txt | awktt '!a[$0]++{gsub("1 ","",$1);gsub ("_tf.tsv",""); print $2,$1}' >200807_auc_chem_tf_kegg_drug_997.tsv

for chem in $(awku 1 200807_auc_chem_tf_kegg_drug_997.tsv)
do
  es=$(cat bychem/$chem"_tf.tsv"|awktt '$4==1 {printf "%s\t", $2" ("$3")"}' )
  echo -e "$chem\t$es"
done >chem_top_es_997.tsv

cat 200807_auc_chem_tf_kegg_drug_997.tsv | awktt '
BEGIN {
  while ("cat ~/ChIP-Atlas/CTD/ChemicalName_ChemicalID.txt" | getline) name[$2]=$1
  while (getline < "chem_top_es_997.tsv") es[$1]=$0
}{
  print $1,name[$1],$2,es[$1]
}' >chem_tf_bychem_997.tsv

================================
# local
mkdir /mnt/c/Users/Zou/"Google Drive"/okilab/ChIP-Atlas/200807/
cd /mnt/c/Users/Zou/"Google Drive"/okilab/ChIP-Atlas/200807/
sget
get /home/zou/ChIP-Atlas/200807/chem_tf_bychem_997.tsv
get /home/zou/ChIP-Atlas/200807/mapping_mushi_kegg_drug_997.tsv
===============================



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
}' |sort -k4 -rn -t$'\t' | awk '!a[$1,$3]++' > chem_tf_diz_es_997.tsv


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
  rocdata <- read.table (file.name, sep="\t")
  pred <- prediction(rocdata[,4], rocdata[,5])
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  auc
  # 0.683486
  
  # PR curve
  file.name="chem_tf_es_answer.tsv"
  rocdata <- read.table (file.name, sep="\t")
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
  while (getline < "200615_curated_chemical_disease.tsv") p[$2,$4]=$6"\t"$5
  while (getline < "toMESHID.tsv") d[$4]=$5
  while ("cat ~/ChIP-Atlas/CTD/ChemicalName_ChemicalID.txt" | getline) c[$2]=$1
} {
  print $0,c[$1],d[$3],p[$1,$3]
}' >200807_chem_tf_diz_result.tsv

# bychem (chem-diz) #
mkdir ~/ChIP-Atlas/200807/bychem_diz
for chem in $(awku 1 chem_tf_es_answer.tsv)
do
  cat chem_tf_es_answer.tsv | awk -F"\t" -v OFS="\t" -v chem="$chem" '$1==chem' >~/ChIP-Atlas/200807/bychem_diz/$chem"_tf_diz.tsv"
done

cd ~/ChIP-Atlas/200807/bychem_diz
# all"1" / all"0" を取り除く #
mkdir ../allposinega_chem_diz
for filename in `ls`
do
  cat $filename | awk -F "\t" '!a[$5]++ {print $5}' | echo -e "$filename\t`wc -l`"
done | awk -F "\t" '$2=="1" {print $1}' >../allposinega_chem_diz/filelist.txt

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
  rocdata <- read.table (file.name, sep="\t")
  pred <- prediction(rocdata[,4], rocdata[,5])
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  str <- paste(auc, file.name,sep="\t")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200807/auc_temp.txt", append=TRUE, quote=FALSE, col.names=F)
}

# top_5_tf (es)
cd ../
cat auc_temp.txt | awktt '!a[$0]++{gsub("1 ","",$1);gsub ("_tf_diz.tsv",""); print $2,$1}' >200807_auc_chem_tf_diz_997.tsv

for chem in $(awku 1 200807_auc_chem_tf_diz_997.tsv)
do
  cat 200807_chem_tf_diz_result.tsv | awk -F"\t" -v OFS="\t" -v chem="$chem" '$1==chem {print $2"("$4")"}' | awktt '!a[$0]++' | head -n 5 | awktt '
  BEGIN {
    ORS="\t"
  } {$1=$1;print $0}' | awk -F"\t" -v OFS="\t" -v chem="$chem" '{print chem,$0}'
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
  curl -s "http://rest.kegg.jp/get/"$i | grep -A 1 "Anatomical Therapeutic Chemical (ATC) classification" | awk 'NR == 2' | awk -v i="$i" '{sub(/^[ \t]+/, "")}1 { print i"\t"$1 }' # ATC
  sleep 1
done >chem_kegg_atc.tsv

# no of TRUE
for file in `ls ~/ChIP-Atlas/200807/bychem_diz/`
do
  cat ~/ChIP-Atlas/200807/bychem_diz/$file | awktt '$5==1' |echo -e "$file\t`wl`"|awktt '{gsub("_tf_diz.tsv",""); print $0}'
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



#############################################
#### DEG-pivoted prediction CTD で精度評価######
##############################################

# CREEDS #
mkdir ~/ChIP-Atlas/creeds
cd ~/ChIP-Atlas/creeds
curl http://amp.pharm.mssm.edu/CREEDS/download/disease_signatures-v1.0.csv >disease_signatures-v1.0.csv  # metadata
curl http://amp.pharm.mssm.edu/CREEDS/download/disease_signatures-v1.0.json >disease_signatures-v1.0.json # disease-DEG

csv2tsv disease_signatures-v1.0.csv >~/ChIP-Atlas/creeds/disease_signatures-v1.0.tsv

# metadata
# 1      2        3         4         5           6       7       8         9        10         11         12
# id  cell_type  ctrl_ids  curator  disease_name  do_id  geo_id  organism  pert_ids  platform   umls_cui    version

cat ~/ChIP-Atlas/creeds/disease_signatures-v1.0.tsv|awktt '{
  gsub (":","",$6)
  if (length($11) > 0) print $5,$11
  else print $5,$6
}' >~/ChIP-Atlas/creeds/Name_ID.tsv

# disease-DEG (json) to tsv
# 1             2       3       4
# umlsID/DOID   dzID    gene    up/down (human 限定)
cat disease_signatures-v1.0.json | awktt '{
  gsub ("{"do_id"","\n{"do_id"")
  gsub ("\"down_genes\"","\n\"down_genes\"")
  gsub ("\"up_genes\"","\n\"up_genes\"")
  gsub ("\"disease_name\"","\n\"disease_name\"")
  gsub ("\"curator\"","\n\"curator\"")
  print $0
}' | awktt 'NR>1'| awktt '{
  a[NR]=$0
  if (NR%5==4) printf "%s\t%s\n%s\t%s\n", a[NR], a[NR-2], a[NR], a[NR-1]
}' | awktt '{
  gsub ("\"","\t",$1)
  gsub ("down_genes","down",$2)
  gsub ("up_genes","up",$2)
  gsub (": ","\t",$2)
  gsub ("\"","",$2)
  gsub ("\\]\\]","]",$2)
  gsub ("\\[","",$2)
  print $0
}' | awktt '$8 == "human" { print $12,$8,$14,$15 }'| awktt '{
  count=split($4,x,"], ")
  for (i=1;i<count;i++) print $1,$2,$3,x[i]
}' | awktt '{ gsub(",","\t",$4); print $0 }' | awktt '{print $1,$4,$3}' | awktt '
BEGIN {
  while (getline < "disease_signatures-v1.0.tsv") {
    umls[$1]=$11
    DOID[$1]=$6
  }
} {
  if (length(umls[$1]) > 0) print umls[$1],$0
  else print DOID[$1],$0
}' >200616_creeds_disease_gene.tsv  # 555 行 # 
# disease 236

mkdir ~/ChIP-Atlas/200707/
cd ~/ChIP-Atlas/200707/

chem-gene  ~/ChIP-Atlas/CTD/200515/chem_gene_change.tsv
#  1 ChemName   10074-G5          
#  2 ChemID     C534883 
#  3 Pubmed     26036281   
#  4 Organism   9606 
#  5 gene       MYC  
#  6 change     Down

diz-gene   ~/ChIP-Atlas/creeds/200616_creeds_disease_gene.tsv
#  1 Dis(creeds)ID   C0028043
#  2 creedno.        dz:325
#  3 gene            BEX1
#  4 change          down

cat ~/ChIP-Atlas/CTD/200515/chem_gene_change.tsv | awktt '$4=="9606" {print $2"_"$3,$5,$6}' > chem_gene_original.tsv
# 1                       2       3
# C534883_26036281        MYC     Down

for i in $(awku 1 chem_gene_original.tsv)
do
  id=$i
  cat chem_gene_original.tsv | awk -F"\t" -v OFS="\t" -v i="$i" '
  BEGIN {
    up=0
    down=0
  } {
    if ($1==i && $3 == "Up") up++
    else if ($1==i && $3 == "Down") down++
  }
  END {
    print i,up,down
  }'
done >chem_gene_count.tsv
# id                      Up      Down
# C534883_26036281        0       1

cat ~/ChIP-Atlas/creeds/200616_creeds_disease_gene.tsv | awktt '{gsub(":",""); print $1"_"$2,$3,$4}' > diz_gene_original.tsv
# 1                   2       3
# C0028043_dz325      BEX1    down

for i in $(awku 1 diz_gene_original.tsv)
do
  id=$i
  cat diz_gene_original.tsv | awk -F"\t" -v OFS="\t" -v i="$i" '
  BEGIN {
    up=0
    down=0
  } {
    if ($1==i && $3 == "up") up++
    else if ($1==i && $3 == "down") down++
  }
  END {
    print i,up,down
  }'
done >diz_gene_count.tsv

cat chem_gene_count.tsv | awktt '$2>=10 && $3>=10 {print $1}'>chemlist_10.tsv   # 900 #
cat diz_gene_count.tsv | awktt '$2>=10 && $3>=10 {print $1}'>dizlist_10.tsv     # 553 #

mkdir chemgene
mkdir dizgene

for i in $(cat chemlist_10.tsv)
do
  cat chem_gene_original.tsv | awk -F"\t" -v OFS="\t" -v i="$i" '$1==i' >chemgene/$i".tsv"
done

for i in $(cat dizlist_10.tsv)
do
  cat diz_gene_original.tsv | awk -F"\t" -v OFS="\t" -v i="$i" '$1==i' >dizgene/$i".tsv"
done

curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz | gunzip> refFlat.txt
awku 1 refFlat.txt |sort |awktt '{print NR,$1}'>refseqgene.tsv

# chem/diz gene profile
# Up/Down/others

mkdir chemgeneprof
mkdir dizgeneprof

for i in $(cat chemlist_10.tsv)
do
  cat refseqgene.tsv | awk -F"\t" -v OFS="\t" -v i="$i" '
  BEGIN {
    file="chemgene/"i".tsv"
    while (getline < file) change[$2]=$3
  } {
    if (change[$2]=="Up") print $0,"up"
    else if (change[$2]=="Down") print $0,"down"
    else print $0,"NaN"
  }' > chemgeneprof/$i"_prof.tsv"
done

for i in $(cat dizlist_10.tsv)
do
  cat refseqgene.tsv | awk -F"\t" -v OFS="\t" -v i="$i" '
  BEGIN {
    file="dizgene/"i".tsv"
    while (getline < file) change[$2]=$3
  } {
    if (change[$2]=="up") print $0,"up"
    else if (change[$2]=="down") print $0,"down"
    else print $0,"NaN"
  }' > dizgeneprof/$i"_prof.tsv"
done

# #               D1
# #           up  down  NaN
# #     up    a   b     c
# # C1  down  d   e     f
# #     NaN   g   h     i

# fisher.test, up vs down #
=======================================
R
files <- list.files(pattern="tsv")
for (file.name in files) {
  data0=read.table(file.name, sep="\t")
  data=data0[1:2,1:2]
  fisher <- fisher.test(data)
  str <- paste(file.name, fisher$p.value, sep="\t")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200707/up_vs_down_result_temp.tsv", append=TRUE, quote=F, col.names=F,)
}
=======================================

# fisher.test, deg vs non-deg #
=======================================
R
files <- list.files(pattern="tsv")
for (file.name in files) {
  data0=read.table(file.name, sep="\t")
  a=data0[1,1]+data0[1,2]+data0[2,1]+data0[2,2]
  b=data0[1,3]+data0[2,3]
  c=data0[3,1]+data0[3,2]
  d=data0[3,3]
  data=matrix(c(a, b, c, d), nrow = 2, ncol = 2)  
  fisher <- fisher.test(data)
  str <- paste(file.name, fisher$p.value, sep="\t")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200707/degvsnondeg_fisher_result_temp.tsv", append=TRUE, quote=F, col.names=F,)
}
=======================================

cd ~/ChIP-Atlas/200707/
cat degvsnondeg_fisher_result_temp.tsv | awktt '{gsub("1 ","",$1);print $1,-log($2)/log(10)}'|sort -n -r -k2 > degvsnondeg_fisher_result.tsv
cat up_vs_down_result_temp.tsv | awktt '{gsub("1 ","",$1);print $1,-log($2)/log(10)}'|sort -n -r -k2 > updown_fisher_result.tsv

answer="200615_curated_chemical_disease.tsv"


for i in updown_fisher
do
  cat $i"_result.tsv" | awktt '{gsub("_x_","\t");gsub("_","\t");gsub(".tsv","");print $0}'| awktt '
  BEGIN {
    while ("cat ~/ChIP-Atlas/200617/DO2umlsID.tsv"| getline) umls["DOID"$1] = $2
    while ("cat ~/ChIP-Atlas/200617/umls2MESHID.tsv"| getline) mesh[$1] = $2
  } {
    if (length(umls[$3]) > 0) print $1,$2,mesh[umls[$3]],$4,$5
    else print $1,$2,mesh[$3],$4,$5
  }' | awktt '$3 != NULL' | awktt '!a[$1,$3]++ {print $1,$3,$5}' | awktt '
  BEGIN {
    while (getline < "'$answer'") x[$2,$4]++
  } {
    if (x[$1,$2] > 0) print $0,1
    else print $0,0
  }' >roc/$i".tsv"
done
cd ~/ChIP-Atlas/200707/roc
========================
R
library(ROCR)

file.name="updown_fisher.tsv"
rocdata <- read.table (file.name, sep="\t")
pred <- prediction(rocdata[,3], rocdata[,4])
perf <- performance(pred, "tpr", "fpr")
auc.tmp <- performance(pred,"auc")
auc <- as.numeric(auc.tmp@y.values)
auc
# AUC = 0.6413368
perf <- performance(pred, "prec", "rec")
aucpr.tmp <- performance(pred,"aucpr")
aucpr <- as.numeric(aucpr.tmp@y.values)
aucpr
# AUPR = 0.05048038
========================

for i in degvsnondeg_fisher
do
  cat $i"_result.tsv" | awktt '{gsub("_x_","\t");gsub("_","\t");gsub(".tsv","");print $0}'| awktt '
  BEGIN {
    while ("cat ~/ChIP-Atlas/200617/DO2umlsID.tsv"| getline) umls["DOID"$1] = $2
    while ("cat ~/ChIP-Atlas/200617/umls2MESHID.tsv"| getline) mesh[$1] = $2
  } {
    if (length(umls[$3]) > 0) print $1,$2,mesh[umls[$3]],$4,$5
    else print $1,$2,mesh[$3],$4,$5
  }' | awktt '$3 != NULL' | awktt '!a[$1,$3]++ {print $1,$3,$5}' | awktt '
  BEGIN {
    while (getline < "'$answer'") x[$2,$4]++
  } {
    if (x[$1,$2] > 0) print $0,1
    else print $0,0
  }' >roc/$i".tsv"
done

cd ~/ChIP-Atlas/200707/roc
========================
R
library(ROCR)

file.name="degvsnondeg_fisher.tsv"
rocdata <- read.table (file.name, sep="\t")
pred <- prediction(rocdata[,7], rocdata[,8])
perf <- performance(pred, "tpr", "fpr")
auc.tmp <- performance(pred,"auc")
auc <- as.numeric(auc.tmp@y.values)
# AUC = 0.6285592
perf <- performance(pred, "prec", "rec")
aucpr.tmp <- performance(pred,"aucpr")
aucpr <- as.numeric(aucpr.tmp@y.values)
aucpr
# AUPR = 0.04606771
========================


###############################################
# CTD, CREEDS, L1000 に入っている遺伝子の比較 #
###############################################

mkdir ~/ChIP-Atlas/200820
cd ~/ChIP-Atlas/200820

cat ~/ChIP-Atlas/200522_TFEA_ROC_pubmed/TFEA_whole_pubmed.tsv | awktt '{print $1"_"$2}' | awktt '!a[$0]++' >chem_pubmed.tsv

cat ~/ChIP-Atlas/CTD/200515/chem_gene_change.tsv | awktt '
BEGIN {
  while ("cat chem_pubmed.tsv" | getline) a[$0]++
} length(a[$2"_"$3]) > 0
' | awku 5 >CTD_TFEA_genelist_all.tsv

curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz |gunzip > lincs_gene.tsv

curl ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt >hgnc_complete_set.txt

cat hgnc_complete_set.txt|awktt '{print $2,$9,$11}'| awktt '{
  gsub ("\"","")
#   gsub ("\|","\t")
  print $1,$2"|"$3
}' | awktt '{
  c=split($2,x,"|")
  for (i=1;i<=c;i++) print $1,x[i]
}'| awktt 'NR>=3 && $2!=""' >HCSC_gene_synonym.tsv

cat lincs_gene.tsv | awktt '
BEGIN {
  while ("cat ~/HCSC_gene_synonym.tsv" | getline) a[$2]=$1
} NR > 1 && $4==1 {
  if(length(a[$2])>0) print a[$2]
  else print $2
}' > lincs_genelist.tsv

cat ~/ChIP-Atlas/200707/refseqgene.tsv | awktt '
BEGIN {
  while (getline < "CTD_TFEA_genelist_all.tsv") ctd[$0]++
  while (getline < "lincs_genelist.tsv") lincs[$0]++
}{
  if (ctd[$2] > 0 && lincs[$2] > 0) print $0,"1","1"
  else if (ctd[$2] > 0) print $0,"1","0"
  else if (lincs[$2] > 0) print $0,"0","1"
  else print $0,"0","0"
}' >refseq_ctd_lincs_gene.tsv

cat CTD_TFEA_genelist_all.tsv | awktt '
BEGIN {
  while ("cat ~/ChIP-Atlas/200707/refseqgene.tsv" | getline) ref[$2]++
  while (getline < "lincs_genelist.tsv") lincs[$0]++
}{
  if (ref[$1] > 0 && lincs[$1] > 0) print $0,"1","1"
  else if (ref[$1] > 0) print $0,"1","0"
  else if (lincs[$1] > 0) print $0,"0","1"
  else print $0,"0","0"
}' >ctd_refseq_lincs_gene.tsv

cat lincs_genelist.tsv | awktt '
BEGIN {
  while ("cat ~/ChIP-Atlas/200707/refseqgene.tsv" | getline) ref[$2]++
  while (getline < "CTD_TFEA_genelist_all.tsv") ctd[$0]++
}{
  if (ref[$1] > 0 && ctd[$1] > 0) print $0,"1","1"
  else if (ref[$1] > 0) print $0,"1","0"
  else if (ctd[$1] > 0) print $0,"0","1"
  else print $0,"0","0"
}' >lincs_refseq_ctd_gene.tsv

cat ~/ChIP-Atlas/CTD/200515/chem_gene_change.tsv | awktt '
BEGIN {
  while ("cat chem_pubmed.tsv" | getline) a[$0]++
} length(a[$2"_"$3]) > 0
' |awktt '{print $5}'>CTD_TFEA_genelist_not_unique.tsv


cat CTD_TFEA_genelist_not_unique.tsv | awktt '
{
  a[$1]++
}
END {
  for (i in a) print i,a[i]
}'>ctd_genelist_count.tsv

cat ctd_genelist_count.tsv | awktt '
BEGIN {
  while (getline < "ctd_refseq_lincs_gene.tsv") lincs[$1]=$3
}{
  print $0,lincs[$1]
}'| sort -k2 -rn -t$'\t' | awktt '{print NR,$0}' >ctd_count_lincs.tsv

mkdir /mnt/c/Users/Zou/"Google Drive"/okilab/ChIP-Atlas/200819/
cd /mnt/c/Users/Zou/"Google Drive"/okilab/ChIP-Atlas/200819/
sget
get /home/zou/ChIP-Atlas/200820/ctd_count_lincs.tsv


library(ggplot2)

ggplot(data[1:100,], aes(X1, X3, fill=factor(X4))) + geom_bar(stat='identity', width=1) + scale_fill_manual(breaks=levels(data$X4), values=c("#E7E6E6", "#FF6600")) + coord_flip() + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill = "transparent",color = NA),
    panel.grid.minor = element_line(color = NA), 
    panel.grid.major = element_line(color = NA),
    plot.background = element_rect(fill = "transparent",color = NA),
    axis.line = element_line(colour="black",size=1)
    ) + scale_x_continuous(trans = "reverse") + scale_y_continuous(position = "right")
    
    
    
cat ~/ChIP-Atlas/200616/200616_creeds_disease_gene.tsv | awku 3 >creeds_genelist_unique.tsv
cat ~/ChIP-Atlas/200616/200616_creeds_disease_gene.tsv | awktt '{print $3}' >creeds_genelist_not_unique.tsv

cat creeds_genelist_not_unique.tsv | awktt '
{
  a[$1]++
}
END {
  for (i in a) print i,a[i]
}'>creeds_genelist_count.tsv

cat creeds_genelist_unique.tsv | awktt '
BEGIN {
  while ("cat ~/ChIP-Atlas/200707/refseqgene.tsv" | getline) ref[$2]++
  while (getline < "lincs_genelist.tsv") lincs[$0]++
}{
  if (ref[$1] > 0 && lincs[$1] > 0) print $0,"1","1"
  else if (ref[$1] > 0) print $0,"1","0"
  else if (lincs[$1] > 0) print $0,"0","1"
  else print $0,"0","0"
}' >creeds_refseq_lincs_gene.tsv

cat creeds_genelist_count.tsv | awktt '
BEGIN {
  while (getline < "creeds_refseq_lincs_gene.tsv") lincs[$1]=$3
}{
  print $0,lincs[$1]
}'| sort -k2 -rn -t$'\t' | awktt '{print NR,$0}' >creeds_count_lincs.tsv


mkdir /mnt/c/Users/Zou/"Google Drive"/okilab/ChIP-Atlas/200819/
cd /mnt/c/Users/Zou/"Google Drive"/okilab/ChIP-Atlas/200819/
sget
get /home/zou/ChIP-Atlas/200820/ctd_count_lincs.tsv
get /home/zou/ChIP-Atlas/200820/creeds_count_lincs.tsv

library(ggplot2)

ggplot(data, aes(X1, X3, fill=factor(X4))) + geom_bar(stat='identity', width=1) + scale_fill_manual(breaks=levels(data$X4), values=c("#E7E6E6")) + coord_flip() + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill = "transparent",color = NA),
    panel.grid.minor = element_line(color = NA), 
    panel.grid.major = element_line(color = NA),
    plot.background = element_rect(fill = "transparent",color = NA),
    axis.line = element_line(colour="black",size=1)
    ) + scale_x_continuous(trans = "reverse") + scale_y_continuous(position = "right")
    
    
library(ggplot2)

ggplot(data, aes(X1, X4, fill=factor(X3))) + geom_bar(stat='identity', width=1) + geom_abline(intercept=0, slope=-1/19717, lwd=2) + scale_fill_manual(breaks=levels(data$X3),values=c("#32CD32")) + coord_flip() + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill = "transparent",color = NA),
    panel.grid.minor = element_line(color = NA), 
    panel.grid.major = element_line(color = NA),
    plot.background = element_rect(fill = "transparent",color = NA),
    axis.line = element_line(colour="black",size=1)
    ) + scale_x_continuous(trans = "reverse") + scale_y_continuous(position = "right") 

get /home/zou/ChIP-Atlas/200820/CTD_TFEA_genelist_all.tsv
get /home/zou/ChIP-Atlas/200820/creeds_genelist_unique.tsv
get /home/zou/ChIP-Atlas/200820/lincs_genelist.tsv

cat CTD_TFEA_genelist_all.tsv >> "temp"
cat creeds_genelist_unique.tsv >> "temp"
cat lincs_genelist.tsv >> "temp"
cat temp |awku 1 >all_genelist.tsv

cat all_genelist.tsv | awktt '
BEGIN {
  while (getline < "CTD_TFEA_genelist_all.tsv") ctd[$1]++
  while (getline < "creeds_genelist_unique.tsv") creeds[$1]++
  while (getline < "lincs_genelist.tsv") lincs[$1]++
} {
  if ($1 in ctd && $1 in creeds && $1 in lincs) print $1,"1","1","1"
  else if ($1 in ctd && $1 in creeds) print $1,"1","1","0"
  else if ($1 in ctd && $1 in lincs) print $1,"1","0","1"
  else if ($1 in creeds && $1 in lincs) print $1,"0","1","1"
  else if ($1 in ctd) print $1,"1","0","0"
  else if ($1 in creeds) print $1,"0","1","0"
  else if ($1 in lincs) print $1,"0","0","1"
  else print $1,"error"
}' >all_gene_ctd_creeds_lincs.tsv
  

get /home/zou/ChIP-Atlas/200820/all_gene_ctd_creeds_lincs.tsv

# Venn_diagram
library(nVennR)
data <- read.csv("all_gene_ctd_creeds_lincs.tsv",sep="\t")
ctd <- subset(data, CTD == "1")$GENE
creeds <- subset(data, CREEDS == "1")$GENE
lincs <- subset(data, LINCS == "1")$GENE
p <- plotVenn(list(CTD=ctd, CREEDA=creeds, LINCS=lincs))
showSVG(p, outFile = "test.svg", labelRegions = F)



data <- read.csv("C:/Users/Zou/Google Drive/okilab/ChIP-Atlas/200804_論文作成/ctd_gene_percentile.csv",header=F)
g1 <- subset(data)$V7
g2 <- subset(data)$V8
wilcox.test(g1, g2)

# W = 713195, p-value = 2.5e-176 (CTD)

data <- read.csv("C:/Users/Zou/Google Drive/okilab/ChIP-Atlas/200804_論文作成/creeds_gene_percentile.csv",header=F)
g1 <- subset(data)$V1
g2 <- subset(data)$V5
g2[!is.na(g2)]
wilcox.test(g1[!is.na(g1)], g2[!is.na(g2)],correct=FALSE)
U.test(g1[!is.na(g1)], g2[!is.na(g2)])
# W = 11888094, p-value = 3.2e-75 (CREEDS)


