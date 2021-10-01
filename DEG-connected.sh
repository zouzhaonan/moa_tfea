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
  gsub ("{"do_id"","
{"do_id"")
  gsub ("\"down_genes\"","
\"down_genes\"")
  gsub ("\"up_genes\"","
\"up_genes\"")
  gsub ("\"disease_name\"","
\"disease_name\"")
  gsub ("\"curator\"","
\"curator\"")
  print $0
}' | awktt 'NR>1'| awktt '{
  a[NR]=$0
  if (NR%5==4) printf "%s	%s
%s	%s
", a[NR], a[NR-2], a[NR], a[NR-1]
}' | awktt '{
  gsub ("\"","	",$1)
  gsub ("down_genes","down",$2)
  gsub ("up_genes","up",$2)
  gsub (": ","	",$2)
  gsub ("\"","",$2)
  gsub ("\]\]","]",$2)
  gsub ("\[","",$2)
  print $0
}' | awktt '$8 == "human" { print $12,$8,$14,$15 }'| awktt '{
  count=split($4,x,"], ")
  for (i=1;i<count;i++) print $1,$2,$3,x[i]
}' | awktt '{ gsub(",","	",$4); print $0 }' | awktt '{print $1,$4,$3}' | awktt '
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
  cat chem_gene_original.tsv | awk -F"	" -v OFS="	" -v i="$i" '
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
  cat diz_gene_original.tsv | awk -F"	" -v OFS="	" -v i="$i" '
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
  cat chem_gene_original.tsv | awk -F"	" -v OFS="	" -v i="$i" '$1==i' >chemgene/$i".tsv"
done

for i in $(cat dizlist_10.tsv)
do
  cat diz_gene_original.tsv | awk -F"	" -v OFS="	" -v i="$i" '$1==i' >dizgene/$i".tsv"
done

curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz | gunzip> refFlat.txt
awku 1 refFlat.txt |sort |awktt '{print NR,$1}'>refseqgene.tsv

# chem/diz gene profile
# Up/Down/others

mkdir chemgeneprof
mkdir dizgeneprof

for i in $(cat chemlist_10.tsv)
do
  cat refseqgene.tsv | awk -F"	" -v OFS="	" -v i="$i" '
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
  cat refseqgene.tsv | awk -F"	" -v OFS="	" -v i="$i" '
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
  data0=read.table(file.name, sep="	")
  data=data0[1:2,1:2]
  fisher <- fisher.test(data)
  str <- paste(file.name, fisher$p.value, sep="	")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200707/up_vs_down_result_temp.tsv", append=TRUE, quote=F, col.names=F,)
}
=======================================

# fisher.test, deg vs non-deg #
=======================================
R
files <- list.files(pattern="tsv")
for (file.name in files) {
  data0=read.table(file.name, sep="	")
  a=data0[1,1]+data0[1,2]+data0[2,1]+data0[2,2]
  b=data0[1,3]+data0[2,3]
  c=data0[3,1]+data0[3,2]
  d=data0[3,3]
  data=matrix(c(a, b, c, d), nrow = 2, ncol = 2)  
  fisher <- fisher.test(data)
  str <- paste(file.name, fisher$p.value, sep="	")
  write.table (str, "/lustre7/home/zou/ChIP-Atlas/200707/degvsnondeg_fisher_result_temp.tsv", append=TRUE, quote=F, col.names=F,)
}
=======================================

cd ~/ChIP-Atlas/200707/
cat degvsnondeg_fisher_result_temp.tsv | awktt '{gsub("1 ","",$1);print $1,-log($2)/log(10)}'|sort -n -r -k2 > degvsnondeg_fisher_result.tsv
cat up_vs_down_result_temp.tsv | awktt '{gsub("1 ","",$1);print $1,-log($2)/log(10)}'|sort -n -r -k2 > updown_fisher_result.tsv

answer="200615_curated_chemical_disease.tsv"


for i in updown_fisher
do
  cat $i"_result.tsv" | awktt '{gsub("_x_","	");gsub("_","	");gsub(".tsv","");print $0}'| awktt '
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
rocdata <- read.table (file.name, sep="	")
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
  cat $i"_result.tsv" | awktt '{gsub("_x_","	");gsub("_","	");gsub(".tsv","");print $0}'| awktt '
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
rocdata <- read.table (file.name, sep="	")
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
