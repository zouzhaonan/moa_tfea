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
#   gsub ("\|","	")
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
}'| sort -k2 -rn -t$'	' | awktt '{print NR,$0}' >ctd_count_lincs.tsv

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
}'| sort -k2 -rn -t$'	' | awktt '{print NR,$0}' >creeds_count_lincs.tsv


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
data <- read.csv("all_gene_ctd_creeds_lincs.tsv",sep="	")
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
