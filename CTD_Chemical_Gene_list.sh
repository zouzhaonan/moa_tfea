#!/bin/sh
#$ -S /bin/sh

DATE=$1

mkdir ~/ChIP-Atlas/CTD/$DATE
cd ~/ChIP-Atlas/CTD/$DATE
# CTD_chem_gene_ixns.tsv
curl http://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz | gunzip >CTD_chem_gene_ixns.tsv
# ChemicalName  ChemicalID  CasRN  GeneSymbol  GeneID  GeneForms  Organism  OrganismID  Interaction  InteractionActions  PubMedIDs

# increased/decreased expression && 欠損値がない record を抽出#
# ChemicalName, ChemicalID, PubMedID, OrganismID, GeneSymbol, Up_Down #
cat CTD_chem_gene_ixns.tsv | grep -v ^# | awk -F "	" '{
 ChemicalName = $1
 ChemicalID = $2
 PubMedID = $11
 OrganismID = $8
 GeneSymbol = $4
 if ($9 == $1" results in increased expression of "$4" mRNA") Up_Down = "Up"
 else if ($9 == $1" results in decreased expression of "$4" mRNA") Up_Down = "Down"
 else if ($9 == $1" results in increased expression of "$4" protein") Up_Down = "Up"
 else if ($9 == $1" results in decreased expression of "$4" protein") Up_Down = "Down"
 else if ($9 == $1" analog results in increased expression of "$4" mRNA") Up_Down = "Up"
 else if ($9 == $1" analog results in decreased expression of "$4" mRNA") Up_Down = "Down"
 else Up_Down = NULL
 printf "%s	%s	%s	%s	%s	%s
", ChemicalName, ChemicalID, PubMedID, OrganismID, GeneSymbol, Up_Down
}' | awk -F "	" '$1 && $2 && $3 && $4 && $5 && $6' | sort | uniq >chem_gene_change_temp1.tsv

# PubMedID "|" 処理 #
cat chem_gene_change_temp1.tsv | awk -F "	" '{
 count = split ($3, PubMedID, /[|]/)
 for (i = 1; i <= count; i++) print $1"	"$2"	"PubMedID[i]"	"$4"	"$5"	"$6
}' | sort | uniq >chem_gene_change.tsv

mkdir ~/ChIP-Atlas/CTD/$DATE/Up
mkdir ~/ChIP-Atlas/CTD/$DATE/Down

# Up or Down #
cd ~/ChIP-Atlas/CTD/$DATE
cat chem_gene_change.tsv | awk -F "	" '{
 if ($6 == "Up") print $0 >>"Up/chem_gene_Up.tsv"
 else print $0 >>"Down/chem_gene_Down.tsv"
}'

# Organism ごとに #
for CHANGE in Up Down
do
 # フォルダ #
 for ORGANISM in hg19 mm9 other
 do
   mkdir ~/ChIP-Atlas/CTD/$DATE/$CHANGE/$ORGANISM
 done
 # genelist #
 cd ~/ChIP-Atlas/CTD/$DATE/$CHANGE/
 cat chem_gene_$CHANGE.tsv | awk -F "	" '{
   ChemicalID = $2
   PubMedID = $3
   OrganismID = $4
   GeneSymbol = $5
   filename = ChemicalID"-"PubMedID"-"OrganismID".txt"
   if (OrganismID == "9606") print GeneSymbol >>"hg19/"filename
   else if (OrganismID == "10090") print GeneSymbol >>"mm9/"filename
   else print GeneSymbol >>"other/"filename
 }'
done
