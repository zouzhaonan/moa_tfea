# Summary 作成 #

DATE=210515
# ChemicalName ChemicalID 辞書 #
cd ~/ChIP-Atlas/CTD/$DATE/
cat CTD_chem_gene_ixns.tsv | grep -v ^# | awk -F "	" '{ print $1"	"$2 }' | sort | uniq >ChemicalName_ID_dictionary.tsv

for ORGANISM in hg19 mm9
do
  # ChIP-Atlas は microRNA に現時点で対応していないので，出力エラーを取り除く． #
  # 成功したもののリスト #
  cd ~/ChIP-Atlas/CTD/$DATE/Results/$ORGANISM
  rm *.bed*
  rm *.tmpForinsilicoChIP*
  for filename in `ls | grep "\.tsv$"`; do
    wc -l $filename 
  done > ../$ORGANISM"_wc_temp.txt"
  cd ~/ChIP-Atlas/CTD/$DATE/Results
  cat $ORGANISM"_wc_temp.txt" | awk -F " " ' $1 > 0 {gsub (/[ \.]/, "	"); print $2}' >$ORGANISM"_TEFA_done_list_temp.txt"
  
  # 失敗したもののリスト #
  cat ../"filelist_"$ORGANISM"_common_above_ten.txt" | awk -F "	" -v ORGANISM="$ORGANISM" '
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
    cat $filename.tsv | awk -F "	" -v ORGANISM="$ORGANISM" -v filename="$filename" '{
      if ($11 > 1) {
        True="../"ORGANISM"_summary_TRUE_temp.tsv"
        print filename"	"$0 >>True
        exit
      }
    }'
    cat $filename.tsv | awk -F "	" -v ORGANISM="$ORGANISM" -v filename="$filename" '{
      if ($11 <= 1) {
        False="../"ORGANISM"_summary_FALSE_temp.tsv"
        print filename"	"$0 >>False
        exit
      }
    }'
  done
    
  # ChemicalName を追加 #
  cd ~/ChIP-Atlas/CTD/$DATE/Results/
  for FE in TRUE FALSE
  do
    cat $ORGANISM"_summary_"$FE"_temp.tsv" | awk -F "	" '
    BEGIN{
      OFS="	"
    }{
      gsub ("-", "	", $1)
      print $0
    }' | awk -F "	" '
    BEGIN {
      OFS="	"
      while ((getline < "../ChemicalName_ID_dictionary.tsv") > 0)
      Name_of[$2] = $1
    }{
      print Name_of[$1]"	"$0
    }' > $ORGANISM"_summary_"$FE".tsv"
  done
done
