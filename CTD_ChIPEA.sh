#############
### TFEA ####
#############

# Chemical_gene_interaction_pattern から #
# AAA results in increased/decreased expression of BBB mRNA #
# AAA results in increased/decreased expression of BBB protein #
# AAA analog results in increased/decreased expression of BBB mRNA #
# という記述のみを抽出 #

# STEP. 1/7
# ========================
DATE="200515"

# スパコン 
qsub -l s_vmem=16G -l mem_req=16G -e /dev/null -o /dev/null ~/ChIP-Atlas/CTD/Chemical_Gene_list.sh $DATE


# insilicoChIP に使える genelist を準備する #
# 1. gene 数 10 以上 #
# 2. Up/Down 両方のデータある #

# STEP. 2/7
# ========================
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


# STEP. 3/7
# =====================
for ORGANISM in hg19 mm9
do
  # Up/Down 両方のデータがある record を抽出 #
  cd ~/ChIP-Atlas/CTD/$DATE/
  cat "filelist_"$ORGANISM"_Up_temp_above_ten.txt" | awk -F "	" -v ORGANISM="$ORGANISM" '    # 外部変数の渡し方に注意！！ #
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
#======================
