#!/bin/sh
#$ -S /bin/sh

DATE=$1 
CHANGE=$2
ORGANISM=$3
# gene 数が 10 未満のやつを $CHANGE/$ORGANISM/genes_under_ten/ に移動 #
mkdir ~/ChIP-Atlas/CTD/$DATE/$CHANGE/$ORGANISM/genes_under_ten/
cd ~/ChIP-Atlas/CTD/$DATE/$CHANGE/$ORGANISM/
for filename in `ls | grep "\.txt$"`
do
 wc -l $filename
done >~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_"$CHANGE"_temp_wc.txt"
cd ~/ChIP-Atlas/CTD/$DATE/
cat "filelist_"$ORGANISM"_"$CHANGE"_temp_wc.txt" | awk -F " " ' $1 < 10 { print $2 }' >"filelist_"$ORGANISM"_"$CHANGE"_temp_under_ten.txt"
for filename in `cat "filelist_"$ORGANISM"_"$CHANGE"_temp_under_ten.txt"`
do
 mv $CHANGE/$ORGANISM/$filename $CHANGE/$ORGANISM/genes_under_ten/
done
# gene 数が 10 以上のやつをリストアップ #
ls $CHANGE/$ORGANISM/ | grep "\.txt$" >~/ChIP-Atlas/CTD/$DATE/"filelist_"$ORGANISM"_"$CHANGE"_temp_above_ten.txt"
