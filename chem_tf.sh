#!/bin/sh
#$ -S /bin/sh
label=$1
cd ~/ChIP-Atlas/200805/
for chem in `awku 2 "chem_temp_"$label".txt"`
do
  cat ~/ChIP-Atlas/200805/chemID_TFno_ES.tsv | awk -F "	" -v OFS="	" -v chem="$chem" '$1==chem' >$chem"_temp.tsv"
  for ((i=1;i<=997;i++))
  do
    cat $chem"_temp.tsv" | awk -F "	" -v OFS="	" -v i="$i" -v chem="$chem" '$1==chem {
      if ($2==i) printf "%s	",$3
    }'
    a=`cat $chem"_temp.tsv" | awk -F "	" -v OFS="	" -v i="$i" -v chem="$chem" '$1==chem {
      if ($2==i) printf "%s	",$3
    }'`
    if [[ $a == "" ]]; then
      echo -n -e "0	"
    fi
    if [[ $i == 997 ]]; then
      echo $chem
    fi
  done
done >"chem_tf_matrix_temp_"$label".tsv"
