for s in $@
do
#    mkdir $s
    cp $s.bed.sepa minus.bed.sepa
done

for i in `ls | grep samples_minus`
do 
  u=`echo $i | cut -f 2 -d "_"`
  echo $i $u
  Rscript --vanilla plotDistancesATG.R minus.bed.sepa $i minus.bed.sepa.csv
  Rscript --vanilla plotDistanceGGplotStyle.R minus.bed.sepa.csv $s.rel.png $s.rel.pdf
done 
