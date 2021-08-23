for i in $@
do
	python buildDataQA.py tair.gff  3000 1000 < $i | python adjustValuesA.py > $i.sepa
done