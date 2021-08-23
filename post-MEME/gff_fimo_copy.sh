#!/bin/bash
export LC_ALL=C

for i in $@
do
#	mkdir "fimo_gffs_out/centered/"$i/
	k=`ls -d ara/centered/$i/fimo_out_*/ | wc -l`
#	echo $k
	for ((j = 1; j <= k; j++))
	do
	nl=`<ara/centered/$i/fimo_out_$j/fimo.gff wc -l`
#	echo $nl
	if [ $nl \> "1" ]
	then
#		echo $nl
#		echo $j
		n=`head -2 "ara/centered/"$i"/fimo_out_"$j"/fimo.gff" | tail -1 | cut -f9 | cut -f 7 -d "="`
#		echo $n
#		mkdir "fimo_gffs_out/centered/"$i/
		cp "ara/centered"/$i/"fimo_out_"$j/"fimo.gff" "fimo_gffs_out/centered/"$i.$n".gff"
	fi
	done
done