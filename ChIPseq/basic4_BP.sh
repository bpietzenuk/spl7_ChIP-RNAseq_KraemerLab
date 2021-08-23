#!/bin/bash
export LC_ALL=C

###exclude Mito & Chloro
ls $1"_peaks.narrowPeak" $2"_peaks.narrowPeak" $3"_peaks.narrowPeak" $4"_peaks.narrowPeak"
for i in $1 $2 $3 $4
do 
   grep -v mito $i"_peaks.narrowPeak" | grep -v chloro > $i"_peaks.narrowPeakq"
done 

exit
#exit
#make the initial overlap 
python overlapSimple.py $1"_peaks.narrowPeakq" $2"_peaks.narrowPeakq" $3"_peaks.narrowPeakq" $4"_peaks.narrowPeakq" > $5".soverlap" ####output from MACs
#exit
#exit

#determine overlap sizes 
python ovSizes.py < $5".soverlap" > $5".ovsizes"

#exit
#exit
#make overlap sizes plot
Rscript --vanilla ov.size.plot.R $5".ovsizes" $5".png" $5
#exit
#filtered overlap on 80% overlap
python d80.py < $5".soverlap" > $5".fsoverlap"  ###Creates cluster out of the 80% overlap rule
#exit
python allClusterized.py filt $1"_peak" $2"_peak" $3"_peak" $4"_peak" < $5".fsoverlap" ###goes for Output of d80.py; checks all replicates and give peaks common name sorts back to replicates
#
#exit
cat $1"_peak.filt.output" $2"_peak.filt.output" $3"_peak.filt.output" $4"_peak.filt.output" | sort | uniq -c | awk '$1 >= 2' | awk '{ print $2 }' | grep -v -i mitoq | grep -v -i chloroq > $5".2timesPeaks" 
cat $1"_peak.filt.output" $2"_peak.filt.output" $3"_peak.filt.output" $4"_peak.filt.output" | sort | uniq -c | awk '$1 >= 2' | awk '{ print $2"\t"$1 }' | python modo.py > $5."peak_occurances.txt"
#exit
#exit

## gives number of replicates where HIT comes appears
python rebuildFromMerged.py <  $5".2timesPeaks" | grep -v -i mito | grep -v -i cql >  $5".merged.bed"
#exit

###
python singlePointPeaks.py < $5".merged.bed" > $5".merged.bedc"
awk '{ print $4"\t"$1"\t"$2"\t"$3"\t+" }' < $5".merged.bed" | grep -v -i mito | grep -v -i cql > $5".peaks.for.extraction"
#exit 

###Get sequence from peaks from genome
python extractRegionsFromFastaC.py Athaliana.fasta < $5".peaks.for.extraction" > $5".extractedPicos.txt"

#exit
####check YOUR R-Version and ADD folder!
Rscript --vanilla tvenn4.R $1"_peak.filt.output" $2"_peak.filt.output" $3"_peak.filt.output" $4"_peak.filt.output" $5".venn.png" $1 $2 $3 $4
#exit
samtools faidx $5".extractedPicos.txt" 
#exit

###Connection between peaks and genes - 
Rscript --vanilla doAnnotationsBaseCircular.R tair.gff $5".merged.bed"
cut -f 2,3,4,5,8,9 $5.merged.bed.associated | sed 's/Peak./Peak-/' > $5.merged.bed.associated.modified
#exit
#exit
###Maps annotation file with peaks and associated peaks with genes
cat $5.merged.bed.associated.modified | python col_mapper.py $5."peak_occurances.txt" -2 | python col_mapper.py newMappable.txt -2 > $5.merged.bed.associated.annotated.txt

###make list of peaks correspoding to feature
grep -i trna $5.merged.bed.associated.annotated.txt | cut -f 5 > trna.peaks.to.exclude
grep -i -v trna $5.merged.bed.associated.annotated.txt | cut -f 5 > promotor.without.trna.to.include

###filter by created feature lists
python excludeSequences.py trna.peaks.to.exclude <  $5".extractedPicos.txt" >  $5".extractedPicos.txt.M" 
python includeSequences.py promotor.without.trna.to.include <  $5".extractedPicos.txt" >  $5".extractedPicos.txt.MP"

###runs meme
meme-chip -neg ../../minus.extractedPicos.woIntersect.txt -db ../../ArabidopsisDAPv1.meme  -meme-maxw 10 -meme-minw 4 $i -oc results/$i"_Discriminative_woIntersect" -meme-mod anr -meme-nmotif 20 > $i.woIntersect.log.txt