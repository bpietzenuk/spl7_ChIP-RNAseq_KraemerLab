library( ChIPpeakAnno )
library( GenomicFeatures )
args = commandArgs(trailingOnly=TRUE)

txdb = makeTxDbFromGFF( args[1], format = "gff3" )

krng = toGRanges( txdb, feature = "gene" ) 
#q()
baseFolder = getwd()
dataSets = c( args[2] ) 

for ( dta in dataSets ) {
peaks = toGRanges( paste( baseFolder, "/", dta, sep = "" ), format = "BED" )

ov = annoPeaks( peaks, krng, bindingType = "fullRange", bindingRegion = c( -3000,1000) )

write.table( ov, file = paste( baseFolder, "/", dta, ".associated", sep = "" ), quote = F, col.names = T, row.names = T, sep = "\t" )

P = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns") 
z = assignChromosomeRegion( peaks, proximal.promoter.cutoff=3000L, immediate.downstream.cutoff=1000L, TxDb = txdb,
precedence=c("Promoters", "immediateDownstream",                                                                                                  "fiveUTRs", "threeUTRs", 
"Exons", "Introns") ) 
 barplot( z$percentage )
 
 
#binOverFeature( peaks, annotationData = krng, radius = 5000, nbins = 20, FUN = length, errFun = 0,
#                ylab = "count", "Around TSS" )
print( z ) 
 write.table( z$percentage, file = paste( baseFolder, "/", dta, ".features.txt", sep = "" ), row.names = T, col.names = T, quote = F, sep = "\t" )

}
