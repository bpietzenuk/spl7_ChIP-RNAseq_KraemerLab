library(GenomicRanges)
ChIP_plus <- read.delim("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/plus.merged.bed.associated.annotated.edit.txt",header = T )
ChIP_plus_sub <- ChIP_plus[c(1:8)]
head(ChIP_plus_sub)
nrow(ChIP_plus_sub)

ChIP_minus <- read.delim("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/minus.merged.bed.associated.annotated.edited.txt",header = T )
ChIP_minus_sub <- ChIP_minus[c(1:8)]
nrow(ChIP_minus_sub)
head(ChIP_minus_sub)


######30% width range
####Intersect by center position of ChIP peak
###calculating the center of that peak and adding 1/10 of peak width to each side
##PLUS
ChIP_plus_sub$peak_center <- (ChIP_plus_sub$width/2) + ChIP_plus_sub$start
ChIP_plus_sub$width_thirty <- ChIP_plus_sub$width * (0.3)
ChIP_plus_sub$start_width_thirty <- ChIP_plus_sub$peak_center - ChIP_plus_sub$width_thirty
ChIP_plus_sub$end_width_thirty <- ChIP_plus_sub$peak_center + ChIP_plus_sub$width_thirty
view(ChIP_plus_sub)

##MINUS
ChIP_minus_sub$peak_center <- (ChIP_minus_sub$width/2) + ChIP_minus_sub$start
ChIP_minus_sub$width_thirty <- ChIP_minus_sub$width * (0.3)
ChIP_minus_sub$start_width_thirty <- ChIP_minus_sub$peak_center - ChIP_minus_sub$width_thirty
ChIP_minus_sub$end_width_thirty <- ChIP_minus_sub$peak_center + ChIP_minus_sub$width_thirty
view(ChIP_minus_sub)

##Merge within the range of start.width and end.width
minus_range_thirty <- with(ChIP_minus_sub, GRanges(ChIP_minus_sub$gene.name, IRanges(start_width_thirty, end_width_thirty, names=gene.name), "+"))
#minus_start <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(start_width, width=1, names=gene.name), "+"))
#minus_end <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(end_width, width=1, names=gene.name), "+"))
plus_range_thirty <-with(ChIP_plus_sub, GRanges(ChIP_plus_sub$gene.name, IRanges(start_width_thirty, end_width_thirty, names=gene.name), "+"))

#test <- mergeByOverlaps(minus_peak, plus_range)
#str(test)

#test_df <- data.frame(test)

olaps_peak_thirty <- findOverlaps(minus_range_thirty, plus_range_thirty)
#olaps_start <- findOverlaps(minus_start, plus_range)
#olaps_end <- findOverlaps(minus_end, plus_range)

queryHits(olaps_peak_thirty)
#queryHits(olaps_start)
#queryHits(olaps_end)

ChIP_intersect_peak_thirty <- cbind(ChIP_minus_sub[queryHits(olaps_peak_thirty),], ChIP_plus_sub[subjectHits(olaps_peak_thirty),])
nrow(ChIP_intersect_peak_thirty)
view(ChIP_intersect_peak)
##Validation test
anno <- read.delim("/Volumes/JARVIS-X/PostDoc/Plant_Genomic_Ressources/00 Functional Description for paper/TAIR10_functional_descriptions_for_Meetings_Mapman_strict_genes_202008_03_bHLHadded.txt", header = T)
names(anno)[1] <- "AGI"
head(anno)

names(ChIP_intersect_peak_thirty)[6] <- "AGI"
test_val_thirty <- merge(ChIP_intersect_peak_thirty, anno, by.x = "AGI", by.y = "AGI", all.x=T)
write.table(ChIP_intersect_peak_thirty, "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/l2FC05/ChIP_intersect_peak(final).txt", sep = "\t", row.names = F)

####
annotation <- read.delim("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/novel.tair.annot", header=F)
intersect <- read.table("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/ChIP_intersect_peak(final).txt", header = T)
#intersect_dedup <- intersect[!duplicated(intersect[ , c("peak","peak.1")]),]
intersect_anno <- merge(intersect, anno, by.x="gene.name", by.y = "AGI", all.x=T)
write.table(intersect_anno,"/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/intersect.merged.bed.associated.annotated.edited.txt", row.names = F, sep="\t")
nrow(intersect[intersect$gene.name==intersect$gene.name.1,]) ###check if gene names = equal in each column)
head(intersect)
newintersect <- intersect[c(1,10,22,5,17,6)]
head(newintersect)
#test2 <- newintersect[duplicated(newintersect$gene.name),]
colnames(newintersect) <- c("Chromosome","center_ChIPminus","center_ChIPplus", "peak_ChIPminus", "peak_ChIPplus","associated_gene")
nrow(newintersect)
newintersectdeduplicated_by_peaks <- newintersect[!duplicated(newintersect[ , c("peak_ChIPminus","peak_ChIPplus")]),]
nrow(newintersectdeduplicated_by_peaks)
head(newintersectdeduplicated_by_peaks)
newintersectdeduplicated_by_peaks$center_ChIPminus <- round(newintersectdeduplicated_by_peaks$center_ChIPminus,digits=0)
newintersectdeduplicated_by_peaks$center_ChIPplus <- round(newintersectdeduplicated_by_peaks$center_ChIPplus,digits=0)
newintersectdeduplicated_by_peaks$neu_center_ChIPminus<-ifelse(newintersectdeduplicated_by_peaks$center_ChIPminus>newintersectdeduplicated_by_peaks$center_ChIPplus,newintersectdeduplicated_by_peaks$center_ChIPplus,newintersectdeduplicated_by_peaks$center_ChIPminus)
newintersectdeduplicated_by_peaks$neu_center_ChIPplus<-ifelse(newintersectdeduplicated_by_peaks$center_ChIPminus>newintersectdeduplicated_by_peaks$center_ChIPplus,newintersectdeduplicated_by_peaks$center_ChIPminus,newintersectdeduplicated_by_peaks$center_ChIPplus)
head(newintersectdeduplicated_by_peaks)
newintersect_deduplicatedByPeaks_start_end_corrected <- newintersectdeduplicated_by_peaks[c(1,7,8,4,5,6)]
head(newintersect_deduplicatedByPeaks_start_end_corrected)
write.table(newintersect_deduplicatedByPeaks_start_end_corrected,"/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/newintersectdeduplicated_by_peaks.txt", row.names = F, col.names = F, sep="\t")