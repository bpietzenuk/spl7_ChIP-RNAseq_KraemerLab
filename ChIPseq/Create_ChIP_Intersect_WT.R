library(GenomicRanges)
###ChIP
ChIP_plus <- read.delim("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/plus.merged.bed.associated.annotated.edit.txt",header = T )
ChIP_plus_sub <- ChIP_plus[c(1:8)]
head(ChIP_plus_sub)
nrow(ChIP_plus_sub)

ChIP_minus <- read.delim("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/minus.merged.bed.associated.annotated.edited.txt",header = T )
ChIP_minus_sub <- ChIP_minus[c(1:8)]
nrow(ChIP_minus_sub)
head(ChIP_minus_sub)

##WT
WT_Cudef <- read.table("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/M_vs.O.merged.bed.associated.annotated.edited.txt" )
WT_Cudef_sub <- WT_Cudef[c(1:3,10:11)]
WT_Cudef_sub$V1 <- sub("^", "Chr", WT_Cudef_sub$V1 )
colnames(WT_Cudef_sub) <- c("chromosome", "start", "end", "peak_name", "gene.name")
WT_Cudef_sub$width <- (WT_Cudef_sub$end - WT_Cudef_sub$start)+1
head(WT_Cudef_sub)
nrow(WT_Cudef_sub)

WT_Cusuf <- read.table("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/N_vs.P.merged.bed.associated.annotated.edited.txt")
WT_Cusuf_sub <- WT_Cusuf[c(1:4,10:11)]
WT_Cusuf_sub$V1 <- sub("^", "Chr", WT_Cusuf_sub$V1 )
colnames(WT_Cusuf_sub) <- c("chromosome", "start", "end", "peak_name", "gene.name")
WT_Cusuf_sub$width <- (WT_Cusuf_sub$end - WT_Cusuf_sub$start)+1
head(WT_Cusuf_sub)
nrow(WT_Cusuf_sub)

######30% overlap with range from center
###calculating the center of that peak and adding 1/30 of peak width to each side
##WT_suf & ChIPplus
ChIP_minus_sub$peak_center <- (ChIP_minus_sub$width/2) + ChIP_minus_sub$start
ChIP_minus_sub$width_thirty <- ChIP_minus_sub$width * (0.3)
ChIP_minus_sub$start_width_thirty <- ChIP_minus_sub$peak_center - ChIP_minus_sub$width_thirty
ChIP_minus_sub$end_width_thirty <- ChIP_minus_sub$peak_center + ChIP_minus_sub$width_thirty
view(ChIP_minus_sub)

##MINUS
WT_Cudef_sub$peak_center <- (WT_Cudef_sub$width/2) + WT_Cudef_sub$start
WT_Cudef_sub$width_thirty <- WT_Cudef_sub$width * (0.3)
WT_Cudef_sub$start_width_thirty <- WT_Cudef_sub$peak_center - WT_Cudef_sub$width_thirty
WT_Cudef_sub$end_width_thirty <- WT_Cudef_sub$peak_center + WT_Cudef_sub$width_thirty
view(WT_Cudef_sub)

##Merge within the range of start.width and end.width
WT_Cudef_sub_range_thirty <- with(WT_Cusuf_sub, GRanges(WT_Cusuf_sub$gene.name, IRanges(start_width_thirty, end_width_thirty, names=gene.name), "+"))
#minus_start <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(start_width, width=1, names=gene.name), "+"))
#minus_end <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(end_width, width=1, names=gene.name), "+"))
minus_range_thirty <-with(ChIP_minus_sub, GRanges(ChIP_minus_sub$gene.name, IRanges(start_width_thirty, end_width_thirty, names=gene.name), "+"))

#test <- mergeByOverlaps(minus_peak, plus_range)
#str(test)

#test_df <- data.frame(test)

olaps_peak_thirty <- findOverlaps(WT_Cudef_sub_range_thirty, minus_range_thirty)
#olaps_start <- findOverlaps(minus_start, plus_range)
#olaps_end <- findOverlaps(minus_end, plus_range)

queryHits(olaps_peak_thirty)
#queryHits(olaps_start)
#queryHits(olaps_end)

CuDef_intersect_peak_thirty <- cbind(WT_Cudef_sub[queryHits(olaps_peak_thirty),], ChIP_minus_sub[subjectHits(olaps_peak_thirty),])
nrow(CuSuf_intersect_peak_thirty) ### =0
view(CuSuf_intersect_peak)

#####WT CuDef
ChIP_plus_sub$peak_center <- (ChIP_plus_sub$width/2) + ChIP_plus_sub$start
ChIP_plus_sub$width_thirty <- ChIP_plus_sub$width * (0.3)
ChIP_plus_sub$start_width_thirty <- ChIP_plus_sub$peak_center - ChIP_plus_sub$width_thirty
ChIP_plus_sub$end_width_thirty <- ChIP_plus_sub$peak_center + ChIP_plus_sub$width_thirty
view(ChIP_plus_sub)

##MINUS
WT_Cusuf_sub$peak_center <- (WT_Cusuf_sub$width/2) + WT_Cusuf_sub$start
WT_Cusuf_sub$width_thirty <- WT_Cusuf_sub$width * (0.3)
WT_Cusuf_sub$start_width_thirty <- WT_Cusuf_sub$peak_center - WT_Cusuf_sub$width_thirty
WT_Cusuf_sub$end_width_thirty <- WT_Cusuf_sub$peak_center + WT_Cusuf_sub$width_thirty
view(WT_Cusuf_sub)

##Merge within the range of start.width and end.width
WT_Cusuf_sub_range_thirty <- with(WT_Cusuf_sub, GRanges(WT_Cusuf_sub$gene.name, IRanges(start_width_thirty, end_width_thirty, names=gene.name), "+"))
#minus_start <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(start_width, width=1, names=gene.name), "+"))
#minus_end <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(end_width, width=1, names=gene.name), "+"))
plus_range_thirty <-with(ChIP_plus_sub, GRanges(ChIP_plus_sub$gene.name, IRanges(start_width_thirty, end_width_thirty, names=gene.name), "+"))

#test <- mergeByOverlaps(minus_peak, plus_range)
#str(test)

#test_df <- data.frame(test)

olaps_peak_thirty <- findOverlaps(WT_Cusuf_sub_range_thirty, plus_range_thirty)
#olaps_start <- findOverlaps(minus_start, plus_range)
#olaps_end <- findOverlaps(minus_end, plus_range)

queryHits(olaps_peak_thirty)
#queryHits(olaps_start)
#queryHits(olaps_end)

CuSuf_intersect_peak_thirty <- cbind(WT_Cusuf_sub[queryHits(olaps_peak_thirty),], ChIP_plus_sub[subjectHits(olaps_peak_thirty),])
nrow(CuSuf_intersect_peak_thirty)  ### =0
view(CuSuf_intersect_peak)