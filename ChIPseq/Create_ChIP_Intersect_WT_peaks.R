ChIP_plus_sub <- read.table("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/plus.merged.bed" )
#ChIP_plus_sub <- ChIP_plus[c(1:8)]
colnames(ChIP_plus_sub) <- c("chromosome", "start", "end", "peak_name_controlCu_ChIP","replicates")
ChIP_plus_sub$width <- (ChIP_plus_sub$end - ChIP_plus_sub$start)+1
median(ChIP_plus_sub$width)
head(ChIP_plus_sub)
nrow(ChIP_plus_sub)

ChIP_minus_sub <- read.table("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/minus.merged.bed")
#ChIP_minus_sub <- ChIP_minus[c(1:8)]
colnames(ChIP_minus_sub) <- c("chromosome", "start", "end", "peak_name_lowCu_ChIP","replicates")
ChIP_minus_sub$width <- (ChIP_minus_sub$end - ChIP_minus_sub$start)+1
median(ChIP_minus_sub$width)
nrow(ChIP_minus_sub)
head(ChIP_minus_sub)

##WT
WT_Cudef <- read.table("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/M_vs.O.merged.bed" )
nrow(WT_Cudef)
count(as.data.frame(table(WT_Cudef$V4)))
sum(duplicated(WT_Cudef$V4))
WT_Cudef_sub <- WT_Cudef[c(1:4)]
WT_Cudef_sub$V1 <- sub("^", "Chr", WT_Cudef_sub$V1 )
colnames(WT_Cudef_sub) <- c("chromosome", "start", "end", "peak_name_lowCuWT")
WT_Cudef_sub$width <- (WT_Cudef_sub$end - WT_Cudef_sub$start)+1
median(WT_Cudef_sub$width)
head(WT_Cudef_sub)
nrow(WT_Cudef_sub)

WT_Cusuf <- read.table("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/N_vs.P.merged.bed")
nrow(WT_Cusuf)
WT_Cusuf_sub <- WT_Cusuf[c(1:4)]
WT_Cusuf_sub$V1 <- sub("^", "Chr", WT_Cusuf_sub$V1 )
colnames(WT_Cusuf_sub) <- c("chromosome", "start", "end", "peak_name_controlCuWT")
WT_Cusuf_sub$width <- (WT_Cusuf_sub$end - WT_Cusuf_sub$start)+1
median(WT_Cusuf_sub$width)
head(WT_Cusuf_sub)
nrow(WT_Cusuf_sub)
######30% width range
####Intersect by center position of ChIP peak
###calculating the center of that peak and adding 1/30 of peak width to each side
##PLUS
ChIP_plus_sub$peak_center <- (ChIP_plus_sub$width/2) + ChIP_plus_sub$start
ChIP_plus_sub$width_thirty <- ChIP_plus_sub$width * (0.3)
ChIP_plus_sub$start_width_thirty <- ChIP_plus_sub$peak_center - ChIP_plus_sub$width_thirty
ChIP_plus_sub$end_width_thirty <- ChIP_plus_sub$peak_center + ChIP_plus_sub$width_thirty
head(ChIP_plus_sub)

##MINUS
WT_Cusuf_sub$peak_center <- (WT_Cusuf_sub$width/2) + WT_Cusuf_sub$start
WT_Cusuf_sub$width_thirty <- WT_Cusuf_sub$width * (0.3)
WT_Cusuf_sub$start_width_thirty <- WT_Cusuf_sub$peak_center - WT_Cusuf_sub$width_thirty
WT_Cusuf_sub$end_width_thirty <- WT_Cusuf_sub$peak_center + WT_Cusuf_sub$width_thirty
head(WT_Cusuf_sub)

##Merge within the range of start.width and end.width
WT_Cusuf_sub_range_thirty <- with(WT_Cusuf_sub, GRanges(WT_Cusuf_sub$chromosome, IRanges(start_width_thirty, end_width_thirty, names=chromosome), "+"))
#minus_start <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(start_width, width=1, names=gene.name), "+"))
#minus_end <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(end_width, width=1, names=gene.name), "+"))
plus_range_thirty <-with(ChIP_plus_sub, GRanges(ChIP_plus_sub$chromosome, IRanges(start_width_thirty, end_width_thirty, names=chromosome), "+"))

#test <- mergeByOverlaps(minus_peak, plus_range)
#str(test)

#test_df <- data.frame(test)

WT_ChIP_Cusuf_olaps_peak_thirty <- findOverlaps(WT_Cusuf_sub_range_thirty, plus_range_thirty)
#olaps_start <- findOverlaps(minus_start, plus_range)
#olaps_end <- findOverlaps(minus_end, plus_range)

queryHits(WT_ChIP_Cusuf_olaps_peak_thirty)
#queryHits(olaps_start)
#queryHits(olaps_end)

WT_Cusuf_intersect_peak_thirty <- cbind(WT_Cusuf_sub[queryHits(WT_ChIP_Cusuf_olaps_peak_thirty),], ChIP_plus_sub[subjectHits(WT_ChIP_Cusuf_olaps_peak_thirty),])
nrow(WT_Cusuf_intersect_peak_thirty)
count(as.data.frame(table(WT_Cusuf_intersect_peak_thirty$peak_name_controlCu_ChIP)))
sum(duplicated(WT_Cusuf_intersect_peak_thirty$peak_name_controlCu_ChIP))
colnames(WT_Cusuf_intersect_peak_thirty) <- c("chromosome_controlCuWT", "start_controlCuWT", "end_controlCuWT", "peak_name_controlCuWT","width_controlCuWT", "peak_center_controlCuWT", "width_30_controlCuWT", "start_width_30_controlCuWT", "end_width_30_controlCuWT", "chromosome_controlCu_ChIP", "start_controlCu_ChIP", "end_controlCu_ChIP" ,"peak_name_controlCu_ChIP", "replicates_controlCu_ChIP", "width_controlCu_ChIP", "peak_center_controlCu_ChIP", "width_30_controlCu_ChIP", "start_width_30_controlCu_ChIP", "end_width_30_controlCu_ChIP")
t1 <- WT_Cusuf_intersect_peak_thirty %>% distinct(peak_name_controlCu_ChIP, .keep_all = TRUE)
count(as.data.frame(table(t1$peak_name_controlCu_ChIP)))
sum(duplicated(t1$peak_name_controlCu_ChIP))
head(t1)
write.table(t1, "/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/WT_intersect/WT_ChIP_Cusuf_intersect_peak_thirty_deduplicated.bed", sep = "\t", row.names = F, col.names = F)

##Minus
ChIP_minus_sub$peak_center <- (ChIP_minus_sub$width/2) + ChIP_minus_sub$start
ChIP_minus_sub$width_thirty <- ChIP_minus_sub$width * (0.3)
ChIP_minus_sub$start_width_thirty <- ChIP_minus_sub$peak_center - ChIP_minus_sub$width_thirty
ChIP_minus_sub$end_width_thirty <- ChIP_minus_sub$peak_center + ChIP_minus_sub$width_thirty
head(ChIP_minus_sub)

##MINUS
WT_Cudef_sub$peak_center <- (WT_Cudef_sub$width/2) + WT_Cudef_sub$start
WT_Cudef_sub$width_thirty <- WT_Cudef_sub$width * (0.3)
WT_Cudef_sub$start_width_thirty <- WT_Cudef_sub$peak_center - WT_Cudef_sub$width_thirty
WT_Cudef_sub$end_width_thirty <- WT_Cudef_sub$peak_center + WT_Cudef_sub$width_thirty
head(WT_Cudef_sub)

##Merge within the range of start.width and end.width
WT_Cudef_sub_range_thirty <- with(WT_Cudef_sub, GRanges(WT_Cudef_sub$chromosome, IRanges(start_width_thirty, end_width_thirty, names=chromosome), "+"))
#minus_start <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(start_width, width=1, names=gene.name), "+"))
#minus_end <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(end_width, width=1, names=gene.name), "+"))
minus_range_thirty <-with(ChIP_minus_sub, GRanges(ChIP_minus_sub$chromosome, IRanges(start_width_thirty, end_width_thirty, names=chromosome), "+"))

#test <- mergeByOverlaps(minus_peak, plus_range)
#str(test)

#test_df <- data.frame(test)

WT_ChIP_Cudef_olaps_peak_thirty <- findOverlaps(WT_Cudef_sub_range_thirty, minus_range_thirty)
#olaps_start <- findOverlaps(minus_start, plus_range)
#olaps_end <- findOverlaps(minus_end, plus_range)

queryHits(WT_ChIP_Cudef_olaps_peak_thirty)
#queryHits(olaps_start)
#queryHits(olaps_end)

WT_ChIP_Cudef_intersect_peak_thirty <- cbind(WT_Cudef_sub[queryHits(WT_ChIP_Cudef_olaps_peak_thirty),], ChIP_minus_sub[subjectHits(WT_ChIP_Cudef_olaps_peak_thirty),])
nrow(WT_ChIP_Cudef_intersect_peak_thirty)
head(WT_ChIP_Cudef_intersect_peak_thirty)
colnames(WT_ChIP_Cudef_intersect_peak_thirty) <- c("chromosome_lowCuWT", "start_lowCuWT", "end_lowCuWT", "peak_name_lowCuWT","width_lowCuWT", "peak_center_lowCuWT", "width_30_lowCuWT", "start_width_30_lowCuWT", "end_width_30_lowCuWT", "chromosome_lowCu_ChIP", "start_lowCu_ChIP", "end_lowCu_ChIP" ,"peak_name_lowCu_ChIP", "replicates_lowCu_ChIP", "width_lowCu_ChIP", "peak_center_lowCu_ChIP", "width_30_lowCu_ChIP", "start_width_30_lowCu_ChIP", "end_width_30_lowCu_ChIP")
t2 <- WT_ChIP_Cudef_intersect_peak_thirty %>% distinct(peak_name_lowCu_ChIP, .keep_all = TRUE)
count(as.data.frame(table(t2$peak_name_lowCu_ChIP)))
sum(duplicated(t2$peak_name_lowCu_ChIP))
head(t2)
write.table(t2, "/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/WT_intersect/WT_ChIP_Cudef_intersect_peak_thirty_dedup.txt", sep = "\t", row.names = F)
 
####Intersect of Intersect using ChIP

t1_range_thirty <- with(t1, GRanges(t1$chromosome_controlCu_ChIP, IRanges(start_width_30_controlCu_ChIP, end_width_30_controlCu_ChIP, names=chromosome_controlCu_ChIP), "+"))
#minus_start <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(start_width, width=1, names=gene.name), "+"))
#minus_end <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(end_width, width=1, names=gene.name), "+"))
t2_range_thirty <-with(t2, GRanges(t2$chromosome_lowCu_ChIP, IRanges(start_width_30_lowCu_ChIP, end_width_30_lowCu_ChIP, names=chromosome_lowCu_ChIP), "+"))

#test <- mergeByOverlaps(minus_peak, plus_range)
#str(test)

#test_df <- data.frame(test)

t1_t2_olaps_peak_thirty <- findOverlaps(t1_range_thirty, t2_range_thirty)
#olaps_start <- findOverlaps(minus_start, plus_range)
#olaps_end <- findOverlaps(minus_end, plus_range)

queryHits(t1_t2_olaps_peak_thirty)
#queryHits(olaps_start)
#queryHits(olaps_end)

t1_t2_intersect_peak_thirty <- cbind(t1[queryHits(t1_t2_olaps_peak_thirty),], t2[subjectHits(t1_t2_olaps_peak_thirty),])
nrow(t1_t2_intersect_peak_thirty)
head(t1_t2_intersect_peak_thirty)
count(as.data.frame(table(t1_t2_intersect_peak_thirty$peak_name_controlCu_ChIP)))
sum(duplicated(t1_t2_intersect_peak_thirty$peak_name_controlCu_ChIP))
count(as.data.frame(table(t1_t2_intersect_peak_thirty$peak_name_lowCu_ChIP)))
sum(duplicated(t1_t2_intersect_peak_thirty$peak_name_lowCu_ChIP))
#colnames(t1_t2_intersect_peak_thirty) <- c("chromosome_WT", "start_WT", "end_WT", "peak_name_controlCuWT","width_WT", "peak_center_WT", "width_30", "start_width_30", "end_width_30", "chromosome_ChIP", "start_ChIP", "end_ChIP" ,"peak_name_controlCu_ChIP", "replicates_ChIP", "width_ChIP", "peak_center_ChIP", "width_30_ChIP", "start_width_30_ChIP", "end_width_30_ChIP")
t3 <- t1_t2_intersect_peak_thirty %>% distinct(peak_name_controlCu_ChIP, .keep_all = TRUE)
count(as.data.frame(table(t3$peak_name_controlCu_ChIP)))
sum(duplicated(t3$peak_name_controlCu_ChIP))
count(as.data.frame(table(t3$peak_name_lowCu_ChIP)))
sum(duplicated(t3$peak_name_lowCu_ChIP))
nrow(t3)
head(t3)
write.table(t3, "/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/WT_intersect/Intersect_of_Intersects_control_lowCu_ChIP_WT.txt", sep = "\t", row.names = F)



           
###WT
WT_Cusuf_sub$peak_center <- (WT_Cusuf_sub$width/2) + WT_Cusuf_sub$start
WT_Cusuf_sub$width_thirty <- WT_Cusuf_sub$width * (0.3)
WT_Cusuf_sub$start_width_thirty <- WT_Cusuf_sub$peak_center - WT_Cusuf_sub$width_thirty
WT_Cusuf_sub$end_width_thirty <- WT_Cusuf_sub$peak_center + WT_Cusuf_sub$width_thirty
head(WT_Cusuf_sub)

##MINUS
WT_Cudef_sub$peak_center <- (WT_Cudef_sub$width/2) + WT_Cudef_sub$start
WT_Cudef_sub$width_thirty <- WT_Cudef_sub$width * (0.3)
WT_Cudef_sub$start_width_thirty <- WT_Cudef_sub$peak_center - WT_Cudef_sub$width_thirty
WT_Cudef_sub$end_width_thirty <- WT_Cudef_sub$peak_center + WT_Cudef_sub$width_thirty
head(WT_Cudef_sub)

##Merge within the range of start.width and end.width
WT_Cudef_sub_range_thirty <- with(WT_Cudef_sub, GRanges(WT_Cudef_sub$chromosome, IRanges(start_width_thirty, end_width_thirty, names=chromosome), "+"))
#minus_start <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(start_width, width=1, names=gene.name), "+"))
#minus_end <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(end_width, width=1, names=gene.name), "+"))
WT_Cusuf_sub_range_thirty <-with(WT_Cusuf_sub, GRanges(WT_Cusuf_sub$chromosome, IRanges(start_width_thirty, end_width_thirty, names=chromosome), "+"))

#test <- mergeByOverlaps(minus_peak, plus_range)
#str(test)

#test_df <- data.frame(test)

WT_olaps_peak_thirty <- findOverlaps(WT_Cudef_sub_range_thirty, WT_Cusuf_sub_range_thirty)
#olaps_start <- findOverlaps(minus_start, plus_range)
#olaps_end <- findOverlaps(minus_end, plus_range)

queryHits(WT_olaps_peak_thirty)
#queryHits(olaps_start)
#queryHits(olaps_end)

WT_olaps_peak_thirty <- cbind(WT_Cudef_sub[queryHits(WT_olaps_peak_thirty),], WT_Cusuf_sub[subjectHits(WT_olaps_peak_thirty),])
nrow(WT_olaps_peak_thirty)
head(WT_olaps_peak_thirty)
write.table(WT_olaps_peak_thirty, "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/l2FC05/RAW_ChIP_intersect_peak_thirty.txt", sep = "\t", row.names = F)

###Intersect vs Intersect
ChIP_intersect_peak_thirty
nrow(ChIP_intersect_peak_thirty)
head(ChIP_intersect_peak_thirty)

names(ChIP_intersect_peak_thirty)[10] <- "start_width_thirty_minus"
names(ChIP_intersect_peak_thirty)[11] <- "end_width_thirty_minus"
names(ChIP_intersect_peak_thirty)[22] <- "start_width_thirty_plus"
names(ChIP_intersect_peak_thirty)[23] <- "end_width_thirty_plus"
names(ChIP_intersect_peak_thirty)[4] <- "peak_name_minus_ChIP"
names(ChIP_intersect_peak_thirty)[15] <- "peak_name_plus_ChIP"

head(WT_olaps_peak_thirty)
names(WT_olaps_peak_thirty)[8] <- "start_width_thirty_minus"
names(WT_olaps_peak_thirty)[9] <- "end_width_thirty_minus"
names(WT_olaps_peak_thirty)[17] <- "start_width_thirty_plus"
names(WT_olaps_peak_thirty)[18] <- "end_width_thirty_plus"
names(WT_olaps_peak_thirty)[4] <- "peak_name_minus_WT"
names(WT_olaps_peak_thirty)[13] <- "peak_name_plus_WT"
##Merge within the range of start.width and end.width
WT_inter_minus_range_thirty <- with(WT_olaps_peak_thirty, GRanges(WT_olaps_peak_thirty$chromosome, IRanges(start_width_thirty_minus, end_width_thirty_minus, names=chromosome), "+"))
WT_inter_plus_range_thirty <- with(WT_olaps_peak_thirty, GRanges(WT_olaps_peak_thirty$chromosome, IRanges(start_width_thirty_plus, end_width_thirty_plus, names=chromosome), "+"))
#minus_start <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(start_width, width=1, names=gene.name), "+"))
#minus_end <- with(ChIP_minus_sub, GRanges("gene.name", IRanges(end_width, width=1, names=gene.name), "+"))
ChIP_inter_minus_range_thirty <-with(ChIP_intersect_peak_thirty, GRanges(ChIP_intersect_peak_thirty$chromosome, IRanges(start_width_thirty_minus, end_width_thirty_minus, names=chromosome), "+"))
ChIP_inter_plus_range_thirty <-with(ChIP_intersect_peak_thirty, GRanges(ChIP_intersect_peak_thirty$chromosome, IRanges(start_width_thirty_plus, end_width_thirty_plus, names=chromosome), "+"))
#test <- mergeByOverlaps(minus_peak, plus_range)
#str(test)

#test_df <- data.frame(test)

WT_ChIP_inter_minus_peak_thirty <- findOverlaps(WT_inter_minus_range_thirty, ChIP_inter_minus_range_thirty)
WT_ChIP_inter_plus_peak_thirty <- findOverlaps(WT_inter_plus_range_thirty, ChIP_inter_plus_range_thirty)
WT_ChIP_inter_mix1_peak_thirty <- findOverlaps(WT_inter_plus_range_thirty, ChIP_inter_minus_range_thirty)
WT_ChIP_inter_mix2_peak_thirty <- findOverlaps(WT_inter_minus_range_thirty, ChIP_inter_plus_range_thirty)
#olaps_start <- findOverlaps(minus_start, plus_range)
#olaps_end <- findOverlaps(minus_end, plus_range)

queryHits(WT_ChIP_inter_minus_peak_thirty)
subsetByOverlaps(WT_inter_minus_range_thirty, ChIP_inter_minus_range_thirty)
queryHits(WT_ChIP_inter_plus_peak_thirty)
queryHits(WT_ChIP_inter_mix1_peak_thirty)
queryHits(WT_ChIP_inter_mix2_peak_thirty)
#queryHits(olaps_start)
#queryHits(olaps_end)

WT_ChIP_inter_minus_peak_thirty_df <- cbind(WT_olaps_peak_thirty[queryHits(WT_ChIP_inter_minus_peak_thirty),], ChIP_intersect_peak_thirty[subjectHits(WT_ChIP_inter_minus_peak_thirty),])
nrow(WT_ChIP_inter_minus_peak_thirty_df)
head(WT_ChIP_inter_minus_peak_thirty_df)
count(as.data.frame(table(WT_ChIP_inter_minus_peak_thirty_df$peak_name_minus_ChIP)))
count(as.data.frame(table(WT_ChIP_inter_minus_peak_thirty_df$peak_name_minus_WT)))
count(as.data.frame(table(WT_ChIP_inter_minus_peak_thirty_df$peak_name_plus_WT)))
sum(duplicated(WT_ChIP_inter_minus_peak_thirty_df$peak_name_minus_ChIP))
sum(duplicated(WT_ChIP_inter_minus_peak_thirty_df$peak_name_minus_WT))
write.table(WT_ChIP_inter_minus_peak_thirty, "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/l2FC05/RAW_ChIP_intersect_peak_thirty.txt", sep = "\t", row.names = F)

WT_ChIP_inter_plus_peak_thirty_df <- cbind(WT_olaps_peak_thirty[queryHits(WT_ChIP_inter_plus_peak_thirty),], ChIP_intersect_peak_thirty[subjectHits(WT_ChIP_inter_plus_peak_thirty),])
nrow(WT_ChIP_inter_plus_peak_thirty_df)
head(WT_ChIP_inter_plus_peak_thirty_df)
count(as.data.frame(table(WT_ChIP_inter_plus_peak_thirty_df$peak_name_plus_ChIP)))
count(as.data.frame(table(WT_ChIP_inter_plus_peak_thirty_df$peak_name_minus_ChIP)))
count(as.data.frame(table(WT_ChIP_inter_plus_peak_thirty_df$peak_name_plus_WT)))
count(as.data.frame(table(WT_ChIP_inter_plus_peak_thirty_df$peak_name_minus_WT)))
sum(duplicated(WT_ChIP_inter_plus_peak_thirty_df$peak_name_plus_ChIP))
sum(duplicated(WT_ChIP_inter_plus_peak_thirty_df$peak_name_plus_WT))
write.table(WT_ChIP_inter_minus_peak_thirty, "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/l2FC05/RAW_ChIP_intersect_peak_thirty.txt", sep = "\t", row.names = F)

WT_ChIP_inter_mix1_peak_thirty_df <- cbind(WT_olaps_peak_thirty[queryHits(WT_ChIP_inter_mix1_peak_thirty),], ChIP_intersect_peak_thirty[subjectHits(WT_ChIP_inter_mix1_peak_thirty),])
nrow(WT_ChIP_inter_mix1_peak_thirty_df)
head(WT_ChIP_inter_mix1_peak_thirty_df)
count(as.data.frame(table(WT_ChIP_inter_mix1_peak_thirty_df$peak_name_minus_ChIP)))
count(as.data.frame(table(WT_ChIP_inter_mix1_peak_thirty_df$peak_name_minus_WT)))
sum(duplicated(WT_ChIP_inter_mix1_peak_thirty_df$peak_name_minus_ChIP))
sum(duplicated(WT_ChIP_inter_mix1_peak_thirty_df$peak_name_minus_WT))
write.table(WT_ChIP_inter_minus_peak_thirty, "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/l2FC05/RAW_ChIP_intersect_peak_thirty.txt", sep = "\t", row.names = F)

WT_ChIP_inter_mix2_peak_thirty_df <- cbind(WT_olaps_peak_thirty[queryHits(WT_ChIP_inter_mix2_peak_thirty),], ChIP_intersect_peak_thirty[subjectHits(WT_ChIP_inter_mix2_peak_thirty),])
nrow(WT_ChIP_inter_mix2_peak_thirty_df)
head(WT_ChIP_inter_mix2_peak_thirty_df)
count(as.data.frame(table(WT_ChIP_inter_mix2_peak_thirty_df$peak_name_minus_ChIP)))
count(as.data.frame(table(WT_ChIP_inter_mix2_peak_thirty_df$peak_name_plus_WT)))
sum(duplicated(WT_ChIP_inter_mix2_peak_thirty_df$peak_name_minus_ChIP))
sum(duplicated(WT_ChIP_inter_mix2_peak_thirty_df$peak_name_plus_WT))
write.table(WT_ChIP_inter_minus_peak_thirty, "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/Integrative_Analysis_CHiP_RNAseq/l2FC05/RAW_ChIP_intersect_peak_thirty.txt", sep = "\t", row.names = F)


diff3 <- WT_ChIP_inter_plus_peak_thirty_df[!WT_ChIP_inter_plus_peak_thirty_df$peak_name_plus_ChIP%in%WT_ChIP_inter_minus_peak_thirty_df$peak_name_plus_ChIP,]
diff4 <- WT_ChIP_inter_plus_peak_thirty_df[!WT_ChIP_inter_plus_peak_thirty_df$peak_name_minus_ChIP%in%WT_ChIP_inter_minus_peak_thirty_df$peak_name_minus_ChIP,]

head(diff3)
diffx <- diff3[,c(4,8,9,13,17,18,22,28,29,33,40,41)]
diffy <- diff4[,c(4,8,9,13,17,18,22,28,29,33,40,41)]

######################
t1 ###control
t2 ###lowCu
t3 ###intersect
minus_anno <- read.delim("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/minus.merged.bed.associated.annotated.edited.txt")
minus_cut <- minus_anno[,c(5,6,7)]
plus_anno <- read.delim("/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/plus.merged.bed.associated.annotated.edit.txt" )
plus_cut <- plus_anno[,c(5,6,7)]
head(plus_cut)
anno <- read.delim("/Volumes/JARVIS-X/PostDoc/Plant_Genomic_Ressources/00 Functional Description for paper/TAIR10_functional_descriptions_for_Meetings_Mapman_strict_genes_202008_03_bHLHadded.txt", header = T)
names(anno)[1] <- "AGI"
head(anno)

head(t1)
t1_anno <- merge(t1,plus_cut, by.x = "peak_name_controlCu_ChIP", by.y="peak")
head(t1_anno)

head(t2)
t2_anno <- merge(t2,minus_cut, by.x = "peak_name_lowCu_ChIP", by.y="peak")
head(t2_anno)

head(t3)
t3_anno <- merge(t3,plus_cut, by.x = "peak_name_controlCu_ChIP", by.y="peak")
head(t3_anno)

t1_anno2 <- merge(t1_anno, anno,by.x ="gene.name", by.y="AGI" )
head(t1_anno2)

t2_anno2 <- merge(t2_anno, anno,by.x ="gene.name", by.y="AGI" )
head(t2_anno2)

t3_anno2 <- merge(t3_anno, anno,by.x ="gene.name", by.y="AGI" )
head(t3_anno2)

write.table(t1_anno2, "/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/WT/WT_Cusuf_intersect_peak_thirty_anno2.txt", sep = "\t", row.names = F)
write.table(t2_anno2, "/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/WT/WT_ChIP_Cudef_intersect_peak_thirty_anno2.txt", sep = "\t", row.names = F)
write.table(t3_anno2, "/Volumes/JARVIS-X/PostDoc/01_Projects/SPL7 - Anna/ChIP/WT/WT_ChIP_inter_plus_peak_thirty_df_anno2.txt", sep = "\t", row.names = F)

