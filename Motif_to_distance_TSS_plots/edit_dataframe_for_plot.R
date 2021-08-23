args = commandArgs(trailingOnly=TRUE)

#anno <- read.delim("/Volumes/JARVIS-X/PostDoc/Plant_Genomic_Ressources/A_thaliana_Genome/01_Araport11_latest_release_2020_09_01/Araport11_genes_only_strandfitted_edited.txt", header = FALSE)
#anno <- read.delim("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/tair_genes_only_strandfitted_edited.txt", header = FALSE)
#anno_edit <- anno[,c(1,4,5,6,7)]
#colnames(anno_edit)[2] <- "gene.start"
#colnames(anno_edit)[3] <- "gene.end"
#colnames(anno_edit)[4] <- "gene.size"
#anno_edit$gene.size <- (anno_edit$gene.end - anno_edit$gene.start)+1

filter <- read.table( file = args[1])

filter_edit <- filter[,c(6,7,8,5,4)]
write.table(filter_edit, file = args[2], quote=FALSE, col.names = FALSE, row.names = FALSE)

#filter_edit

#filter_edit$V5 = paste("chr", filter_edit$V5, sep="")


#filter_edit <- merge(filter, anno_edit, by.x = "V1", by.y = "V7", all.x = T)
#write.table(filter_edit, file = args[2] )

#filter_edit$distance <- ifelse((filter_edit$V9=="includeFeature"|filter_edit$V9=="inside"), (filter_edit$V7 - filter_edit$gene.start)*2000/filter_edit$gene.size, filter_edit$V7 - filter_edit$gene.start)

#density <- density(filter_edit$distance)

#x <- density$x
#y <- density$y
#bw <- density$bw
#n <- density$n

#density_sub <- data.frame( x = x, y = y, bw = bw, n = n)
#A_GTACTaC_density_sub <- A_GTACTaC_density_sub[A_GTACTaC_density_sub$x>-3000,]
#write.csv( density_sub, file = args[3] )








