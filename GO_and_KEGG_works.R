if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ChIPseeker","clusterProfiler", "AnnotationDbi", "org.At.tair.db" ))
install.packages("ggplot2")
# Load libraries
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(org.At.tair.db)
library(ggplot2)

ls("package:org.At.tair.db")
keytypes(org.At.tair.db)
head(keys(org.At.tair.db, keytype = "TAIR"))

#############ChIPall
infile_all <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP//test_GO_KEGG/all.merged.bed.associated.noHeader.txt"
all.table <- read.delim(infile_all, stringsAsFactors = F, header = F)
#DE.table$peak <- as.factor(DE.table$peak)
#DE.table_single <- data.frame()

#for(level%in%DE.table$peak)
#{DE.table_sub <-  DE.table[DE.table$peak == level,]
#DE.table_single <- rbind(DE.table_single, DE.table_sub[DE.table_sub$distance == min(DE.table_sub$distance),])

#}

all <- all.table$V9

str(all)

# Run GO enrichment analysis 
class(all)
ego_all <- enrichGO(gene = all, 
                     keyType = "TAIR", 
                     OrgDb = org.At.tair.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

# Dotplot visualization
barplot <- ggplot(ego_all, aes(x=Count, y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ego_all, showCategory=50)

ekegg_all <- enrichKEGG(gene = all,
                         organism = 'ath',
                         pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_all, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ekegg_all)
