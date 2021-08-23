if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ChIPseeker","TxDb.Hsapiens.UCSC.hg19.knownGene", "EnsDb.Hsapiens.v75","clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "org.At.tair.db" ))
install.packages("ggplot2")
# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.At.tair.db)
library(ggplot2)

ls("package:org.At.tair.db")
keytypes(org.At.tair.db)
head(keys(org.At.tair.db, keytype = "TAIR"))

# Load data
samplefiles <- list.files("data/idr-bed", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Nanog", "Pou5f1")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)


plotAnnoBar(peakAnnoList)

plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")

nanog_annot <- data.frame(peakAnnoList[["Nanog"]]@anno)

# Get the entrez IDs
infile_minus <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/test_GO_KEGG/all.merged.bed.sepa.edit.txt"
all.table <- read.delim(infile_minus, stringsAsFactors = F, header = F)
all_ID <- all.table$V4


# Run GO enrichment analysis 
class(all_ID)
ego_all <- enrichGO(gene = all_ID, 
                keyType = "TAIR", 
                OrgDb = org.At.tair.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

#Visualization
barplot <- ggplot(ego_all, aes(x=Count, y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()

dotplot(ego_all, showCategory=50)

ekegg_all <- enrichKEGG(gene = all_ID,
                    organism = 'ath',
                    pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_all, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()

dotplot(ekegg_all)


#####plus
infile_plus <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP//test_GO_KEGG/plus.merged.bed.sepa.edit.txt"
plus.table <- read.delim(infile_plus, stringsAsFactors = F, header = F)
#DE.table$peak <- as.factor(DE.table$peak)
#DE.table_single <- data.frame()

#for(level%in%DE.table$peak)
#{DE.table_sub <-  DE.table[DE.table$peak == level,]
#DE.table_single <- rbind(DE.table_single, DE.table_sub[DE.table_sub$distance == min(DE.table_sub$distance),])
  
#}

plus <- plus.table$V4

str(plus)

# Run GO enrichment analysis 
class(plus)
ego_plus <- enrichGO(gene = plus, 
                keyType = "TAIR", 
                OrgDb = org.At.tair.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Dotplot visualization
barplot <- ggplot(ego_plus, aes(x=Count, y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ego_plus, showCategory=50)

ekegg_plus <- enrichKEGG(gene = plus,
                    organism = 'ath',
                    pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_plus, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ekegg_plus)

################################################
infile_minus <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/test_GO_KEGG/minus.merged.bed.associated"
minus.table <- read.delim(infile_minus, stringsAsFactors = F, header = T)
minus_ID <- minus.table$feature


# Run GO enrichment analysis 
class(minus_ID)
ego_minus <- enrichGO(gene = minus_ID, 
                      keyType = "TAIR", 
                      OrgDb = org.At.tair.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

#Visualization
barplot <- ggplot(ego_minus, aes(x=Count, y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()

dotplot(ego_minus, showCategory=50)

ekegg_minus <- enrichKEGG(gene = minus_ID,
                          organism = 'ath',
                          pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_minus, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()

dotplot(ekegg_minus)


#####plus
infile_plus <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP//test_GO_KEGG/plus.merged.bed.associated"
plus.table <- read.delim(infile_plus, stringsAsFactors = F, header = T)
#DE.table$peak <- as.factor(DE.table$peak)
#DE.table_single <- data.frame()

#for(level%in%DE.table$peak)
#{DE.table_sub <-  DE.table[DE.table$peak == level,]
#DE.table_single <- rbind(DE.table_single, DE.table_sub[DE.table_sub$distance == min(DE.table_sub$distance),])

#}

plus <- plus.table$feature

str(plus)

# Run GO enrichment analysis 
class(plus)
ego_plus <- enrichGO(gene = plus, 
                     keyType = "TAIR", 
                     OrgDb = org.At.tair.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

# Dotplot visualization
barplot <- ggplot(ego_plus, aes(x=Count, y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ego_plus, showCategory=50)

ekegg_plus <- enrichKEGG(gene = plus,
                         organism = 'ath',
                         pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_plus, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ekegg_plus)

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

#####all_sepa
infile_all <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP//test_GO_KEGG/all.merged.bed.sepa.edit.txt"
all.table <- read.delim(infile_all, stringsAsFactors = F, header = F)
#DE.table$peak <- as.factor(DE.table$peak)
#DE.table_single <- data.frame()

#for(level%in%DE.table$peak)
#{DE.table_sub <-  DE.table[DE.table$peak == level,]
#DE.table_single <- rbind(DE.table_single, DE.table_sub[DE.table_sub$distance == min(DE.table_sub$distance),])

#}

all <- all.table$V4

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
dotplot(ego_all, showCategory=20)

ekegg_all <- enrichKEGG(gene = all,
                        organism = 'ath',
                        pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_all, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ekegg_all)

#####SPL7_upregulated low Cu (TO DO)
infile_all <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP//test_GO_KEGG/spl7depact_lowCu.txt"
low_Cu_up.table <- read.delim(infile_all, stringsAsFactors = F, header = F)
#DE.table$peak <- as.factor(DE.table$peak)
#DE.table_single <- data.frame()

#for(level%in%DE.table$peak)
#{DE.table_sub <-  DE.table[DE.table$peak == level,]
#DE.table_single <- rbind(DE.table_single, DE.table_sub[DE.table_sub$distance == min(DE.table_sub$distance),])

#}

low_Cu_up <- low_Cu_up.table$V1

str(low_Cu_up)

# Run GO enrichment analysis 
class(low_Cu_up)
ego_low_Cu_up <- enrichGO(gene = low_Cu_up, 
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
dotplot(ego_low_Cu_up, showCategory=20)

ekegg_low_Cu_up <- enrichKEGG(gene = low_Cu_up,
                        organism = 'ath',
                        pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_all, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ekegg_low_Cu_up)

#####SPL7_repressed low Cu (TO DO)
infile_all <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP//test_GO_KEGG/spl7deprepr_lowCu.txt"
low_Cu_dwn.table <- read.delim(infile_all, stringsAsFactors = F, header = F)
#DE.table$peak <- as.factor(DE.table$peak)
#DE.table_single <- data.frame()

#for(level%in%DE.table$peak)
#{DE.table_sub <-  DE.table[DE.table$peak == level,]
#DE.table_single <- rbind(DE.table_single, DE.table_sub[DE.table_sub$distance == min(DE.table_sub$distance),])

#}

low_Cu_dwn <- low_Cu_dwn.table$V1

str(low_Cu_dwn)

# Run GO enrichment analysis 
class(low_Cu_dwn)
ego_low_Cu_dwn <- enrichGO(gene = low_Cu_dwn, 
                          keyType = "TAIR", 
                          OrgDb = org.At.tair.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = TRUE)

# Dotplot visualization
barplot <- ggplot(ego_low_Cu_dwn, aes(x=Count, y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ego_low_Cu_dwn, showCategory=1)

ekegg_low_Cu_dwn <- enrichKEGG(gene = low_Cu_dwn,
                              organism = 'ath',
                              pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_low_Cu_dwn, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ekegg_low_Cu_dwn)

#####SPL7_regulated low Cu (TO DO)
infile_all <- "/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP//test_GO_KEGG/spl7depregulated_lowCu.txt"
low_Cu_reg.table <- read.delim(infile_all, stringsAsFactors = F, header = F)
#DE.table$peak <- as.factor(DE.table$peak)
#DE.table_single <- data.frame()

#for(level%in%DE.table$peak)
#{DE.table_sub <-  DE.table[DE.table$peak == level,]
#DE.table_single <- rbind(DE.table_single, DE.table_sub[DE.table_sub$distance == min(DE.table_sub$distance),])

#}

low_Cu_reg <- low_Cu_reg.table$V1

str(low_Cu_reg)

# Run GO enrichment analysis 
class(low_Cu_reg)
ego_low_Cu_reg <- enrichGO(gene = low_Cu_reg, 
                           keyType = "TAIR", 
                           OrgDb = org.At.tair.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = TRUE)

# Dotplot visualization
barplot <- ggplot(ego_low_Cu_reg, aes(x=Count, y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ego_low_Cu_reg, showCategory=1)

ekegg_low_Cu_reg <- enrichKEGG(gene = low_Cu_reg,
                               organism = 'ath',
                               pvalueCutoff = 0.05)

barplot <- ggplot(ekegg_low_Cu_reg, aes(x=log10(pvalue), y=Description)) + 
  geom_bar(position=position_dodge(), stat="identity")
barplot + theme_classic()
dotplot(ekegg_low_Cu_reg)


