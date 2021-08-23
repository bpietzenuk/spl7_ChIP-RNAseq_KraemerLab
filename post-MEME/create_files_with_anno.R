setwd("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/01_Motif_analysis/filters_gff/")
anno <- read.delim("/Volumes/JARVIS-X/PostDoc/Plant_Genomic_Ressources/00 Functional Description for paper/TAIR10_functional_descriptions_for_Meetings_Mapman_strict_genes_202008_03_bHLHadded.txt", header = T)

filelist <- list.files(pattern="*.noHeader.txt")

for (i in 1:length(filelist)) 
  {assign(filelist[i], read.table(text = gsub (";","\t", readLines(filelist[i])),header = F, fill = T))
}

for (i in 1:length(filelist))
{assign(paste(filelist[i],"_d"),get(filelist[i])[,c(13,11,12,8,1:3,6,5)])
} 

for (i in 1:length(filelist))
{assign(paste(filelist[i],"_anno"),merge(get(paste(filelist[i],"_d")), anno,by.x = "V5", by.y = "Ath_ID"))
}

for (i in 1:length(filelist))
       {write.table(get(paste(filelist[i],"_anno")), paste(filelist[i], "_anno.txt", sep = ""), sep = "\t", row.names = F, quote = F)}
