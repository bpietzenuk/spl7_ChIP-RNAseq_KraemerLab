setwd("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/candidate_motifs/")
      
A <- read.table("spl7depact_lowCu.txt", header = F)
C_uncorrected <- read.table("spl7depact_controlCu_uncorrected.txt", header = F)
D_corrected <- read.table("spl7deprepr_controlCu_corrected.txt", header = F)

setwd("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/candidate_motifs/filters_G_Q/")
K_corrected <- read.table("K_corrected.txt", header = F)
K_uncorrected <- read.table("K_uncorrected.txt", header = F)
O_corrected <- read.table("O_corrected.dedup.txt", header = F)
P_corrected <- read.table("P_corrected.dedup.txt", header = F)
P_uncorrected <- read.table("P_uncorrected.dedup.txt", header = F)
Q_corrected <- read.table("Q_corrected.dedup.txt", header = F)
Q_uncorrected <- read.table("Q_uncorrected.dedup.txt", header = F)

setwd("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/01_Motif_analysis/filters_gff/filter_results/")

A_AGAGAAGA <- read.delim("spl7depact_lowCu.txt.minus.wIntersectpicos.AGAGAAGA.associated.modified.noHeader.txt_anno.txt", header = T)
A_GTACTGC <- read.delim("spl7depact_lowCu.txt.minus.wIntersectpicos.GTACTGC.associated.modified.noHeader.txt_anno.txt", header = T)
C_ATAACAAC <- read.delim("spl7depact_controlCu_uncorrected.txt.plus.wIntersectpicos.ATAACAAC.associated.modified.noHeader.txt_anno.txt", header = T)
D_GTGGGTTA <- read.delim("spl7deprepr_controlCu_corrected.txt.plus.wIntersectpicos.GTGGGTTA.associated.modified.noHeader.txt_anno.txt", header = T)
K_corrected_AGAGAAGA <- read.delim("K_corrected_AGAGAAGA.associated.modified.noHeader.txt_anno.txt", header = T)
K_uncorrected_AGAGAAGA <- read.delim("K_uncorrected_AGAGAAGA.associated.modified.noHeader.txt_anno.txt", header = T)
O_corrected_ACGCGCTACC <- read.delim("O_corrected_ACGCGCTACC.associated.plus.modified.noHeader.txt_anno.txt", header = T)
O_corrected_CGCGCTACCA <- read.delim("O_corrected_CGCGCTACCA.associated.plus.modified.noHeader.txt_anno.txt", header = T)
P_corrected_CGCGCTACCA <- read.delim("P_corrected_CGCGCTACCA.associated.plus.modified.noHeader.txt_anno.txt", header = T)
P_uncorrected_CGTGGAGTTG <- read.delim("P_uncorrected_CGTGGAGTTG.associated.plus.modified.noHeader.txt_anno.txt", header = T)
Q_corrected_CCATTGCGCC <- read.delim("Q_corrected_CCATTGCGCC.associated.plus.modified.noHeader.txt_anno.txt", header = T)
Q_corrected_GCTCAAATGG <- read.delim("Q_corrected_GCTCAAATGG.associated.plus.modified.noHeader.txt_anno.txt", header = T)
Q_corrected_GGGTTCGATT <- read.delim("Q_corrected_GGGTTCGATT.associated.plus.modified.noHeader.txt_anno.txt", header = T)
Q_uncorrected_GCTCAAATGG <- read.delim("Q_uncorrected_GCTCAAATGG.associated.plus.modified.noHeader.txt_anno.txt", header = T)

setwd("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/01_Motif_analysis/filters_gff/GQ/Last_candidates/")
O_corrected_CGCGCTACCA <- read.delim("O_corrected_CGCGCTACCA.associated.plus.modified.noHeader.txt_anno.txt", header = T)
O_uncorrected_TGGTAGCGCG <- read.delim("O_uncorrected_TGGTAGCGCG.associated.plus.modified.noHeader.txt_anno.txt", header = T)
P_uncorrected_GCTCAGTTAA <- read.delim("P_uncorrected_GCTCAGTTAA.associated.plus.modified.noHeader.txt_anno.txt", header = T)
Q_uncorrected_CGAGTCCCGG <- read.delim("Q_uncorrected_CGAGTCCCGG.associated.plus.modified.noHeader.txt_anno.txt", header = T)
Q_uncorrected_GCTCAAATGG <- read.delim("Q_uncorrected_GCTCAAATGG.associated.plus.modified.noHeader.txt_anno.txt", header = T)

#######Merge
A_AGAGAAGA_filter <- merge(A, A_AGAGAAGA, by.x = "V1", by.y = "V5")
A_GTACTGC_filter <- merge(A, A_GTACTGC, by.x = "V1", by.y = "V5")

C_ATAACAAC_filter <- merge(C_uncorrected, C_ATAACAAC, by.x = "V1", by.y = "V5")

D_GTGGGTTA_filter <- merge(D_corrected, D_GTGGGTTA, by.x = "V1", by.y = "V5")

K_corrected_AGAGAAGA_filter <- merge(K_corrected, K_corrected_AGAGAAGA, by.x = "V1", by.y = "V5")
K_uncorrected_AGAGAAGA_filter <- merge(K_uncorrected, K_uncorrected_AGAGAAGA, by.x = "V1", by.y = "V5")

O_corrected_ACGCGCTACC_filter <- merge(O_corrected, O_corrected_ACGCGCTACC, by.x = "V1", by.y = "V5")
O_corrected_CGCGCTACCA_filter <- merge(O_corrected, O_corrected_CGCGCTACCA, by.x = "V1", by.y = "V5")
O_corrected_CGCGCTACCA_filter <- merge(O_corrected, O_corrected_CGCGCTACCA, by.x = "V1", by.y = "V5")
O_uncorrected_TGGTAGCGCG_filter <- merge(O_corrected, O_uncorrected_TGGTAGCGCG, by.x = "V1", by.y = "V5")

P_corrected_CGCGCTACCA_filter <- merge(P_corrected, P_corrected_CGCGCTACCA, by.x = "V1", by.y = "V5")
P_uncorrected_CGTGGAGTTG_filter <- merge(P_uncorrected, P_uncorrected_CGTGGAGTTG, by.x = "V1", by.y = "V5")
P_uncorrected_GCTCAGTTAA_filter <- merge(P_uncorrected, P_uncorrected_GCTCAGTTAA, by.x = "V1", by.y = "V5")

Q_corrected_CCATTGCGCC_filter <- merge(Q_corrected, Q_corrected_CCATTGCGCC, by.x = "V1", by.y = "V5")
Q_corrected_GCTCAAATGG_filter <- merge(Q_corrected, Q_corrected_GCTCAAATGG, by.x = "V1", by.y = "V5")
Q_corrected_GGGTTCGATT_filter <- merge(Q_corrected, Q_corrected_GGGTTCGATT, by.x = "V1", by.y = "V5")
Q_uncorrected_CGAGTCCCGG_filter <- merge(Q_uncorrected, Q_uncorrected_CGAGTCCCGG, by.x = "V1", by.y = "V5")
Q_uncorrected_GCTCAAATGG_filter <- merge(Q_uncorrected, Q_uncorrected_GCTCAAATGG, by.x = "V1", by.y = "V5")

####save
write.table(A_AGAGAAGA_filter, "A_AGAGAAGA_filter.txt", col.names = F, row.names = F, sep = "\t")
write.table(A_GTACTGC_filter, "A_GTACTGC_filter.txt", col.names = F, row.names = F, sep = "\t")
write.table(C_ATAACAAC_filter, "C_ATAACAAC_filter.txt", col.names = F, row.names = F, sep = "\t")
write.table(D_GTGGGTTA_filter, "D_GTGGGTTA_filter.txt",col.names = F, row.names = F, sep = "\t")

write.table(K_corrected_AGAGAAGA_filter, "K_corrected_AGAGAAGA_filter.txt",col.names = F, row.names = F, sep = "\t")
write.table(K_uncorrected_AGAGAAGA_filter, "K_uncorrected_AGAGAAGA_filter.txt",col.names = F, row.names = F, sep = "\t")

write.table(O_corrected_ACGCGCTACC_filter, "O_corrected_ACGCGCTACC_filter.txt",col.names = F, row.names = F, sep = "\t")
write.table(O_corrected_CGCGCTACCA_filter, "O_corrected_CGCGCTACCA_filter.txt",col.names = F, row.names = F, sep = "\t")
write.table(O_corrected_CGCGCTACCA_filter, "O_corrected_CGCGCTACCA_filter.txt",col.names = F, row.names = F, sep = "\t")
write.table(O_uncorrected_TGGTAGCGCG_filter, "O_uncorrected_TGGTAGCGCG_filter.txt",col.names = F, row.names = F, sep = "\t")

write.table(P_corrected_CGCGCTACCA_filter, "P_corrected_CGCGCTACCA_filter.txt", col.names = F,row.names = F, sep = "\t")
write.table(P_uncorrected_CGTGGAGTTG_filter, "P_uncorrected_CGTGGAGTTG_filter.txt",col.names = F, row.names = F, sep = "\t")
write.table(P_uncorrected_GCTCAGTTAA_filter, "P_uncorrected_GCTCAGTTAA_filter.txt",col.names = F, row.names = F, sep = "\t")

write.table(Q_corrected_CCATTGCGCC_filter, "Q_corrected_CCATTGCGCC_filter.txt", col.names = F,row.names = F, sep = "\t")
write.table(Q_corrected_GCTCAAATGG_filter, "Q_corrected_GCTCAAATGG_filter.txt",col.names = F, row.names = F, sep = "\t")
write.table(Q_corrected_GGGTTCGATT_filter, "Q_corrected_GGGTTCGATT_filter.txt", col.names = F,row.names = F, sep = "\t")
write.table(Q_uncorrected_GCTCAAATGG_filter, "Q_uncorrected_GCTCAAATGG_filter.txt",col.names = F, row.names = F, sep = "\t")
write.table(Q_uncorrected_CGAGTCCCGG_filter, "Q_uncorrected_CGAGTCCCGG_filter.txt", col.names = F,row.names = F, sep = "\t")
write.table(Q_uncorrected_GCTCAAATGG_filter, "Q_uncorrected_GCTCAAATGG_filter.txt",col.names = F, row.names = F, sep = "\t")
