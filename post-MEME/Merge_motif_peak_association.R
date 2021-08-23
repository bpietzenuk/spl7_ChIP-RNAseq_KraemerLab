#####Quick initial file modification
####Add peak occourances to modified peak file WITHOUT functional annotation
minus.associated <- read.table("../../minus.merged.bed.associated.modified", header=T)
minus.occourances <- read.table("../../minus.peak_occurances.txt", header = F)
minus.modified <- merge(minus.associated,minus.occourances, by.x="feature", by.y = "V1")
colnames(minus.modified)[8] <- "occurance"
write.table(minus.modified,"../../minus.merged.bed.associated.modified.score", row.names = F, sep = "\t")

plus.associated <- read.table("../../plus.merged.bed.associated.modified", header = T)
plus.occourances <- read.table("../../plus.peak_occurances.txt", header = F)
plus.modified <- merge(plus.associated,plus.occourances, by.x="feature", by.y = "V1")
colnames(plus.modified)[8] <- "occurance"
write.table(plus.modified,"../../plus.merged.bed.associated.modified.score", row.names = F, sep = "\t")

###Load Files into environment
##General files 
setwd("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/ChIP/01_Motif_analysis/candidates_gff/")

minus.associated <- read.table("../../minus.merged.bed.centered.associated.modified", header=T)
plus.associated <- read.table("../../plus.merged.bed.centered.associated.modified", header = T)

###candidates Integrated Analysis
AGAGAAGA <- read.table("Ara_memechip_spl7_Cudependent_Cu-_active.txt.AGAGAAGA;.gff", header = F)
GTACTGC <- read.table("Ara_memechip_spl7_Cudependent_Cu-_active.txt.GTACTGC;.gff", header = F)

AGAGAAGA.associated <- merge(minus.associated,AGAGAAGA, by.x ="feature", by.y = "V1")
GTACTGC.associated <- merge(minus.associated,GTACTGC, by.x ="feature", by.y = "V1")

AGAGAAGA.associated$motif.start.genome.centered <- (AGAGAAGA.associated$start.center + AGAGAAGA.associated$V4)
AGAGAAGA.associated$motif.end.genome.centered <- (AGAGAAGA.associated$motif.start.genome.centered + ((AGAGAAGA.associated$V5 - AGAGAAGA.associated$V4)+1))

GTACTGC.associated$motif.start.genome.centered <- (GTACTGC.associated$start.center + GTACTGC.associated$V4)
GTACTGC.associated$motif.end.genome.centered <- (GTACTGC.associated$motif.start.genome.centered + ((GTACTGC.associated$V5 - GTACTGC.associated$V4)+1))

AGAGAAGA.associated.modified <- AGAGAAGA.associated[,c(2,18,19,1,5:7,17,8,9)]
AGAGAAGA.associated.gff <- AGAGAAGA.associated[,c(2,1,5,18,19,6,15,16,17)]
AGAGAAGA.associated.gff$chromosome <- str_remove_all(AGAGAAGA.associated.gff$chromosome, "Chr")

GTACTGC.associated.modified <- GTACTGC.associated[,c(2,18,19,1,5:7,17,8,9)]
GTACTGC.associated.gff <- GTACTGC.associated[,c(2,1,5,18,19,6,15,16,17)]
GTACTGC.associated.gff$chromosome <- str_remove_all(GTACTGC.associated.gff$chromosome, "Chr")

write.table(AGAGAAGA.associated.modified, "../result_data/AGAGAAGA.associated.modified", row.names = F, sep = "\t", quote = F)

write.table(AGAGAAGA.associated.gff, "../result_data/AGAGAAGA.associated.gff", row.names = F,col.names = F, sep = "\t", quote = F)

write.table(GTACTGC.associated.modified, "../result_data/GTACTGC.associated.modified", row.names = F, sep = "\t", quote = F)

write.table(GTACTGC.associated.gff, "../result_data/GTACTGC.associated.gff", row.names = F,col.names = F, sep = "\t", quote = F)

###additional filters
####ChIPsum O-Q
plus.associated$feature_plus <- paste0(plus.associated$feature, "_plus") 
minus.associated$feature_minus <- paste0(minus.associated$feature, "_minus") 

O_corrected_ACGCGCTACC <- read.table("O_corrected.dedup.txt_centered.merged.named.fasta.ACGCGCTACC;.gff", header = F)
O_corrected_ACGCGCTACC.associated.minus <- merge(O_corrected_ACGCGCTACC,minus.associated, by.x ="V1", by.y = "feature_minus", all.x = T)
O_corrected_ACGCGCTACC.associated.plus <- merge(O_corrected_ACGCGCTACC.associated.minus,plus.associated, by.x ="V1", by.y = "feature_plus", all.x = T)
O_corrected_ACGCGCTACC.associated.plus = as.data.frame(t(apply(O_corrected_ACGCGCTACC.associated.plus,1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )))
O_corrected_ACGCGCTACC.associated.plus <- O_corrected_ACGCGCTACC.associated.plus[,c(1:18)]
O_corrected_ACGCGCTACC.associated.plus$V12 <- as.numeric(as.character(O_corrected_ACGCGCTACC.associated.plus$V12))
O_corrected_ACGCGCTACC.associated.plus$V4 <- as.numeric(as.character(O_corrected_ACGCGCTACC.associated.plus$V4))
O_corrected_ACGCGCTACC.associated.plus$V5 <- as.numeric(as.character(O_corrected_ACGCGCTACC.associated.plus$V5))

O_corrected_ACGCGCTACC.associated.plus$motif.start.genome.centered <- (O_corrected_ACGCGCTACC.associated.plus$V12 + O_corrected_ACGCGCTACC.associated.plus$V4)
O_corrected_ACGCGCTACC.associated.plus$motif.end.genome.centered <- (O_corrected_ACGCGCTACC.associated.plus$motif.start.genome.centered + ((O_corrected_ACGCGCTACC.associated.plus$V5 - O_corrected_ACGCGCTACC.associated.plus$V4)+1))

O_corrected_ACGCGCTACC.associated.plus.modified <- O_corrected_ACGCGCTACC.associated.plus[,c(11,19,20,1,14,15,16,9,17,18)]
colnames(O_corrected_ACGCGCTACC.associated.plus.modified) = c("chromosome",	"motif.start.genome.centered",	"peak_name",	"associated.gene",	"distance.to.site",	"occourance",	"attribute",	"original.peak.start",	"original.peak.end")
O_corrected_ACGCGCTACC.associated.plus.gff <- O_corrected_ACGCGCTACC.associated.plus[,c(11,1,14,19,20,15,7,8,9)]
colnames(O_corrected_ACGCGCTACC.associated.plus.gff) = c("chromosome",	"peak_name",	"associated.gene",	"motif.start.genome.centered",	"motif.end.genome.centered",	"distance.to.site",	"strand",	"phase",	"attribute")
O_corrected_ACGCGCTACC.associated.plus.gff$chromosome <- str_remove_all(O_corrected_ACGCGCTACC.associated.plus.gff$chromosome, "Chr")

write.table(O_corrected_ACGCGCTACC.associated.plus.modified, "O_corrected_ACGCGCTACC.associated.plus.modified", row.names = F, sep = "\t", quote = F)
write.table(O_corrected_ACGCGCTACC.associated.plus.gff, "O_corrected_ACGCGCTACC.associated.plus.gff", row.names = F,col.names = F, sep = "\t", quote = F)

O_corrected_CGCGCTACCA <- read.table("O_corrected.dedup.txt_centered.merged.named.fasta.CGCGCTACCA;.gff", header = F)
O_corrected_CGCGCTACCA.associated.minus <- merge(O_corrected_CGCGCTACCA,minus.associated, by.x ="V1", by.y = "feature_minus", all.x = T)
O_corrected_CGCGCTACCA.associated.plus <- merge(O_corrected_CGCGCTACCA.associated.minus,plus.associated, by.x ="V1", by.y = "feature_plus", all.x = T)
O_corrected_CGCGCTACCA.associated.plus = as.data.frame(t(apply(O_corrected_CGCGCTACCA.associated.plus,1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )))
O_corrected_CGCGCTACCA.associated.plus <- O_corrected_CGCGCTACCA.associated.plus[,c(1:18)]
O_corrected_CGCGCTACCA.associated.plus$V12 <- as.numeric(as.character(O_corrected_CGCGCTACCA.associated.plus$V12))
O_corrected_CGCGCTACCA.associated.plus$V4 <- as.numeric(as.character(O_corrected_CGCGCTACCA.associated.plus$V4))
O_corrected_CGCGCTACCA.associated.plus$V5 <- as.numeric(as.character(O_corrected_CGCGCTACCA.associated.plus$V5))

O_corrected_CGCGCTACCA.associated.plus$motif.start.genome.centered <- (O_corrected_CGCGCTACCA.associated.plus$V12 + O_corrected_CGCGCTACCA.associated.plus$V4)
O_corrected_CGCGCTACCA.associated.plus$motif.end.genome.centered <- (O_corrected_CGCGCTACCA.associated.plus$motif.start.genome.centered + ((O_corrected_CGCGCTACCA.associated.plus$V5 - O_corrected_CGCGCTACCA.associated.plus$V4)+1))

O_corrected_CGCGCTACCA.associated.plus.modified <- O_corrected_CGCGCTACCA.associated.plus[,c(11,19,20,1,14,15,16,9,17,18)]
colnames(O_corrected_CGCGCTACCA.associated.plus.modified) = c("chromosome",	"motif.start.genome.centered",	"peak_name",	"associated.gene",	"distance.to.site",	"occourance",	"attribute",	"original.peak.start",	"original.peak.end")
O_corrected_CGCGCTACCA.associated.plus.gff <- O_corrected_CGCGCTACCA.associated.plus[,c(11,1,14,19,20,15,7,8,9)]
colnames(O_corrected_CGCGCTACCA.associated.plus.gff) = c("chromosome",	"peak_name",	"associated.gene",	"motif.start.genome.centered",	"motif.end.genome.centered",	"distance.to.site",	"strand",	"phase",	"attribute")
O_corrected_CGCGCTACCA.associated.plus.gff$chromosome <- str_remove_all(O_corrected_CGCGCTACCA.associated.plus.gff$chromosome, "Chr")

write.table(O_corrected_CGCGCTACCA.associated.plus.modified, "O_corrected_CGCGCTACCA.associated.plus.modified", row.names = F, sep = "\t", quote = F)
write.table(O_corrected_CGCGCTACCA.associated.plus.gff, "O_corrected_CGCGCTACCA.associated.plus.gff", row.names = F,col.names = F, sep = "\t", quote = F)

P_uncorrected_CGTGGAGTTG <- read.table("P_uncorrected.dedup.txt_centered.merged.named.fasta.CGTGGAGTTG;.gff", header = F)
P_uncorrected_CGTGGAGTTG.associated.minus <- merge(P_uncorrected_CGTGGAGTTG,minus.associated, by.x ="V1", by.y = "feature_minus", all.x = T)
P_uncorrected_CGTGGAGTTG.associated.plus <- merge(P_uncorrected_CGTGGAGTTG.associated.minus,plus.associated, by.x ="V1", by.y = "feature_plus", all.x = T)
P_uncorrected_CGTGGAGTTG.associated.plus = as.data.frame(t(apply(P_uncorrected_CGTGGAGTTG.associated.plus,1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )))
P_uncorrected_CGTGGAGTTG.associated.plus <- P_uncorrected_CGTGGAGTTG.associated.plus[,c(1:18)]
P_uncorrected_CGTGGAGTTG.associated.plus$V12 <- as.numeric(as.character(P_uncorrected_CGTGGAGTTG.associated.plus$V12))
P_uncorrected_CGTGGAGTTG.associated.plus$V4 <- as.numeric(as.character(P_uncorrected_CGTGGAGTTG.associated.plus$V4))
P_uncorrected_CGTGGAGTTG.associated.plus$V5 <- as.numeric(as.character(P_uncorrected_CGTGGAGTTG.associated.plus$V5))

P_uncorrected_CGTGGAGTTG.associated.plus$motif.start.genome.centered <- (P_uncorrected_CGTGGAGTTG.associated.plus$V12 + P_uncorrected_CGTGGAGTTG.associated.plus$V4)
P_uncorrected_CGTGGAGTTG.associated.plus$motif.end.genome.centered <- (P_uncorrected_CGTGGAGTTG.associated.plus$motif.start.genome.centered + ((P_uncorrected_CGTGGAGTTG.associated.plus$V5 - P_uncorrected_CGTGGAGTTG.associated.plus$V4)+1))

P_uncorrected_CGTGGAGTTG.associated.plus.modified <- P_uncorrected_CGTGGAGTTG.associated.plus[,c(11,19,20,1,14,15,16,9,17,18)]
colnames(P_uncorrected_CGTGGAGTTG.associated.plus.modified) = c("chromosome",	"motif.start.genome.centered",	"peak_name",	"associated.gene",	"distance.to.site",	"occourance",	"attribute",	"original.peak.start",	"original.peak.end")
P_uncorrected_CGTGGAGTTG.associated.plus.gff <- P_uncorrected_CGTGGAGTTG.associated.plus[,c(11,1,14,19,20,15,7,8,9)]
colnames(P_uncorrected_CGTGGAGTTG.associated.plus.gff) = c("chromosome",	"peak_name",	"associated.gene",	"motif.start.genome.centered",	"motif.end.genome.centered",	"distance.to.site",	"strand",	"phase",	"attribute")
P_uncorrected_CGTGGAGTTG.associated.plus.gff$chromosome <- str_remove_all(P_uncorrected_CGTGGAGTTG.associated.plus.gff$chromosome, "Chr")

write.table(P_uncorrected_CGTGGAGTTG.associated.plus.modified, "P_uncorrected_CGTGGAGTTG.associated.plus.modified", row.names = F, sep = "\t", quote = F)
write.table(P_uncorrected_CGTGGAGTTG.associated.plus.gff, "P_uncorrected_CGTGGAGTTG.associated.plus.gff", row.names = F,col.names = F, sep = "\t", quote = F)


P_uncorrected_TCGAACCCGG <- read.table("P_uncorrected.dedup.txt_centered.merged.named.fasta.TCGAACCCGG;.gff", header = F)
P_uncorrected_TCGAACCCGG.associated.minus <- merge(P_uncorrected_TCGAACCCGG,minus.associated, by.x ="V1", by.y = "feature_minus", all.x = T)
P_uncorrected_TCGAACCCGG.associated.plus <- merge(P_uncorrected_TCGAACCCGG.associated.minus,plus.associated, by.x ="V1", by.y = "feature_plus", all.x = T)
P_uncorrected_TCGAACCCGG.associated.plus = as.data.frame(t(apply(P_uncorrected_TCGAACCCGG.associated.plus,1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )))
P_uncorrected_TCGAACCCGG.associated.plus <- P_uncorrected_TCGAACCCGG.associated.plus[,c(1:18)]
P_uncorrected_TCGAACCCGG.associated.plus$V12 <- as.numeric(as.character(P_uncorrected_TCGAACCCGG.associated.plus$V12))
P_uncorrected_TCGAACCCGG.associated.plus$V4 <- as.numeric(as.character(P_uncorrected_TCGAACCCGG.associated.plus$V4))
P_uncorrected_TCGAACCCGG.associated.plus$V5 <- as.numeric(as.character(P_uncorrected_TCGAACCCGG.associated.plus$V5))

P_uncorrected_TCGAACCCGG.associated.plus$motif.start.genome.centered <- (P_uncorrected_TCGAACCCGG.associated.plus$V12 + P_uncorrected_TCGAACCCGG.associated.plus$V4)
P_uncorrected_TCGAACCCGG.associated.plus$motif.end.genome.centered <- (P_uncorrected_TCGAACCCGG.associated.plus$motif.start.genome.centered + ((P_uncorrected_TCGAACCCGG.associated.plus$V5 - P_uncorrected_TCGAACCCGG.associated.plus$V4)+1))

P_uncorrected_TCGAACCCGG.associated.plus.modified <- P_uncorrected_TCGAACCCGG.associated.plus[,c(11,19,20,1,14,15,16,9,17,18)]
colnames(P_uncorrected_TCGAACCCGG.associated.plus.modified) = c("chromosome",	"motif.start.genome.centered",	"peak_name",	"associated.gene",	"distance.to.site",	"occourance",	"attribute",	"original.peak.start",	"original.peak.end")
P_uncorrected_TCGAACCCGG.associated.plus.gff <- P_uncorrected_TCGAACCCGG.associated.plus[,c(11,1,14,19,20,15,7,8,9)]
colnames(P_uncorrected_TCGAACCCGG.associated.plus.gff) = c("chromosome",	"peak_name",	"associated.gene",	"motif.start.genome.centered",	"motif.end.genome.centered",	"distance.to.site",	"strand",	"phase",	"attribute")
P_uncorrected_TCGAACCCGG.associated.plus.gff$chromosome <- str_remove_all(P_uncorrected_TCGAACCCGG.associated.plus.gff$chromosome, "Chr")

write.table(P_uncorrected_TCGAACCCGG.associated.plus.modified, "P_uncorrected_TCGAACCCGG.associated.plus.modified", row.names = F, sep = "\t", quote = F)
write.table(P_uncorrected_TCGAACCCGG.associated.plus.gff, "P_uncorrected_TCGAACCCGG.associated.plus.gff", row.names = F,col.names = F, sep = "\t", quote = F)