install.packages("phyper")
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#where q=size of overlap-1; m=number of upregulated genes in experiment #1;
#n=(total number of genes on platform-m); k=number of upregulated genes in experiment #2.

lowCu <- read.table("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/DAPseq comparisons/DAPseq_lowCu.txt",header = T)


##lowCu versus SPL1
q <- 231-1
m <- 2026
n <- 28497 - m
k <- 3731
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (3731/28497)
y=x*28497

##lowCu versus SPL5
q <- 264-1
m <- 2026
n <- 28497 - m
k <- 4301
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (4301/28497)
y=x*28497

##lowCu versus SPL9
q <- 809-1
m <- 2026
n <- 28497 - m
k <- 13831
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (13831/28497)
y=x*28497

##lowCu versus SPL9
q <- 83-1
m <- 2026
n <- 28497 - m
k <- 1061
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (1061/28497)
y=x*28497

##lowCu versus SPL14
q <- 167-1
m <- 2026
n <- 28497 - m
k <- 2263
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (2263/28497)
y=x*28497


controlCu <- read.table("/Volumes/JARVIS-X/PostDoc/Corona_Björn/SPL7 - Anna/DAPseq comparisons/DAPseq_controlCu.txt",header = T)
##controlCu versus SPL1
q <- -207
m <- 1901
n <- 28497 - m
k <- 3731
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (3731/28497)
y=x*28497

##controlCu versus SPL5
q <- 242-1
m <- 1901
n <- 28497 - m
k <- 4301
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (4301/28497)
y=x*28497

##controlCu versus SPL9
q <- 755-1
m <- 1901
n <- 28497 - m
k <- 13831
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (13831/28497)
y=x*28497

##controlCu versus SPL13
q <- 72-1
m <- 1901
n <- 28497 - m
k <- 1061
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (1061/28497)
y=x*28497

##controlCu versus SPL14
q <- 146-1
m <- 1901
n <- 28497 - m
k <- 2263
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (2263/28497)
y=x*28497

##lowCu versus SPL1_5_9_13_14
q <- 55-1
m <- 2026
n <- 28497 - m
k <- 672
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (672/28497)
y=x*28497

##control versus SPL1_5_9_13_14
q <- 50-1
m <- 1901
n <- 28497 - m
k <- 672
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (672/28497)
y=x*28497

##lowCu versus SPL1_5_9_14_15
q <- 123-1
m <- 2026
n <- 28497 - m
k <- 1665
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (1665/28497)
y=x*28497

##control versus SPL1_5_9_14_15
q <- 110-1
m <- 1901
n <- 28497 - m
k <- 1665
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (1665/28497)
y=x*28497


##lowCu versus SPL9_13_15
q <- 80-1
m <- 2026
n <- 28497 - m
k <- 1046
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (1046/28497)
y=x*28497

##control versus SPL9_13_15
q <- 69-1
m <- 1901
n <- 28497 - m
k <- 1046
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (1046/28497)
y=x*28497


###hOverlap DAPseq/our study vs SPL7 up, down, regulated
####lowcu SPL7 regulated
##lowCu
q <- 2-1
m <- 55
n <- 28497 - m
k <- 1264
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(55/28497) * (1264/28497)
y=x*28497

##controlCu
q <- 1-1
m <- 50
n <- 28497 - m
k <- 1264
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(50/28497) * (1264/28497)
y=x*28497


####lowcu SPL7 upregulated
##lowCu
q <- 0-1
m <- 55
n <- 28497 - m
k <- 675
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(55/28497) * (675/28497)
y=x*28497

##controlCu
q <- 0-1
m <- 50
n <- 28497 - m
k <- 675
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(50/28497) * (675/28497)
y=x*28497


####lowcu SPL7 down-regulated
##lowCu
q <- 2-1
m <- 55
n <- 28497 - m
k <- 590
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(55/28497) * (590/28497)
y=x*28497

##controlCu
q <- 1-1
m <- 50
n <- 28497 - m
k <- 590
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(50/28497) * (590/28497)
y=x*28497


###DAPseq SPLamp versus our study
##lowCu
##lowCu SPL1amp
q <- 377-1
m <- 2026
n <- 28497 - m
k <- 6829
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (6829/28497)
y=x*28497

##lowCu SPL5amp
q <- 815-1
m <- 2026
n <- 28497 - m
k <- 14598
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (14598/28497)
y=x*28497

##lowCu SPL9amp
q <- 820-1
m <- 2026
n <- 28497 - m
k <- 14075
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (14075/28497)
y=x*28497


##lowCu SPL15amp
q <- 780-1
m <- 2026
n <- 28497 - m
k <- 13631
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (13631/28497)
y=x*28497

##lowCu SPL1_5_9_15amp
q <- 321-1
m <- 2026
n <- 28497 - m
k <- 5862
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (5862/28497)
y=x*28497


##controlCu
##controlCu SPL1amp
q <- 339-1
m <- 1901
n <- 28497 - m
k <- 6829
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (6829/28497)
y=x*28497

##controlCu SPL5amp
q <- 760-1
m <- 1901
n <- 28497 - m
k <- 14598
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (14598/28497)
y=x*28497

##lowCu SPL9amp
q <- 756-1
m <- 1901
n <- 28497 - m
k <- 14075
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (14075/28497)
y=x*28497


##controlCu SPL15amp
q <- 725-1
m <- 1901
n <- 28497 - m
k <- 13631
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (13631/28497)
y=x*28497

##controlCu SPL1_5_9_15amp
q <- 282-1
m <- 1901
n <- 28497 - m
k <- 5862
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (5862/28497)
y=x*28497


##filter_amp
###SPL9_15_amp
##lowCu
q <- 706-1
m <- 2026
n <- 28497 - m
k <- 12217
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(2026/28497) * (12217/28497)
y=x*28497

##controlCu SPL1_5_9_15amp
q <- 649-1
m <- 1901
n <- 28497 - m
k <- 12217
phyper(q, m, n, k, lower.tail = F, log.p = FALSE)

#Grenzwert
x=(1901/28497) * (12217/28497)
y=x*28497

