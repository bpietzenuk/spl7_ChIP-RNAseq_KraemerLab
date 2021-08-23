#!/usr/bin/env Rscript
####################################################################Into a pipe

# get the input passed from the shell script
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args[1]))
}

####run
## program...
df = read.table(args[1], header=TRUE)
names(df)[1]<-paste("count")  
names(df)[2]<-paste("peak")  
#names(args)[1] <- "count"
#names(args)[2] <- "peak"
plus <- read.delim("plus.merged.bed.associated.annotated.edit.txt")
test <- merge(plus, df, by = "peak", all.x = TRUE)