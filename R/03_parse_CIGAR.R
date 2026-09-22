library(reshape2)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

df <- read.table(paste0(args[1],"/ReadCIGAR.tsv"))
cigar <- df$V2

dat <- data.frame()
del_c <- c()
for (v in 1:length(cigar)) {
  matches <- gregexpr("(\\d+)([MIDNSHP=X])", cigar[v], perl = TRUE)
  ops <- regmatches(cigar[v], matches)[[1]]
  
  #Get large indel
  del <- paste0(c(3:100),"D")
  if (any(del %in% ops)) {
    del_c <- c(del_c,df$V1[v])
  }
  
  # Split into numbers and letters
  counts <- as.numeric(sub("([0-9]+)([MIDNSHP=X])", "\\1", ops))
  types <- sub("([0-9]+)([MIDNSHP=X])", "\\2", ops)
  
  # Summarize
  summary <- tapply(counts, types, sum)
  dat_ <- as.data.frame(summary)
  dat_$type <- rownames(dat_)  
  rownames(dat_) <- NULL
  dat_$read <- df$V1[v]
  dat <- rbind(dat,dat_)
}

dat |> dcast(read~type,value.var = "summary") -> dat_final
dat_final[is.na(dat_final)] <- 0
dat_final$total <- dat_final[, !(names(dat_final) %in% c("read","D"))] |> apply(MARGIN = 1,FUN = sum)

dat_final$Sprop <- dat_final$S/dat_final$total
dat_final$Dprop <- dat_final$D/dat_final$total
dat_final$Edel <- ifelse(dat_final$read %in% del_c,"fail","pass")
dat_final$res <- ifelse(dat_final$Edel == 'fail' | dat_final$Sprop >= 0.2 | dat_final$Dprop >= 0.3,"exclude","keep")


dat_final |> write.table(paste0(args[1],"/CIGAR_Tab.tsv"),quote = F,row.names = F,col.names = T) 