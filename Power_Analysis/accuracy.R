#originally created on: December 7th, 2017
#by: Nick Sard

#ABOUT: This script was written to 

#loading necessary libraries
library(tidyverse)

#setting working directory and reading in data
setwd("C:/Users/knsar/Google Drive/R/Data analysis/2017/Hunt Creek/timing_runs_concordance4/")

#loading my own functions
# - none -

#getting the list of files
my.files <- list.files(path = "Input/",pattern = "run")
my.files


#reading in location information
locs <- read.table(file = "Input/Hunt_Creek_all_genotypes_all_loci.txt",header = T,sep = "\t",stringsAsFactors = F)
head(locs)

#getting rid of the offp
locs <- filter(locs,  type == "Adults")

#fixing names
table(locs$locs)
locs$loc[locs$loc == "EB_AuSable"] <- "ASR"
locs$loc[locs$loc == "Houghton_Creek"] <- "HGC"
locs$loc[locs$loc == "Gilchrist_Creek"] <- "GCC"
locs$loc[grepl(pattern = "Hunt_Creek",x = locs$loc)] <- "HTC"
table(locs$loc)

#getting just the columns I want
locs <- locs %>% select(loc,id)
head(locs)

i <- 1
i <- NULL
tmp <- NULL
for(i in 1:length(my.files)){
  df <- read.table(file = paste0("Input/",my.files[i]),header = T,sep = "\t",stringsAsFactors = F,strip.white = T)
  colnames(df) <- c("off","dad","mom")
  df <- select(df, off,mom,dad)
  df$off1 <- df$off
  df <- separate(data = df,col = off1,into = c("off1","tmom","tdad"),sep = "@")
  df$mom_test <- ifelse(test = df$mom == df$tmom,0,1)
  df$dad_test <- ifelse(test = df$dad == df$tdad, 0,1)
  df$mom_test[df$mom == "UK"] <- NA
  df$dad_test[df$dad == "UK"] <- NA
  df
  uk.mom <- length(df$mom_test[is.na(df$mom_test)])
  uk.dad <- length(df$dad_test[is.na(df$dad_test)])
  mom.disc <- sum(df$mom_test,na.rm = T)
  dad.disc <- sum(df$dad_test,na.rm = T)
  noff <- nrow(df)
  tmp1 <- data.frame(uk.mom,uk.dad,
                     mom.disc,dad.disc,
                     noff,stringsAsFactors = F)
  tmp1$mom.p.disc <- round((mom.disc/(noff-uk.mom)),2)
  tmp1$dad.p.disc <- round((dad.disc/(noff-uk.dad)),2)
  tmp1$file <- my.files[i]
  tmp <- rbind(tmp,tmp1)
}
tmp

#write to file
write.table(tmp, file = "Output/accuracy.summary.txt",append = F, quote = F, sep = "\t",row.names = F, col.names = T)

#fin!




