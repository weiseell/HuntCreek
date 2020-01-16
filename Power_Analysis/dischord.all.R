#originally created on: December 7th, 2017
#by: Nick Sard

#ABOUT: This script was written to 

#loading necessary libraries
library(tidyverse)

#setting working directory and reading in data
setwd("F:/Scribner_Lab/Hunt Creek Analysis/New Analysis/halfsib_sims/")

#loading my own functions
source("F:/Scribner_Lab/Hunt Creek Analysis/New Analysis/halfsib_sims/source.scripts/unique.pairs.R")

my.files <- list.files(path = "Input/",pattern = "BestConfig")
my.files <- unique.pairs(ids = my.files)
my.files$Var1 <- as.character(my.files$Var1)
my.files$Var2 <- as.character(my.files$Var2)
head(my.files)

i <- 1
i <- NULL
tmp <- NULL
for(i in 1:nrow(my.files)){
 
   #reading in files
  run1 <- read.table(paste0("Input/",my.files$Var1[i]),header = T,sep = "\t",stringsAsFactors = F)
  run2 <- read.table(paste0("Input/",my.files$Var2[i]),header = T,sep = "\t",stringsAsFactors = F)
  colnames(run1) <- c("off1", "dad1", "mom1")
  colnames(run2) <- c("off2", "dad2", "mom2")
  
  #####################
  ###run 1 and run 2###
  #####################
  df <- cbind(run1,run2)
  head(df)
  
  table(df$off1 == df$off2)
  df$dad.logic <- ifelse(df$dad1 == df$dad2,yes = 0, no = 1)
  df$mom.logic <- ifelse(df$mom1 == df$mom2,yes = 0, no = 1)
  head(df)
  
  dsum <- sum(df$dad.logic)
  msum <- sum(df$mom.logic)
  noff <- nrow(df)
  
  df1 <- data.frame(dad = dsum, mom = msum, noff = noff)
  df1$compare <- paste(my.files$Var1[i],my.files$Var2[i])
  tmp <- rbind(tmp, df1)
}

#calculating dischordance
tmp$disc.dad <- round(tmp$dad/tmp$noff,2)
tmp$disc.mom <- round(tmp$mom/tmp$noff,2)
tmp$disc <- round((tmp$dad + tmp$mom)/(tmp$noff * 2),2)
tmp$c2 <- gsub(pattern = "_run[0-9].txt",replacement = "",x = tmp$compare)
head(tmp)

#summary.stats
tmp1 <- tmp %>%
  group_by(c2) %>%
  summarize(count = n(),
            means = round(mean(disc),2),
            meds = median(disc),
            sds = round(sd(disc),2)) %>%
  separate(col = c2,into = c("r1","r2"),sep = " ") %>%
  mutate(both = r1 == r2) %>%
  arrange(both)
tmp
tmp1

#write to file
write.table(tmp, file = "Output/dischordance.long.txt",append = F, quote = F, sep = "\t",row.names = F, col.names = T)
write.table(tmp1, file = "Output/dischordance.summary.txt",append = F, quote = F, sep = "\t",row.names = F, col.names = T)

#fin!




