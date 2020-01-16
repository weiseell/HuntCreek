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

i <- 1
i <- NULL
tmp <- NULL
for(i in 1:length(my.files)){
  df <- read.table(file = paste0("Input/",my.files[i]),header = T,sep = "\t",stringsAsFactors = F,strip.white = T)
  colnames(df) <- c("off","dad","mom")
  df <- select(df, off,mom,dad)
  df$off1 <- df$off
  df <- separate(data = df,col = off1,into = c("off1","tmom","tdad"),sep = "@")
  head(df)
  
  #true parents rs estimate
  tmom <- unique(df$tmom)
  tdad <- unique(df$tdad)
  df_t <- data.frame(par = c(tmom,tdad),
                     sex = c(rep("F",times=length(tmom)),rep("M",time=length(tdad))),
                     stringsAsFactors = F)
  head(df_t)
  j <- 1
  j <- NULL
  for(j in 1:nrow(df_t)){
  if(df_t$sex[j] == "F"){df_t$rs[j] <- nrow(df[df$tmom == df_t$par[j],])}
    if(df_t$sex[j] == "M"){df_t$rs[j] <- nrow(df[df$tdad == df_t$par[j],])}
  }
  head(df_t)
  
  #inferred parents rs estimate
  mom <- unique(df$mom)
  dad <- unique(df$dad)
  df_i <- data.frame(par = c(mom,dad),
                     sex = c(rep("F",times=length(mom)),rep("M",time=length(dad))),
                     stringsAsFactors = F)
  head(df_i)
  j <- 1
  j <- NULL
  for(j in 1:nrow(df_i)){
    if(df_i$sex[j] == "F"){df_i$rs[j] <- nrow(df[df$mom == df_i$par[j],])}
    if(df_i$sex[j] == "M"){df_i$rs[j] <- nrow(df[df$dad == df_i$par[j],])}
  }
  head(df_i)
  
  df_i$type <- "Inferred"
  df_t$type <- "True"
  df1 <- rbind(df_t,df_i)
  tmp1 <- df1 %>% group_by(type,sex) %>% summarize(means = round(mean(rs),2))
  tmp1$run <- my.files[i]
  tmp <- rbind(tmp,tmp1)
}
tmp

#write to file
write.table(tmp, file = "Output/accuracy.summary.txt",append = F, quote = F, sep = "\t",row.names = F, col.names = T)

#fin!




