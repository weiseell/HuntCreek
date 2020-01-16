#originally created on: February 5th, 2018
#by: Ellie Weise

#ABOUT: This script was written to 

#loading necessary libraries
library(tidyverse)

#setting working directory and reading in data
setwd("H:/Scribner_Lab/Hunt Creek Analysis/New Analysis/halfsib_sims/")

#loading my own functions
# - none -

#getting the list of files
my.files <- list.files(path = "Input/",pattern = "BestConfig")
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
counts.dad <- NULL
counts.mom <- NULL
for(i in 1:length(my.files)){
  
  #reading in the pedigree and making columns for the real parents
  df <- read.table(file = paste0("Input/",my.files[i]),header = T,sep = "\t",stringsAsFactors = F,strip.white = T)
  colnames(df) <- c("off","dad","mom")
  df <- select(df, off,mom,dad)
  df$off1 <- df$off
  df <- separate(data = df,col = off1,into = c("off1","tmom","tdad"),sep = "@")
  head(df)
  
  #cleaning up df
  df$dad <- str_replace_all(string=df$dad, pattern="mp.dad.*", repl="UK")
  df$mom <- str_replace_all(string=df$mom, pattern="mp.mom.*", repl="UK")
  df$dad <- str_replace_all(string=df$dad, pattern="\\s", repl="")
  df$mom <- str_replace_all(string=df$mom, pattern="\\s", repl="")
  
  #counting the number of UK in df
  df1 <- as.data.frame(table(df$dad))
  head(df1)
  df2 <- df1 %>% 
    filter(Var1 == "UK")
  counts.dad <- rbind(counts.dad,df2)
  
  df1 <- as.data.frame(table(df$mom))
  head(df1)
  df2 <- df1 %>% 
    filter(Var1 == "UK")
  counts.mom <- rbind(counts.mom,df2)
}

mean.uk.dad <- sum(counts$Freq)/length(my.files)
mean.uk.mom <- sum(counts$Freq)/length(my.files)


