#!/usr/bin/env Rscript

#originally created on: December 7th, 2017
#by: Nick Sard and Ellie Weise

#ABOUT: This script was written to simulate offspring genotypes created by making from disparate (sig Fst diff) populations

### note ###
#need to run this before using this script
#module purge
#module use /opt/modules/Others/all/
#module load R/3.5.1-foss-2018a-X11-20180131
### note ###
 
#loading necessary libraries
library(related)
library(stringr)
library(tidyverse)

#setting working directory and reading in data
setwd("Translocation_Assessment/")

#loading my own functions
source("source/hunt.creek.functions.R")
source("source/share.R")
source("source/brd.mat.functions1.R")

#reading in the parental genotypes
df <- read.table("input/Hunt_Creek_all_genotypes_all_loci.txt",sep = "\t",header = T,stringsAsFactors = F)
#head(df)

#getting the count of offspring
noff <- nrow(df[df$type == "Offspring",])

#selecting just the adults
df <- df[df$type == "Adults",]
df$type <- NULL
df <- df[df$Sex != "UK",]
#head(df)

#using a for loop to prep the data for genepop format
#making two character genotypes into three and changing
#missing genotypes to 000
for(i in 4:ncol(df)){
  for(j in 1:nrow(df)){
    if(str_length(df[j,i]) == 2){
      df[j,i] <- paste0("0",df[j,i])
    }
    if(df[j,i] == 0){
      df[j,i] <- "000"
    }
  }
}
#head(df)

#saving the vector containing sex for later
my.sex.vector <- df$Sex
df$Sex <- NULL

#making genotypes into six character format
df <- six.char(population = df,missing.data = 0,one.pop = F)
#head(df)

#returning the sex vector
df$sex <- my.sex.vector
df <- df[,c(1,2,ncol(df),3:(ncol(df)-1))]

#renaming locations
df$loc[df$loc == "Hunt_Creek_Res"] <- "A"
df$loc[df$loc == "Houghton_Creek"] <- "B"
df$loc[df$loc == "EB_AuSable"] <- "C"
df$loc[df$loc == "Gilchrist_Creek"] <- "D"
df$loc[df$loc == "Hunt_Creek_TrP"] <- "A"
#table(df$loc)

#more renaming
names(df)[1] <- "Pops"
names(df)[2] <- "Sample.name"

#making sample names
df$Sample.name <- paste(df$Pops,df$Sample.name,sep="_")
#head(df)

#getting the genotypes in to six character states for moms and dads
moms <- df[df$sex == "F",c(-1,-3),]
#head(moms)

#dads
dads <- df[df$sex == "M",c(-1,-3),]
#head(dads)

#getting counts of each to make the breeding matrix
n.mom <- nrow(moms)
n.dad <- nrow(dads)

#making a for loop to do the get a distribution for expected random relatedness, A, and He estimates for a randomization
out <- NULL
i <- 1
i <- NULL
for(i in 1:250){

  #making the breeding matrix
  mat.str <- brd.mat(moms = n.mom,dads = n.dad,min.mates = 1,max.mates = 1,min.fert = 1,max.fert = 5)
  #head(mat.str)

  #getting summary stats
  ms1 <- mat.stats(mat = mat.str)
  #ms1

  #making that into a df and assigning colnames to moms and making a new col for dads
  tmp <- data.frame(mat.str)
  colnames(tmp) <- sample(x = moms[,1],size = length(moms[,1]),replace=F)
  tmp <- cbind(sample(x = dads[,1],size = length(dads[,1]),replace=F),tmp)
  names(tmp)[1] <-"dads"
  row.names(tmp) <- NULL
  #head(tmp)

  #making this into long form and removing all mate pairs that didnt make babies
  tmp1 <- gather(data = tmp,key = "moms",value = "noff",-dads) %>% filter(noff != 0)
  tmp1$dads <- as.character(tmp1$dads)
  #head(tmp1)

  #for each mate pair I will produce offspring genotypes using the laws of independent segregation and assortment, then converting back to
  #three character format
  off <- make.gts(tmp1 = tmp1,moms = moms,dads = dads)
  off <- three.char(tmp = off,one.pop = T)
  #head(off)

  #subsetting for r&d
  my.off <- sample(x = off$Sample.name,size = noff,replace = F)
 # my.off <- sample(x = off$Sample.name,size = 10,replace = F)
  off1 <- off[off$Sample.name %in% my.off,]
  #head(off1)

  #preping for strataG conversion
  off2 <- data.frame(loc = rep("loc1",times=nrow(off1)),off1,stringsAsFactors = F)
  off2 <- off2[,c(2,1,3:ncol(off2))]
  colnames(off2) <- gsub(pattern = "\\.1\\.",replacement = "_",x = colnames(off2))
  colnames(off2) <- gsub(pattern = "\\.",replacement = "_",x = colnames(off2))
  #head(off2)

  #converting to strataG format
  tmp4 <- gt.summary(population = off1,one.pop = T,pop.sum = T)

  ###################################
  ### Now calculating relatedness ###
  ###################################
  outfile <- suppressMessages(coancestry(off1,trioml = 1,wang = 0,lynchli = 0,lynchrd = 0,ritland = 0,quellergt = 0,dyadml = 0,working.directory = getwd(),output.file = FALSE))
  tmp4$related <- mean(outfile$relatedness$trioml)
  tmp4$sim <- i + 0
  #tmp4

  #saving it to out
  out  <- rbind(out,tmp4)
}

#writing to file
write.table(x = out,file = "output/randomized.hunt.creek.offspring.0.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)































print("hello")
