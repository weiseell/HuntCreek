#originally created on: December 7th, 2017
#by: Nick Sard

#ABOUT: This script was written to simulate offspring genotypes created by making from disparate (sig Fst diff) populations

#loading necessary libraries
library(tidyverse)

#setting working directory
setwd("Power_Analysis/")

#loading my own functions
source("source.scripts/power_functions.R")

#reading in the parental genotypes
df <- read.table("Input/Hunt_Creek_all_genotypes_all_loci.txt",sep = "\t",header = T,stringsAsFactors = F)
head(df)

noff <- df %>% 
  filter(type == "Offspring") %>% 
  nrow()
df <- df %>% 
  filter(type == "Adults") %>% 
  select(-type) %>% 
  filter(Sex != "UK")
head(df)
#using a for loop to prep the data for genepop format
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
head(df)

vector1 <- df$Sex
df <- df %>% 
  select(-Sex)

df <- six.char(population = df,missing.data = 0,one.pop = F)
head(df)
df$sex <- vector1
df <- df[,c(1,2,ncol(df),3:(ncol(df)-1))]

df$loc[df$loc == "Hunt_Creek_Res"] <- "aHCR"
df$loc[df$loc == "Houghton_Creek"] <- "aHOC"
df$loc[df$loc == "EB_AuSable"] <- "aEAS"
df$loc[df$loc == "Gilchrist_Creek"] <- "aGLC"
df$loc[df$loc == "Hunt_Creek_TrP"] <- "aHCR"
table(df$loc)
names(df)[1] <- "Pops"
names(df)[2] <- "Sample.name"

df$Sample.name <- paste(df$Pops,df$Sample.name,sep="_")
head(df)

#getting the genotypes in to six character states for moms and dads
moms <- df[df$sex == "F",c(-1,-3),]
head(moms)

#dads
dads <- df[df$sex == "M",c(-1,-3),]
head(dads)

#getting counts of each to make the breeding matrix
n.mom <- nrow(moms)
n.dad <- nrow(dads)
n.par <- n.mom*n.dad

#making the matring matrix
prob_vals = c(0.90,0.10)
max.iterations = 1000
i <- 1
while(i < max.iterations){
  print(i)

  #making a breeding a matrix and filling with informaiton
  mat.str <- matrix(data = rpois(n =  n.mom*n.dad,lambda = 4),nrow = n.dad,ncol = n.mom)
  mat.str
  
  #using a small function to determine with are kept or removed
  rs.keep <- function(){
    no <- prob_vals[1]
    yes <- prob_vals[2]
    return(sample(x = c("R","K"),size = 1,replace = T,prob = c(no,yes)))
  }
  my.keeps <- replicate(n.mom*n.dad,expr= rs.keep())
  my.keeps
  
  #making a new mat.str of the name size and filling with that information
  mat.str1 <- matrix(data = my.keeps,nrow = n.dad,ncol = n.mom)
  mat.str[which(mat.str1=="R")] <- 0
  test.dad <- which(rowSums(mat.str)==0)
  test.mom <- which(colSums(mat.str)==0)
  test.dad_L <- length(test.dad)
  test.mom_L <- length(test.mom)
  
  if(test.dad_L != 0 | test.mom_L != 0 ) {
    
    
    #if I get to the max interations I am going to force it to work
    if(i == max.iterations-1){
      
      #mothers first        
      if(test.mom_L != 0){
        k <- 1
        k <- NULL
        for(k in 1:length(test.mom)){
          
          #randomly picking a mom to mate with
          my.dad.pick <- sample(x = 1:n.dad,size = 1,replace = F)
          
          #making sure I get a non-zero RS value here
          my.mp.rs <- 0
          while(my.mp.rs < 1){my.mp.rs <- rpois(n = 1,lambda = 4)}
          mat.str[my.dad.pick,test.mom[k]] <- my.mp.rs
        } #end of mom loop
      }
      
      
      #now dads
      if(test.dad_L != 0){
        k <- 1
        k <- NULL
        for(k in 1:length(test.dad)){
          
          #randomly picking a mom to mate with
          my.mom.pick <- sample(x = 1:n.mom,size = 1,replace = F)
          
          #making sure I get a non-zero RS value here
          my.mp.rs <- 0
          while(my.mp.rs < 1){my.mp.rs <- rpois(n = 1,lambda = 4)}
          mat.str[test.dad[k],my.mom.pick] <- my.mp.rs
          mat.str[test.dad[k],]
        } #end of dad loop
      }
      
      #checking to make sure these things are true now
      test.dad <- which(rowSums(mat.str)==0)
      test.mom <- which(colSums(mat.str)==0)
      test.dad_L <- length(test.dad)
      test.mom_L <- length(test.mom)
      if(test.dad_L != 0 | test.mom_L != 0 ) {
        stop("Counld not fix the mat.str issue",j)
      } else {
        i <- max.iterations+199 
      }
    } #end of hard-fix statement
    i <- i + 1
    
  } else {
    i <- max.iterations+100
  }
}
head(mat.str)

#making that into a df and assigning colnames to moms and making a new col for dads
tmp <- data.frame(mat.str)
colnames(tmp) <- sample(x = moms[,1],size = length(moms[,1]),replace=F)
tmp <- cbind(sample(x = dads[,1],size = length(dads[,1]),replace=F),tmp)
names(tmp)[1] <-"dads"
head(tmp)

#making this into long form and removing all mate pairs that were unsuccessful
tmp1 <- gather(data = tmp,key = "moms",value = "noff",-dads) %>% filter(noff != 0)
tmp1$dads <- as.character(tmp1$dads)
str(tmp1)
head(tmp1)

#for each mate pair I will produce offspring genotypes using the laws of independent segregation and assortment
off <- make.gts(tmp1 = tmp1,moms = moms,dads = dads)
head(off)

#need to randomly subset the dataset so that it represents the same number of offpsring we have
myoff <- sample(x = off$Sample.name,size = noff,replace = F)
off <- off[off$Sample.name %in% myoff,]

#converting back to three.char for each set of gts
moms <- three.char(tmp = moms,one.pop = T)
dads <- three.char(tmp = dads,one.pop = T)
off <- three.char(tmp = off,one.pop = T)

#making a markers file
markers <- matrix(data = 0,nrow = 4,ncol = (ncol(off)-1)/2)
markers[1,] <- colnames(off)[seq(from =2,to=ncol(off)-1,by = 2)]
markers[2,] <- 2
markers[3,] <- 0.02
markers[4,] <- 0.001
markers

#creating colony dat file
colony.dat.create.R(moms = moms,dads = dads,kids = off,markers = markers,update.alfs = 0,spp.type = 2,
                    inbreeding = 0,ploidy = 0,fem.gamy = 1,mal.gamy = 1,clone = 0,sib.scale = 0,sib.prior = 0,
                    known.alfs = 0,run.number = 1,run.length = 2,monitor = 0,windows.version = 0,
                    full.likelihood = 1,likelihood.precision = 2,prob.mom = 1.0,prob.dad = 1.0)

#fin!