#originally created on: December 7th, 2017
#by: Nick Sard

#ABOUT: This script was written to 

#loading necessary libraries
library(tidyverse)

#setting working directory and reading in data
setwd("H:/Scribner_Lab/Hunt Creek Analysis/New Analysis/halfsib_sims")

#loading my own functions
source("H:/Scribner_Lab/Hunt Creek Analysis/New Analysis/actual_data_analysis/source.scripts/unique.pairs.R")

#getting the list of files
my.files <- list.files(path = "Input/",pattern = "Best")

#reading in location informatoni
locs <- read.table(file = "Input/Hunt_Creek_all_genotypes_all_loci.txt",header = T,sep = "\t",stringsAsFactors = F)
head(locs)

#getting rid of the offp
locs <- filter(locs,  type == "Adults")

#fixing names
table(locs$loc)
locs$loc[locs$loc == "EB_AuSable"] <- "ASR"
locs$loc[locs$loc == "Houghton_Creek"] <- "HGC"
locs$loc[locs$loc == "Gilchrist_Creek"] <- "GCC"
locs$loc[grepl(pattern = "Hunt_Creek",x = locs$loc)] <- "HTC"
table(locs$loc)

#getting just the columns I want
locs <- locs %>% select(loc,id)
head(locs)

i <- 5
i <- NULL
tmp <- NULL
conf <- NULL
uks <- NULL
for(i in 1:length(my.files)){
  print(i)
  #reading each pedigree
  df <- read.table(file = paste0("Input/",my.files[i]),header = T,sep = "\t",stringsAsFactors = F,strip.white = T)
  colnames(df) <- c("off","dad","mom")
  df <- select(df, off,mom,dad)
  head(df)
  #extracting true parent information from offspring names, while preserving teh actual offspring name
  df$off1 <- df$off
  df <- separate(data = df,col = off1,into = c("off1","tmom","tdad"),sep = "@")
  df$mom <- gsub(pattern = "_",replacement = "-",x = df$mom)
  df$dad <- gsub(pattern = "_",replacement = "-",x = df$dad)
  df$tmom <- gsub(pattern = "_",replacement = "-",x = df$tmom)
  df$tdad <- gsub(pattern = "_",replacement = "-",x = df$tdad)
  head(df)
  
  #getting the location information for both inferred and true moms and dads
  df <- merge(x = df,y = locs,by.x = "mom",by.y = "id",all.x = T)
  names(df)[ncol(df)] <- "mom_loc"
  df <- merge(x = df,y = locs,by.x = "dad",by.y = "id",all.x = T)
  names(df)[ncol(df)] <- "dad_loc"
  df <- merge(x = df,y = locs,by.x = "tmom",by.y = "id",all.x = T)
  names(df)[ncol(df)] <- "tmom_loc"
  df <- merge(x = df,y = locs,by.x = "tdad",by.y = "id",all.x = T)
  names(df)[ncol(df)] <- "tdad_loc"
  
  #doing teh logical test and accounting for UK parents
  df$mom_test <- ifelse(test = df$mom_loc == df$tmom_loc,0,1)
  df$dad_test <- ifelse(test = df$dad_loc == df$tdad_loc, 0,1)
  df$mom_test[df$mom == "UK"] <- NA
  df$dad_test[df$dad == "UK"] <- NA
  head(df)
  
  if(nrow( df[is.na(df$mom_test) | is.na(df$dad_test),7:10]) > 0){
    tmp3 <- df[is.na(df$mom_test) | is.na(df$dad_test),7:10]
    tmp3$file <- my.files[i] 
  }


  #getting summary stats and recording them in tmp1
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
  
  #saving them for later
  tmp <- rbind(tmp,tmp1)
  
  df2 <- data.frame(tmom_loc = rep("A",12),
             mom_loc = rep(letters[1:4],times=3),stringsAsFactors = F)
  as.data.frame(table(df2$tmom_loc,df2$mom_loc))
  
  if(nrow(tmp3)>0){
  uks <- rbind(uks,tmp3)
  }
  #looking at the confusion matrix now
  moms <- as.data.frame(table(df$tmom_loc,df$mom_loc))
  moms$type <- "mom"
  dads <- as.data.frame(table(df$tdad_loc,df$dad_loc))
  dads$type <- "dad"
  conf1 <- rbind(moms,dads)
  ggplot(conf1,aes(x = Var1,y = Var2,fill= Freq, label = Freq))+
    facet_wrap(~type)+
    geom_tile(color="black",stat = "identity")+
    geom_text(size=14)+
    scale_fill_gradient2(high = "grey") +
    labs(x = "True parent Population",y="Inferred parent population")
}
mean(tmp$mom.p.disc)
mean(tmp$dad.p.disc)

head(uks)
uks$type <- "moms"
uks$type[is.na(uks$dad_loc)] <- "dad"
table(uks$type)
moms1 <- uks %>%
  filter(type == "moms")%>%
  group_by(type,tmom_loc) %>%
  count()
moms1
moms1$n/ sum(moms1$n)
dads1 <- uks %>%
  filter(type != "moms")%>%
  group_by(type,tmom_loc) %>%
  count()
dads1$n/ sum(dads1$n)
#write to file
write.table(tmp, file = "Output/accuracy_loc.summary.txt",append = F, quote = F, sep = "\t",row.names = F, col.names = T)

#fin!


