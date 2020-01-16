#originally created on: October 18th, 2017
#by: Nick Sard

#ABOUT: This script was written to great a set of randomly assembled location matrices representing the number of offspring
#assigned to a given type of location mate pair. 

#loading necessary libraries
library(tidyverse)
library(gridExtra)
library(grid)

#setting working directory
setwd("~/Desktop/Cloned_repositories/HuntCreek/Randomization/")

#listing files
list.files("Input/")

#reading gts for adults and offspring
df <- read.table(file = "Input/Hunt_Creek_all_genotypes_all_loci.txt",header = T,sep = "\t",stringsAsFactors = F)
tmp <- read.table(file = "Input/obs.values.run.2.txt",header = T, sep = "\t",stringsAsFactors = F)
head(df)
tmp
#getting a feel for the data
table(df$type,df$Sex)
table(df$loc)

#making three letter Abbreviations for each locatin
df$loc[df$loc == "EB_AuSable"] <- "EAS"
df$loc[df$loc == "Gilchrist_Creek"] <- "GCC"
df$loc[df$loc == "Houghton_Creek"] <- "HOC"
df$loc[df$loc == "Hunt_Creek_Res"] <- "HTC"
df$loc[df$loc == "Hunt_Creek_TrP"] <- "HTC"
table(df$loc)

#separating offspring and parents from df
par <- df[df$loc != "Hunt_Creek_Off",]
off <- df[df$loc == "Hunt_Creek_Off",]

#I need to create 1000 pedigrees where I randomly assign an offspring to two locations, one for dads and one for Moms

#getting vector of dads and moms locations - Including adults without a sex assigned in each
moms <- par[par$Sex != "M",]
dads <- par[par$Sex != "F",]
loc.mom <- moms$loc
loc.dad <- dads$loc
locs <- unique(par$loc)
out <-  NULL
out1 <- NULL
for(i in 1:1000){
  off.ped <- off[,c(2,3)]
  off.ped$mom <- sample(x = loc.mom,size = nrow(off.ped),replace = T)
  off.ped$dad <- sample(x = loc.dad,size = nrow(off.ped),replace = T)
  off.ped$type <- NULL
  off.ped$mp <- paste0(off.ped$dad,off.ped$mom)
  head(off.ped)
  
  #getting houses of each type
  off.ped.counts <- as.data.frame(table(off.ped$mp),stringsAsFactors = F)
  off.ped.counts
  
  #making a table of all possible assignments
  assignments <- as.data.frame(expand.grid(locs,locs,stringsAsFactors = F))
  assignments$mate.pairs <- paste0(assignments$Var1,assignments$Var2)
  assignments
  
  #merging the two together
  assignments <- merge(x = assignments,y = off.ped.counts,by.x = "mate.pairs",by.y = "Var1",all.x = T)
  assignments$Freq[is.na(assignments$Freq)] <-  0
  assignments
  
  #calculate chi-squared values for assignments
  mat <- assignments %>% 
    select(-mate.pairs) %>% 
    spread(key = Var2,value = Freq)
  row.names(mat) <- mat$Var1
  mat$Var1 <- NULL
  mat
  x <- chisq.test(mat)

  #making chi-squared table for each run
  x1 <- (x$observed)
  x1 <- as.data.frame(x1)
  x1$Var2<- row.names(x1)
  x1 <- x1 %>%gather(key = Var1,value = Val,-Var2)
  names(x1)[3] <- "obs"
  x1
  x2 <- (x$expected)
  x2 <- as.data.frame(x2)
  x2$Var2<- row.names(x2)
  x2 <- x2 %>%gather(key = Var1,value = Val,-Var2)
  names(x2)[3] <- "exp"
  head(x2)
  x1$exp <- x2$exp
  x1$diff <- x1$obs-x1$exp 
  x1$chi <- x1$diff^2/x1$exp
  out1 <- rbind(out1,x1)
  #saving the total chi-squared value for each run
  x <- as.numeric(x$statistic)
  out <- c(out,x)

}

#histogram of chi-squared values
hist(out)

#histograms of chi-squared values by location
head(out1)

head(tmp)
tmp
tmp$mp <- paste0(tmp$Female,tmp$Male)
out1$mp <- paste0(out1$Var2,out1$Var1)
head(out1)
tiff(filename = "Output/expected_histograms_by_location.tiff",width = 14,height = 11,units = "in",res = 300)
ggplot(out1, aes(x=diff))+
  geom_histogram()+
  facet_grid(Var1~Var2)
dev.off()
#calculate p-values
out[out > 65.5]/length(out)
tmp$pvalue <- 0
i <- 2
i <- NULL
for(i in 1:nrow(tmp)){
  out2 <- out1$chi[out1$mp == tmp$mp[i]]
  if(length(out2[out2 > tmp$chi_sq[i]]) != 0){
    tmp$pvalue[i] <- length(out2[out2 > tmp$chi_sq[i]])/length(out2)
  }

}
tmp$pval2 <- p.adjust(p = tmp$pvalue,method = "BH")
tmp

tmp$sig <- ifelse(tmp$pval2 < 0.05, "*","")
tmp$obs2 <- paste0(tmp$Obs,tmp$sig)
tmp$sign <- ifelse(tmp$Diff > 0, "+", "-")
tmp$sign[tmp$sig != "*"] <- NA
tmp
#ggplot
tiff(filename = "Output/expected_with_diff_heatmap.tiff",width = 14,height = 11,units = "in",res = 300)
ggplot(tmp,aes(x = Male, y = Female, fill=Diff, label = obs2))+
  geom_bin2d(color = "black")+
  geom_text(size = 5)+
  scale_fill_gradient2(high = "dark grey", low = "dark grey",mid = "white",midpoint = 0)+
  theme_classic()+
  labs(x = "Male Location", y = "Female Location",fill = "Obs.-Exp.")+
  theme(text=element_text(size=10,  family="serif"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(margin = margin(r = 0)),
        axis.text = element_text(size=12,colour = "black"),
        legend.title= element_text(size=12),
        legend.text = element_text(size =10))
dev.off()
