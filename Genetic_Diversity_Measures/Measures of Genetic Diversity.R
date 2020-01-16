#Basic Population Analysis for Hunt Creek Data

#Created by: Ellie Weise
#Originally Created on: 10/16/17

#Objectives:
#1. Summarize allele frequencies, expected and observed heterozygosities, and the number of alleles per locus
#2. Break summary down by location
#3. Display summaries in tables
#4. Display summaries in ggplot
#5: Summarize allele frequencies by location by strata in a graph

#loading libraries
library(strataG)
library(adegenet)
library(ape)
library(tidyverse)
library(hierfstat)

#set wd
setwd("~/Desktop/Cloned_repositories/HuntCreek/Genetic_Diversity_Measures/")

#read in functions/source scripts:
source("source/unique.pairs.R")

#read in data
tmp <- read.csv("Input/Hunt_Creek_all_genotypes.csv",stringsAsFactors = FALSE)

#clean up data and convert to gtypes format####
head(tmp)
#abbreviating all locations
tmp$loc[tmp$loc == "EB_AuSable"] <- "EAS"
tmp$loc[tmp$loc == "Houghton_Creek"] <- "HOC"
tmp$loc[tmp$loc == "Gilchrist_Creek"] <- "GCC"
tmp$loc[tmp$loc == "Hunt_Creek_Res"] <- "HTC"
tmp$loc[tmp$loc == "Hunt_Creek_Off"] <- "HCO"
tmp$loc[tmp$loc == "Hunt_Creek_TrP"] <- "HTC"

#creating a file to put into gf2gtypes
pars <- tmp[-2]
pars$Sex <- NULL
head(pars)
#running df2gtypes (strataG)
msats <- df2gtypes(pars, ploidy = 2, id.col = 2,strata.col = 1,loc.col = 3)
msats

#making a data frame to hold all stats calculated
stats <- NULL

##individual ANOVA####
#calculating number of alleles and heterozygosity for each individual
tmp_geno <- tmp %>% 
  select(Omy301.1:One13.2)
#get the number of unique genotypes for each allele
tmp_geno$num_all <- apply(tmp_geno,1,function(x)length(unique(x)))
#calculate observed heterozygosity for each individual
tmp_geno <- tmp_geno %>% 
  mutate(Ho = num_all/16)
#combine stats and genotypes back with individual ids
tmp_id <- tmp %>% 
  select(loc:Sex)
tmp <- cbind(tmp_id,tmp_geno)

#subsetting tmp by location
tmp_stat <- tmp %>% select(loc,num_all,Ho)
EAS <- subset(tmp_stat,tmp_stat$loc == "EAS")
HOC <- subset(tmp_stat,tmp_stat$loc == "HOC")
GCC <- subset(tmp_stat,tmp_stat$loc == "GCC")
HTC <- subset(tmp_stat,tmp_stat$loc == "HTC")

all_allele <- tmp %>% 
  select(loc,num_all)
#subsetting for only adult populations
all_allele <- subset(all_allele,all_allele$loc != "HCO")
all_Ho <- tmp %>% 
  select(loc,Ho)
all_Ho <- subset(all_Ho,all_Ho != "HCO")
#getting rid of rows with NAs
all_Ho <- na.omit(all_Ho)

#ANOVA over all populations
fm1 <- aov(num_all~loc,data=all_allele)
anova(fm1)
fm2 <- aov(Ho~loc,data=all_Ho)
anova(fm2)

#doing t-test between each source and Hunt Creek####
#random sample of hunt creek residents to equalize sample size
HTC_all_samp <- sample(x = HTC$num_all, size = 14, replace = F)
HTC_Ho_samp <- sample(x = HTC$Ho, size = 14, replace = F)

#t-test for number of alleles
t.test(HTC_all_samp,EAS$num_all)
t.test(HTC_all_samp,HOC$num_all)
t.test(HTC_all_samp,GCC$num_all)
#t-test for heterozygosity
t.test(HTC_Ho_samp,EAS$Ho)
t.test(HTC_Ho_samp,HOC$Ho)
t.test(HTC_Ho_samp,GCC$Ho)

#calculating number of indiv per location(strata)####
num <- as.data.frame(table(strata(msats)))
colnames(num) <- c("loc","N")
rbind(stats,num)

#calculating obs and exp heterozygosity####
summ <- as.data.frame(summarizeLoci(msats,by.strata = TRUE))
summ$locus <- row.names(summ)
summ <- summ %>% 
  gather(key = "info",value = "val",-locus) %>%
  mutate(info = sub(pattern = "\\.",replacement = "_",x = info)) %>%
  separate(info, into = c("loc","stat"),sep = "_")

#finishing table of summary stats ####
#rename stat headings
summ$stat2 <- "NONE"
summ$stat2[summ$stat == "allelic.richness"] <- "Ar"
summ$stat2[summ$stat == "exptd.heterozygosity"] <- "He"
summ$stat2[summ$stat == "FIS"] <- "Fis"
summ$stat2[summ$stat == "num.alleles"] <- "A"
summ$stat2[summ$stat == "num.genotyped"] <- "Ngt"
summ$stat2[summ$stat == "obsvd.heterozygosity"] <- "Ho"
summ$stat2[summ$stat == "prop.genotyped"] <- "Pgt"
summ$stat2[summ$stat == "prop.unique.alleles"] <- "Puniq"
table(summ$stat2)

#Calculate FIS by location####
pop.stats <- summ %>% 
  select(-stat) %>%
  spread(key = stat2,value = val) %>%
  mutate(FIS = (1-(Ho/He))) %>%
  gather(key = "stat",value = "val",A:FIS) %>%
  group_by(loc,stat) %>%
  summarize(means = mean(val),
            sds = sd(val),
            ses = sds/sqrt(n())) %>%
  select(loc:means) %>%
  spread(key=stat,value=means)
pop.stats

#calculating Fis distrbutions overall loci and pops####
msats1 <- gtypes2genind(msats)
str(msats1)

x <- basic.stats(msats1)
x
y <- as.data.frame(x$Fis)
str(y)
y
pop.stats
y1 <- as.data.frame.list(colMeans(y))

y1 <- y1 %>% 
  gather()
colnames(x = y1) <- c("loc","Fis")
pop.stats <- merge(pop.stats,y1)

#Correcting p-values for heterozygosity####
#nonparametric test for heterozygosity
tmp2 <- summ %>% 
  select(-stat) %>%
  spread(key = stat2,value = val) %>%
  mutate(FIS = (1-(Ho/He)))
head(tmp2)
tmp3 <- unique.pairs(ids = unique(tmp$loc))
tmp3$stat <- NA
tmp3$pval <- NA
i <- 1
i <- NULL
for (i in 1:nrow(tmp3)){
  x <- tmp2 %>% filter(loc == tmp3$Var1[i])
  y <- tmp2 %>% filter(loc == tmp3$Var2[i])
  my.stat <- wilcox.test(x = x$He,y = y$He,paired = T)
  str(my.stat)
  tmp3$stat[i] <- my.stat$statistic 
  tmp3$pval[i] <- my.stat$p.value
  
}
tmp3

#test corrections
tmp3$pval2 <- p.adjust(p = tmp3$pval, method = "BH")
tmp3$type <- "He"
out <- tmp3

#correcting p-values for number of alleles####
#parametric test for allele number
tmp3 <- unique.pairs(ids = unique(tmp$loc))
tmp3$stat <- NA
tmp3$pval <- NA
i <- 1
i <- NULL
for (i in 1:nrow(tmp3)){
  x <- tmp2 %>% filter(loc == tmp3$Var1[i])
  y <- tmp2 %>% filter(loc == tmp3$Var2[i])
  my.stat <- wilcox.test(x = x$A,y = y$A,paired = T)
  str(my.stat)
  tmp3$stat[i] <- my.stat$statistic 
  tmp3$pval[i] <- my.stat$p.value
  
}
tmp3

#test corrections
tmp3$pval2 <- p.adjust(p = tmp3$pval, method = "BH")
tmp3$type <- "A"
tmp3

out <- rbind(out,tmp3)
out

#Fst calculation using popStructTest from strataG -- NOT FAST ~ 4 min####
fst = popStructTest(msats, nrep = 100, type = "both", stats = "fst", quietly = TRUE)
fst1 <- as.data.frame(fst$pairwise$result)
fst1 <- fst1 %>% 
  select(strata.1,strata.2,Fst,Fst.p.val) %>% 
  rename(locA = strata.1,locB = strata.2,p.val = Fst.p.val)
#making all values less than 0 just 0
fst1$value2 <- fst1$Fst
fst1$value2[fst1$Fst <0]<-0

#test corrections
fst1$pval2 <- p.adjust(p = fst1$p.val, method = "fdr")
fst1$type <- "Fst"
fst2 <- fst1 %>% select(locA,locB,Fst,p.val,pval2,type)
colnames(fst2) <- colnames(out)
out <- rbind(out,fst2)
out
#Fst Plot####
#the plot
tiff(filename = "Output/fst.plot.tiff",width = 14,height = 11,units = "in",res = 300)
ggplot(fst1, aes(x = locA, y=locB, label=round(value2,3)))+
  geom_tile(color="black",fill = "white",size = 1)+
  geom_text(size=10)+
  scale_y_discrete(position = "left",
                   labels = c("Gilchrist\n Creek","Hunt Creek\n Offspring","Houghton\n Creek","Hunt \n Creek"))+
  scale_x_discrete(position = "top",
                   labels = c("East Au\n Sable","Gilchrist\n Creek","Hunt Creek\n Offspring","Houghton\n Creek"))+
  theme_bw()+
  labs(x = "",y="", fill="Fst")+
  theme(text=element_text(size=20,  family="serif"),
        axis.text = element_text(size=24,colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
dev.off()

#write to file
write.table(x = pop.stats,file = "Output/pop.gen.statistics.txt",append = F,sep = "\t",row.names = F,quote = F,col.names = T)
write.table(x = out,file = "Output/wilcox.He.A.txt",append = F,sep = "\t",row.names = F,quote = F,col.names = T)
