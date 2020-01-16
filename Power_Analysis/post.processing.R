#originally created on: December 7th, 2017
#by: Nick Sard

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
locs

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


i <- 2
i <- NULL
ind_accuracy <- NULL
loc_accuracy <- NULL
rs_accuracy <- NULL
for(i in 1:length(my.files)){
  
  #reading in the pedigree and making columns for the real parents
  df <- read.table(file = paste0("Input/",my.files[i]),header = T,sep = "\t",stringsAsFactors = F,strip.white = T)
  colnames(df) <- c("off","dad","mom")
  df <- select(df, off,mom,dad)
  df$off1 <- df$off
  df <- separate(data = df,col = off1,into = c("off1","tmom","tdad"),sep = "@")
  head(df)

  #accuracy by individual
  df$mom_test <- ifelse(test = df$mom == df$tmom,0,1)
  df$dad_test <- ifelse(test = df$dad == df$tdad, 0,1)
  df$mom_test[grepl(pattern = "mp", x = df$mom)] <- NA
  df$dad_test[grepl(pattern = "mp", x = df$dad)] <- NA
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
  tmp1
  
  #saving individual analysis for later
  ind_accuracy <- rbind(ind_accuracy,tmp1)
  
  #accurary by location
  #need for fix names for merge
  df$mom <- gsub(pattern = "_",replacement = "-",x = df$mom)
  df$dad <- gsub(pattern = "_",replacement = "-",x = df$dad)
  df$tmom <- gsub(pattern = "_",replacement = "-",x = df$tmom)
  df$tdad <- gsub(pattern = "_",replacement = "-",x = df$tdad)
  head(df)
  head(locs)
  #getting the location information for both inferred and true moms and dads
  df <- merge(x = df,y = locs,by.x = "mom",by.y = "id",all.x = T)
  names(df)[ncol(df)] <- "mom_loc"
  df <- merge(x = df,y = locs,by.x = "dad",by.y = "id",all.x = T)
  names(df)[ncol(df)] <- "dad_loc"
  df <- merge(x = df,y = locs,by.x = "tmom",by.y = "id",all.x = T)
  names(df)[ncol(df)] <- "tmom_loc"
  df <- merge(x = df,y = locs,by.x = "tdad",by.y = "id",all.x = T)
  names(df)[ncol(df)] <- "tdad_loc"
  
  #doing the logical test and accounting for UK parents
  df$mom_test <- ifelse(test = df$mom_loc == df$tmom_loc,0,1)
  df$dad_test <- ifelse(test = df$dad_loc == df$tdad_loc, 0,1)
  df$mom_test[grepl(pattern = "mp", x = df$mom)] <- NA
  df$dad_test[grepl(pattern = "mp", x = df$dad)] <- NA
  
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
  tmp1
  
  #saving loc analysis for later
  loc_accuracy <- rbind(loc_accuracy,tmp1)
  
  #need for fix names for merge
  df$tmom <- gsub(pattern = "_",replacement = "-",x = df$tmom)
  df$tdad <- gsub(pattern = "_",replacement = "-",x = df$tdad)
  
  #true parents rs estimate
  tmom <- unique(df$tmom)
  tdad <- unique(df$tdad)
  df_t <- data.frame(par = c(tmom,tdad),
                     sex = c(rep("F",times=length(tmom)),rep("M",time=length(tdad))),
                     stringsAsFactors = F)
  head(df_t)
  
  #getting location information
  head(locs)
  df_t <- merge(x = df_t,y = locs,by.x = "par",by.y = "id",all.x = T)
  head(df_t)
  
  j <- 1
  j <- NULL
  for(j in 1:nrow(df_t)){
    if(df_t$sex[j] == "F"){df_t$rs[j] <- nrow(df[df$tmom == df_t$par[j],])}
    if(df_t$sex[j] == "M"){df_t$rs[j] <- nrow(df[df$tdad == df_t$par[j],])}
  }
  #inferred parents rs estimate
  mom <- unique(df$mom)
  dad <- unique(df$dad)
  df_i <- data.frame(par = c(mom,dad),
                     sex = c(rep("F",times=length(mom)),rep("M",time=length(dad))),
                     stringsAsFactors = F)
  head(df_i)
  
  #getting location information
  head(locs)
  df_i <- merge(x = df_i,y = locs,by.x = "par",by.y = "id",all.x = T)
  head(df_i)
  
  j <- 1
  j <- NULL
  for(j in 1:nrow(df_i)){
    if(df_i$sex[j] == "F"){df_i$rs[j] <- nrow(df[df$mom == df_i$par[j],])}
    if(df_i$sex[j] == "M"){df_i$rs[j] <- nrow(df[df$dad == df_i$par[j],])}
  }
  df_i$type <- "Inferred"
  df_t$type <- "True"
  df1 <- rbind(df_t,df_i)
  head(df1)
  df1 <- df1 %>%
    filter(!is.na(loc)) %>%
    group_by(type,loc) %>%
    summarize(means = round(mean(rs),2)) %>%
    spread(key = "type",value = "means",fill = NA)
  df1$ratio <- round(df1$Inferred/df1$True,2)
  df1$run <- my.files[i]
  df1
  #saving for later
  rs_accuracy <- rbind(rs_accuracy,df1)
  
}
head(ind_accuracy)
head(loc_accuracy)
head(rs_accuracy)

ind_sum <- ind_accuracy %>% 
  gather(key = 'sex',value = 'disc',mom.p.disc:dad.p.disc) %>% 
  mutate(sex = ifelse(sex == "mom.p.disc","Mother","Father")) %>% 
  group_by(sex) %>% 
  summarise(means = mean(disc), 
            sds = sd(disc),
            counts = n(),
            ses = sds/sqrt(counts)) %>% 
  mutate(type = "Individual")

loc_sum <- loc_accuracy %>% 
  gather(key = 'sex',value = 'disc',mom.p.disc:dad.p.disc) %>% 
  mutate(sex = ifelse(sex == "mom.p.disc","Mother","Father")) %>% 
  group_by(sex) %>% 
  summarise(means = mean(disc), 
            sds = sd(disc),
            counts = n(),
            ses = sds/sqrt(counts)) %>% 
  mutate(type = "Location")
loc_sum

all_sum <- rbind(ind_sum,loc_sum)
all_sum
#making a figure of accuracy by individual and by location
tiff(filename = "Output/discordance_indiv_loc.tiff",width = 14,height = 11,units = "in",res = 300)
ggplot(all_sum, aes(x=type,y=means,fill=sex))+
  geom_bar(stat="identity",color="black",position="dodge")+
  geom_errorbar(aes(ymin=means-ses,ymax=means+ses),position = position_dodge(0.9),width = .2)+
  scale_fill_grey()+
  theme_bw()+
  labs(x="",y="Mean Prop. of Offspring Misidentified",fill="")+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14, color="black"),
        strip.text = element_text(size = 14),
        legend.position = "none")
dev.off()

#write to file
write.table(ind_accuracy, file = "Output/accuracy.summary.indiv.txt",append = F, quote = F, sep = "\t",row.names = F, col.names = T)
write.table(loc_accuracy, file = "Output/accuracy.summary.loc.txt",append = F, quote = F, sep = "\t",row.names = F, col.names = T)
write.table(rs_accuracy, file = "Output/accuracy.summary.rs.txt",append = F, quote = F, sep = "\t",row.names = F, col.names = T)
#fin!




