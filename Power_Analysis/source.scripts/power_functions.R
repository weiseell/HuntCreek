#################
### six.char() ###

#Description
# This fuction takes genotypes in three character form and puts them into six character form

#Input Parameters:

# population - a table with sample.name in column one, and loci after that. All loci must be in three character form. Examples 000, 111, 999
six.char <- function(population, missing.data = "0", one.pop=T, my.sep = "") {
  
  missing.data <- as.numeric(missing.data)
  missing.data1 <- paste0(missing.data,missing.data,missing.data)
  missing.data1 <- paste0(missing.data1,my.sep,missing.data1)
  missing.data <- paste0(missing.data,my.sep,missing.data)
  cnames <- colnames(population)
  
  if(one.pop==T){  
    pop.names <- population[,1]
    population=population[,-1]  
  } else{
    pop.names <- population[,1]
    id.names <- population[,2]
    population=population[,-c(1,2)]
  }
  
  L<-1:ncol(population)
  locus.locs <- L[seq(1,ncol(population),2)]
  locus.names <- colnames(population)
  
  OUT <- NULL
  for (x in locus.locs) {
    output=paste(population[,x],population[,x+1],sep=my.sep)
    OUT <- cbind(OUT,output)
  }
  
  locus.names <- locus.names[seq(1,length(locus.names),2)]
  
  if(one.pop==T){
    OUT <- data.frame(pop.names,OUT)
    colnames(OUT) <- c(cnames[1],locus.names)
  } else {
    OUT <- data.frame(pop.names,id.names,OUT)
    colnames(OUT) <- c(cnames[1:2],locus.names)
  }
  
  OUT <- data.frame(lapply(OUT, as.character), stringsAsFactors=FALSE)
  OUT[OUT == missing.data] <- missing.data1
  return(OUT)
}
#################
### three.char() ###

#Description
# This fuction takes genotypes in six character form and puts them into three character form
#Input Parameters:

# population - a table with sample.name in column one, and loci after that. All loci must be in six character form. Examples 000000, 111111, 999999

three.char <- function(tmp, one.pop=T){
  
  if(one.pop==T){
    
    #saving sample names into a seperate file
    OUT <- data.frame(tmp[,1])
    colnames(OUT) <- colnames(tmp)[1]
    
    #just grabbing all the genotypes with no sample names
    gts <- tmp[,-1]  
  } else{
    
    #saving sample names into a seperate file
    OUT <- data.frame(tmp[,1:2])
    colnames(OUT) <- colnames(tmp)[1:2]
    
    #just grabbing all the genotypes with no sample names
    gts <- tmp[,-1:-2]
  }
  
  
  #counting the number of loci and saving their names
  L <- ncol(gts)
  Locus.names <- colnames(gts)
  
  #using a small for loop to take each column and split it into two columns and cbinding that to OUT, which originally just had the sample names in it
  for(i in 1:ncol(gts)){
    locus <- data.frame(gts[,i])
    colnames(locus) <- Locus.names[i]
    
    locus$x1 <- substr(locus[,1],1,3)
    locus$x2 <- substr(locus[,1],4,6)
    locus[,1] <- NULL
    
    colnames(locus) <- c(paste(Locus.names[i],".1", sep=""),paste(Locus.names[i],".2",sep=""))
    OUT <- cbind(OUT,locus)
  }
  OUT[,1] <- as.character(OUT[,1])
  return(OUT)
}

#############
###make.gts and ind.seg###

#Description
#produce offspring genotypes using the laws of independent segregation and assortment
#ind.seg is a sampling function that is used in make.gts
ind.seg <- function(mgt,dgt){
  mom_a1 <- substr(mgt,1,3)
  mom_a2 <- substr(mgt,4,6)
  dad_a1 <- substr(dgt,1,3)
  dad_a2 <- substr(dgt,4,6) 
  mom_a <- sample(x = c(mom_a1,mom_a2),size = 1)
  dad_a <- sample(x = c(dad_a1,dad_a2),size = 1)
  off_gt <- paste0(sort(c(mom_a,dad_a)),collapse = "")
  return(off_gt)
}

make.gts <- function(tmp1, moms, dads){
  i <- 1
  i <- NULL
  out <- NULL
  for(i in 1:nrow(tmp1)){
    #identifying the female and male in the mate pair to make the off genotypes with 
    momgt <- filter(moms,Sample.name == tmp1$moms[i])
    dadgt <- filter(dads, Sample.name == tmp1$dads[i])
    #determining the number of offspring they are going to produce
    myoff <- tmp1$noff[i]
    #for each offspring to be simulated per mate pair
    j <- NULL
    for(j in 1:myoff){
      #making their genotypes
      k <- NULL
      offgt <- NULL
      for(k in 2:ncol(momgt)){
        #print(k)
        offgt1 <- ind.seg(mgt = momgt[1,k],dgt = dadgt[1,k])
        offgt <- cbind(offgt,offgt1)
      } #end of gt assembly for (k) loop
      #making the off name and getting the columns right
      off_name <- paste(paste("off",j,sep="_"),momgt[1,1],dadgt[1,1],sep="@")
      offgt <- data.frame(cbind(off_name,offgt),stringsAsFactors = F)
      colnames(offgt) <- colnames(momgt)
      out <- rbind(out,offgt)
    } #end of noff for (j) loop
  }  #enf of matepair for (i) loop
  return(out)
} #end of function

###########
###colony.dat.create###

#Description
#creates a colony input file. requires all parameters to be input into the function as well as several input tables.
colony.dat.create.R <- function(moms, dads, kids, markers, update.alfs = 0,
                                spp.type = 2, inbreeding = 0, ploidy = 0, fem.gamy = 1, mal.gamy = 1,
                                clone = 0, sib.scale = 1, sib.prior = 0, known.alfs = 0, run.number = 1,
                                run.length = 2, monitor = 0, windows.version = 0, full.likelihood = 1,
                                likelihood.precision = 2, prob.mom = 1.0, prob.dad = 1.0){
  
  #getting current working directory and fixing slashes for running on linux
  my.wd <- "'Input.data'"
  my.wd2 <- "'Output.data'" 
  
  #getting the number of kids, moms, and dads
  noff <- nrow(kids); nmoms <- nrow(moms); ndads <- nrow(dads)
  noff <- paste(noff,"! Number of offspring in the sample",sep = "\t")
  
  #getting the number of loci
  nloci <- ncol(markers)
  nloci <- paste(nloci, "! Number of loci",sep = "\t") 
  
  #setting a random number seed
  random.seed <- round(runif(n = 1,min = 1,max = 9999))
  random.seed <- paste(random.seed, "! Seed for random number generator",sep = "\t")
  
  #adding comments to input parameters
  update.alfs <- paste(update.alfs,"! Not upate/update allele frequency",sep = "\t")
  spp.type <- paste(spp.type, "! 2/1=Dioecious/Monoecious",sep = "\t")
  inbreeding <- paste(inbreeding,"! 0/1=No inbreeding/inbreeding", sep = "\t")
  ploidy <- paste(ploidy, "! 0/1=Diploid species/HaploDiploid species",sep="\t")
  gamy <- paste(paste(mal.gamy,fem.gamy),"! 0/1=Polygamy/Monogamy for males & females",sep = "\t")
  clone <- paste(clone,"! 0/1=Clone inference =No/Yes",sep="\t")
  sib.scale <- paste(sib.scale,"! 0/1=Scale full sibship=No/Yes",sep = "\t")
  sib.prior <- paste(sib.prior,"! 0/1/2/3=No sibship prior/Weak sibship prior/Medium sibship prior/Strong sibship prior",sep = "\t")
  known.alfs <- paste(known.alfs,"! 0/1=Unknown/Known population allele frequency",sep = "\t")
  run.number <- paste(run.number,"! Number of runs",sep = "\t")
  run.length <- paste(run.length,"! Length of run",sep = "\t")
  monitor <- paste(monitor,"! 0/1=Monitor method by Iterate#/Time in second",sep = "\t")
  monitor.interval <- paste(100000,"! Monitor interval in Iterate# / in seconds",sep = "\t")
  windows.version <- paste(windows.version,"! Windows version",sep = "\t")
  full.likelihood <- paste(full.likelihood,"! Fulllikelihood",sep = "\t")
  likelihood.precision <- paste(likelihood.precision,"! 1/2/3=low/medium/high Precision for Fulllikelihood",sep = "\t")
  
  #collating info for parents
  if(prob.dad == 0 & prob.mom == 0){
    probs <- c("0.0","0.0")
  } else {
    probs <- c(prob.dad,prob.mom)
  }
  
  npars <- c(nrow(dads),nrow(moms))
  
  my.value <- paste0(0,"\n")
  
  #making the actual file
  cat(my.wd,my.wd2,noff,nloci,random.seed, update.alfs,spp.type,
      inbreeding,ploidy,gamy,clone,sib.scale,sib.prior,
      known.alfs,run.number,run.length,monitor,monitor.interval,
      windows.version,full.likelihood,likelihood.precision,file = "colony2.dat",sep = "\n",append = T)
  cat("\n",file = "colony2.dat",append = T)
  write.table(x = markers,file = "colony2.dat",append = T,quote = F,sep = ",",row.names = F,col.names = F)
  cat("\n",file = "colony2.dat",append = T)
  write.table(x = kids,file = "colony2.dat",append = T,quote = F,sep = " ",row.names = F,col.names = F)
  cat("\n",file = "colony2.dat",append = T)
  cat(probs,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(npars,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  write.table(x = dads,file = "colony2.dat",append = T,quote = F,sep = " ",row.names = F,col.names = F)
  cat("\n",file = "colony2.dat",append = T)
  write.table(x = moms,file = "colony2.dat",append = T,quote = F,sep = " ",row.names = F,col.names = F)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
}

