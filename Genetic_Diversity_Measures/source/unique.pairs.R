#################
### unique.pairs
#################

#Description
#This fuction populations a data frame with all unique pairs of IDs (dyads)

#Input Parameters:
#ids <- a vector containing 2 or more unique id names

#first defining a function I will use later
unique.pairs <- function(ids){
  ids <- as.character(ids)
  x <- expand.grid(ids,ids)
  x <- x[ifelse(x$Var1==x$Var2,T,F)==F,]
  x$both <- paste(pmin(as.character(x$Var1), as.character(x$Var2)),
                  pmax(as.character(x$Var1), as.character(x$Var2)), sep="_")
  x <- x[duplicated(x$both)==F,]
  x$both <- NULL
  row.names(x) <- NULL
  x <- x[,c(2,1)]
  colnames(x) <- c("Var1","Var2")
  return(x)
} # end of unique.pairs