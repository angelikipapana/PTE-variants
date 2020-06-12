##########################################################################
# generating the optimal subset of conditioning variables based on the 
# random forest method for a multivariate variable in K variables (xM)
# For each variable, we determine the most 'important' determinants from
# the ensemble of the remaining variables
# - x : time series data, size n x K
# Code written by Angeliki Papana & Ariadni Papana-Dagiasis
#########################################################################

#libraries
#install.packages("randomForest")
#install.packages("randomForestSRC")
#install.packages("ggrandomForestSRC")
#install.packages("dplyr")
library(randomForest)
library(randomForestSRC)
library(dplyr)

rm(list=ls())


optimalSubsetGenerator <- function(x) {
  
  # Setting column names
 
  column_size  = dim(x)[2]
  column_names = paste("var", 1:column_size, sep="")
  colnames(x)  = c(column_names)
  
  
  #run random forest in regression mode------------------------
  #out1 = matrix(0L, nrow=column_size, ncol=column_size)
  out1 = matrix(0L, nrow=column_size*(column_size-1)/2, ncol=column_size)
  
  #out2 = matrix(0L, nrow=column_size, ncol=column_size+1)
  out2 = matrix(0L, nrow=column_size*(column_size-1)/2, ncol=1+column_size)
  
  counter = 0
  
  K = column_size
  
  K1 = K-1
  for (i in 1:K1) {
      
    i1 = i+1
    for (j in i1:K) {
        
      variabl1 = paste("var",i,sep="")
      variabl2 = paste("var",j,sep="")

      rf.full = rfsrc(as.formula(paste("Multivar(",variabl1,",",variabl2,") ~ .")), x, ntree=1000, importance=F, nodedepth = 5, nodesize = 40, nsplit = 10)
        
      print(rf.full)
      vs=var.select(rf.full, method=c("md"), conservative="high")
      
      vs$topvars
      mydepth=vs$varselect[,1]
      order(mydepth)
      
      my.vars=as.numeric(sapply(strsplit( vs$topvars , "var"), "[[",2))
      
      #create and save
      mydepth.out1= c(i, j, my.vars)
      mydepth.out2=c(i, j, vs$varselect[,1], vs$md.obj$threshold)
      
      counter = counter + 1
      out1[counter,1:length(mydepth.out1)]=mydepth.out1
      out2[counter,1:length(mydepth.out2)]=mydepth.out2
    
    }
  }
  list(out1, out2)
}
    
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("2 argument must be supplied. Data file and output file.", call.=FALSE)
}

datafile = args[1]

# todo remove
#datafile = "C:/angela/Publications/2020IJBC/Application/IntMarketK16/returns.dat"


#read data   
x=read.table(datafile, header=F)


output = optimalSubsetGenerator(x)


# table with ranks
filename1 = args[2]
write.table((output[1]), filename1, row.names=F, col.names=F, quote=F)

# table with depths
filename2 = args[3]
write.table((output[2]), filename2, row.names=F, col.names=F, quote=F)
