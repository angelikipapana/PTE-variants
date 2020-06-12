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
    out1 = matrix(0L, nrow=column_size, ncol=column_size)
    out2 = matrix(0L, nrow=column_size, ncol=column_size+1)
    
    i = 0
    for (column_name in colnames(x)) {
      i=i+1
      print(column_name)
      
      
      #variabl1 = paste("var",i,sep="")
      #variabl2 = paste("var",j,sep="")
      
      #rf.full = rfsrc(as.formula(paste("Multivar(",variabl1,",",variabl2,") ~ .")), x, ntree=1000, importance=F, nodedepth = 5, nodesize = 40, nsplit = 10)
      
      
      
      rf.full = rfsrc(as.formula(paste(column_name," ~ .")), x, ntree=1000, importance=F, nodedepth = 5, nodesize = 40, nsplit = 10)
      
      print(rf.full)
      
      #vs=var.select(rf.full, method=c("md"))
      vs=var.select(rf.full, method=c("md"), conservative="high")
      
      vs$topvars
      mydepth=vs$varselect[,1]
      order(mydepth)
      
      my.vars=as.numeric(sapply(strsplit( vs$topvars , "var"), "[[",2))
      

      mydepth.out1=c(i, my.vars)
      mydepth.out2=c(i, vs$varselect[,1], vs$md.obj$threshold)
      
      #out1[i,]=mydepth.out1
      #out2[i,]=mydepth.out2
      out1[i,1:length(mydepth.out1)]=mydepth.out1
      out2[i,1:length(mydepth.out2)]=mydepth.out2
      
    }
    
    list(out1, out2)
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("2 argument must be supplied. Data file and output file.", call.=FALSE)
}

datafile = args[1]
#read data   
x=read.table(datafile, header=F)


output = optimalSubsetGenerator(x)

# table with ranks
filename1 = args[2]
write.table((output[1]), filename1, row.names=F, col.names=F, quote=F)

# table with depths
filename2 = args[3]
write.table((output[2]), filename2, row.names=F, col.names=F, quote=F)


