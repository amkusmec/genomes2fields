#### Additional File 2
# Teschendorff et al. (2018) Tensorial blind source separation for improved 
# analysis of multi-omic data. Genome Biology, 19: 76.

# Summary: this file contains R-scripts for implementing the tPCA, tWFOBI, 
# tWJADE, JIVE, CCA, sCCA, iCluster and PARAFAC algorithms (as used for the 
# simulation model and application to real omic data). We also provide a 
# function to estimate sensitivity.

# Libraries needed
library(isva); # necessary for the RMT function which estimates number of significant components of variation
library(tensorBSS); # implements tPCA, tFOBI and tJADE


#### INPUT Arguments

# Some of the input objects are common to all methods, others are unique to 
# each method. The common ones are described here:

# data: this is the data object, which contains the multi-way data in two 
# different formats. The "A" entry of data (data$A) gives the array or 
# data-tensor format, whereas the "L" entry of data (data$L) gives the data in 
# list format. In the former case, and in our specific applications, the first 
# mode defines the tissue, cell or data-type, the second mode defines the samples 
# and the third mode the features (e.g. CpGs or genes). In the latter case, each 
# list entry corresponds to the cell/tissue or data-type and consists of the 
# data-matrix which rows representing features and columns representing samples. 

# dim: this is a vector which contains the number of significant components of 
# each data matrix to search for, and is typically obtained by applying RMT to 
# each separate data/tissue-type matrix (i.e. to the individual entries of data$L 
# above).

# dJ: the number of significant components of joint variation across data/tissue
# types. Typically we define this to be the RMT estimate applied on the joint 
# matrix, obtained by stacking up the matrices in each entry of data$L.

# tpJV.idx: this is an index vector labeling the true positive features 
# associated with a factor of interest, which we know drives joint variation in 
# the data, and which therefore the algorithm should capture. These indices must 
# be in the order of the entries in the 3rd mode of the data$A object, or 
# alternatively in the same order as the rows of the individual data matrices in 
# data$L.

# topN: the number of top-ranked features to select from each inferred 
# component. It must be specified and by default it equals the number of true 
# positives.

# maxiter: The maximum number of iterations used in the algorithms.


# Additional input:

# rankJ: An integer giving the number of components of significant joint 
# variation to use when implementing JIVE, if known. If not given, this will be 
# calculated using the chosen method. If the method is "given" then the default 
# is 1.

# rankA: A vector giving the individual ranks (i.e. number of significant 
# components of individual variation of each data matrix), if known. If not 
# given, this will be calculated using the chosen method. If the method is 
# "given" then the default is rep(1, length(data)).

# method: A string with the method to use for rank selection. Possible options 
# are "given", "perm", and "bic". The default is "perm". If ranks are known, you 
# should use "given".

# npermEstDim: this is the number of permutations to use to determine 
# significance of the amount of variance carried by each of the canonical vectors 
# (when using CCA/sCCA).


#### R-SCRIPTS

### Auxiliary function to estimate dimensionality parameters
### Estimate dim and dJ by RMT
EstDim <- function(data){
  data.l <- data$L;
  nt <- length(data.l);
  ### estimate joint variation
  data.m <- data.l[[1]];
  for(i in 1:(nt-1)){
    data.m <- rbind(data.m,data.l[[i+1]]);
  }
  sd.v <- sqrt(apply(data.m,1,var));
  d <- EstDimRMT((data.m - rowMeans(data.m))/sd.v,plot=FALSE)$dim
  dim <- vector();
  for(j in 1:nt){
    dim[j] <- EstDimRMT( data.l[[j]]-rowMeans(data.l[[j]]), plot = FALSE)$dim;
  }
  return(list(dJ=d, dim=dim));
}

DoTPCA <- function(data,dim){

    data.a <- data$A;
    nt <- dim(data.a)[1];
    ng <- dim(data.a)[3];
    dim.v <- c(nt,max(dim));
    tpca.o <- tPCA(data.a,d=dim.v);

    projS.lm <- list();
    for(t in 1:nt){
     projS.lm[[t]] <- tpca.o$S[t,,];
    }
    
    return(list(projS=projS.lm,S=tpca.o$S,U=tpca.o$U,nt=nt,ng=ng));

}

#### tensorial ICA: tWFOBI and tWJADE

DoTICA <- function(data,dim,method=c("FOBI","JADE")){
    data.a <- data$A;
    nt <- dim(data.a)[1];
    ng <- dim(data.a)[3];
    dim.v <- c(nt,max(dim));
    cdata.a <- tensorCentering(data.a)
    tpca.o <- tPCA(cdata.a,d=dim.v);
    pdata.a <- tensorTransform(cdata.a,t(tpca.o$U[[2]]),2); ## whiten the data

    if(method=="FOBI"){    
      tica.o <- tFOBI(pdata.a);
    }
    else {
      tica.o <- tJADE(pdata.a);
    }
    
    projS.lm <- list();
    for(t in 1:nt){    
      projS.lm[[t]] <- tica.o$S[t,,];
    }

    return(list(projS=projS.lm,S=tica.o$S,U=tica.o$W,nt=nt,ng=ng));

}
