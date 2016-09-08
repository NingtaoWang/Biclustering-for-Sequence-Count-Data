## A Statistical Biclustering Method for Sequence Count Data

Here is an simple example.
    
    #The function BM.RunEM biclusters the data. 
    #Inputs of BM.RunEM
    #data: a matrix of sequence counts
    #K: number of row clusters
    #L: number of column clusters
    source("BlockMixtureModel_NegativeBinomial_Rcpp.R")
    mu <- array(1:9, dim = c(3, 3))
    r <- 1
    p <- c(0.3, 0.4, 0.3)
    q <- c(0.2, 0.5, 0.3)
    data <- BM.SimuData(n = 500, m = 500, mu, r, p, q)
    result <- BM.RunEM(data, K = 3, L = 3)
