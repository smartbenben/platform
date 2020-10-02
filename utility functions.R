###########################################################################
# arms: number of treatment arms
# size: number of maximum sample sizes in each arm
enumerate=function(arms, size){
  tmp <- list()
  
  for(j in 1:arms){
    tmp[[j]] <- 0:size[j]
  }
  
  return(as.matrix(expand.grid(tmp)))
}
#############################################
##generate all the possible allolations (n) and outcomes(x)
gen_nx <- function(N, R){
  ## all the possible allocations with fixed N
  ## there are a total of choose (N + R - 1, R - 1)  
  ## allocations
  allocs <- enumerate(arms = R, size = rep(N, R))
  allocs <- allocs[rowSums(allocs)==N, ]
  colnames(allocs) <- LETTERS[1:R]
  rownames(allocs) <- NULL
  
  mat <- NULL
  for(i in seq_len(nrow(allocs))) {
    ## nout: number of possible outcomes
    ## given an allocation vector (n1, n2)
    ## the number of possible outcomes is (n1+1)*(n2+1)
    nout <- prod(allocs[i, ]+1)
    lvec <- rep(allocs[i, ], nout)
    mat1 <- matrix(lvec, nrow = nout, byrow = T)
    mat2 <- enumerate(arms = R, size = allocs[i, ])
    mat <- rbind(mat, cbind(mat1, mat2))
  }
  
  mat
}

#old <- gen_nx(N = 2, R = 2)
#new <- gen_nx(N = 3, R = 2)

#############################################
## find orgins of a new (n,x) in the most recent past
## new_v: a row of the new[dataframe] of all (n,x)
## old[dataframe] of all (n,x) 
tracer <- function(new_v, old, R) {
  nvec <- new_v[1:R]
  xvec <- new_v[-c(1:R)]
  
  ## each nvec has R origins in the lastest past
  ## the new patient to arm r might response or not
  ori <- foreach(r=1:R,.combine = "rbind")%do%{
    v1 <- rep(0, R)
    v1[r] <- 1
    rbind( c(nvec-v1, xvec-rep(0, R)), 
           c(nvec-v1, xvec-v1) )
  }
  rownames(ori) <- NULL
  
  ori
}

#tracer(new_v = new[10, ], old = old, R = 2)
