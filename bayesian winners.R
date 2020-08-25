# alphas, betas: parameters of beta distributions
# alphas: numbers of prior successes in each arm (vector)
# betas: numbers of prior failures in each arm (vector)

#################################################################
# compute Prob(p(A)>p(B))
pAB <- function(alphas, betas){
  i <- 0:(alphas[1]-1)
  num <- beta(alphas[2]+i,sum(betas))
  den <- (betas[1]+i) * beta(1+i, betas[1]) * 
    beta(alphas[2],betas[2])
  sum(num/den)
}
#############################################################
# compute Prob(A is the best among A, B, C)
pABC <- function(alphas, betas){
  mat <- expand.grid( 0:(alphas[2]-1), 0:(alphas[3]-1) )
  i <- mat[,1]
  j <- mat[,2]
  num <- beta(i+j+alphas[1], sum(betas))
  den <- (betas[2]+i)*(betas[3]+j)*
    beta(1+i, betas[2])*
    beta(1+j, betas[3])*
    beta(alphas[1], betas[1])
  sums <- sum(num/den)
  1-pAB(alphas[c(2,1)],betas[c(2,1)])-
    pAB(alphas[c(3,1)],betas[c(3,1)])+sums
}
#############################################################
# calculate probability the 1st beta population is best
win1 <- function(alphas=NULL, betas=NULL){  
  if(length(alphas)==2){
    p1_best <- pAB(alphas, betas)
  }
  if(length(alphas)==3){
    p1_best <- pABC(alphas, betas)
  }
  if(length(alphas)>3){
    #  objective function for integral
    obj <- function(x) {             
      obj <- NULL  
      #  vector x
      for(i in 1 : length(x)) { 
        obj <- c(obj, dbeta(x[i], alphas[1], betas[1]) *              
                   prod(pbeta(x[i], alphas[-1], betas[-1]))) 
      }  
      obj
    }
    p1_best <- integrate(obj, lower = 0, upper = 1)$value
  }
   return(p1_best)
}           

# # # # # #  test programs
#win1(c(3, 4), c(1, 1))
#win1(c(3, 4, 5), c(1, 1, 3))


######################################################################
#  Probability the i-th beta population is the largest
probwin <- function(alphas, betas) {
  # check validity of parmameters
  if((length(alphas) != length(betas)) ||  
     (min(alphas) <= 0) || (min(betas) <= 0)) 
    warning("Incompatible parameters")
  #  accumulate probabilities
  probwin <- NULL         
  #  number of populations
  la <- 1 : length(alphas) 
  if(length(la)==1){
    probwin <- 1
  }else{
    for(i in la) {     
      #  compare each to all others
      probwin <- c(probwin, win1(alphas[c(i, la[-i])], betas[c(i, la[-i])]))
    }
  }
  probwin
}

# # # # # test programs
# alphas <- c(4, 4, 8, 12, 16)       
# betas <-  c(3, 3, 6,  9, 12)
#  all the means are the same
# alphas / (alphas + betas)   
# variances are different
# (alphas*betas)/( ((alphas+betas)^2) * (alphas+betas+1))
# the one with larger variance is more likely to be the best
# probwin(alphas, betas)              
# sum(probwin(alphas, betas))   
##asign patients to an arm with the most promising beta distrubtion
RBayeswin <- function(n, x, priora=NULL, priorb=NULL) {  
  #  Bayesian play the winner           
  alphas <- x + priora              
  betas <- n - x + priorb           
  return(probwin(alphas, betas))     
}

