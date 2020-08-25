# calculate the probability of allocating the next patient

# to each arm, given the current n and p
# n: the number of patients in each arm (vector)
# nmax: the maximum number of patients each arm 
# is allowed to accrue (vector)
# p: the probability of response in each arm (vector)
# early.stop: stop early for futility
# 1- lambda: posterior desicion bound at last patient
# gam: posterior bound is a function of n^gam
ERN1 <- function(n, p, 
                 early.stop, 
                 priora=NULL, priorb=NULL){
  if(min(n) < 0 ||  length(p) != length(n)
     ||  min(p) < 0  ||  max(p) > 1)
    warning("Incompatible parameters") 
  
  combns <- list()
  for(i in 1:length(n)){
    combns[[i]] <- 0:n[i]
  }
  # matrix containing all possile outcomes in each arm
  mat <- expand.grid(combns)
  
  p_next <- rep(0, length(n))

  #  don't count until each arm has this many pts 
  #  arms with too few patients
  #islow <- (n < min_n) 
  #  are there arms with too few patients?
  #if(sum(islow) > 0) { 
    # randomize to those arms with too few pts
   # p_next[islow] <- 1 / sum(islow) 
    #return(p_next)
  #}        
  for(i in 1:nrow(mat)){
    x <- unlist(mat[i,])
    r <- RBayeswin(n, x, priora, priorb)
    w <- prod(dbinom(x, n, p))
    
    if(early.stop){
      #theta_e <- 1-lambda*(n/nmax)^gam
      #theta_e <- ifelse(n<nmax, 1, 1-lambda)
      #theta_f <- (1-lambda)*(n/nmax)^gam
      theta_f <- 0.4
      stp <- which(r<theta_f)
      con <- which(r>theta_f)
      #stp <- which(r>theta_e | r<theta_f)
      #con <- which(r<theta_e & r>theta_f)
      if(length(stp)>0){
        r[stp] <- 0
        if(length(con)>0){
          r[con] <- RBayeswin(n[con], x[con], 
                              priora[con], 
                              priorb[con])
        }
      }
    }
    p_next <- p_next+r*w
  }
  return(p_next)
}

#ERN1(n=c(2,3,4), p=c(.1, .2, .15),  
 #   priora=rep(1,3), priorb=rep(1,3))
#ERN1(n = c(55, 70), p=c(.65, .6), 
#    priora=rep(1,2), priorb=rep(1,2))

# # # # #
#  build next level branch via n[] from EN to ENR1
branch <- function(EN, er, n, EN1) {
  #  number of treatment arms
  k <- length(n)                           
  if(length(n) != length(er)) 
    warning("Incompatible parameters")    #  error check
  I <- diag(rep(1, k))                     #  identity matrix
  for (i in 1 : k) {                       #  loop through all tree branches
    np <- n + I[ , i]                      #  location of tree branch
    npi <- rbind(np + 1)                   #  branch expressed as a subscript
    EN1[npi] <- EN1[npi] + er[i] * EN[rbind(n + 1)]    # accumulate branch
  }
  return(EN1)                              #  return updated EN1
}

#  allocation EN tree for first patient
firstEN <- function(p, early.stop, 
                    priora=NULL, priorb=NULL) {
  #  number of treatment arms
  k <- length(p) 
  #  plan for 1st patient
  e <- ERN1(rep(0, k), p, 
            early.stop, 
            priora, priorb)          
  #  build design matrix for one patient
  EN <- array(0, dim = rep(2, k))     
  #  over all treatment arms
  for(j in 1 : k) {                   
    n <- rep(0, k)
    n[j] <- 1
    EN[rbind(n + 1)] <- e[j] 
  }
  return(EN)
}
#allocate to 1st arm EN[2,1,1]
#allocate to 2ed arm EN[1,2,1]
#allocate to 3rd arm EN[,1,1,2]
#.....
######################################################################
#   Allocation tree for the next patient
nextEN <- function(EN, p, 
                   early.stop, 
                   priora=NULL, priorb=NULL) {
  #  number of treatment arms
  k <- length(dim(EN))  
  #  number of patients already allocated
  N <- dim(EN)[1] - 1     
  # tree for patient #  N + 1
  EN1 <- array(0, dim = rep(N + 2, k))  
  n <- rep(0, k)
  repeat {                        #  loop on every n[] in EN
    n <- cme(n, N, k)             #  generate the next n[]
    if(max(n) == 0)break          #  are we done for this N of patients?
    er <- ERN1(n, p, 
               early.stop, 
               priora, priorb)           #  allocation prob of next patient at n[]
    EN1 <- branch(EN, er, n, EN1) #  build next level of branches
  }                               #  loop on all n[]
  return(EN1)                     #  tree for next patient
}

Ea <- function(tab) {             #  expected allocation to each treatment
  k <- dim(tab)[2] - 1            #  number of treatment arms
  Ea <- rep(k, 0)                       
  for (i in 1 : k) { 
    Ea[i] <- sum( tab[ , i] * tab[ , k + 1])
  }
  return(Ea)
}