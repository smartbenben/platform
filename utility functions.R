###########################################################################
# check if inputs are correct
# n: patients enrolled in each arm [vector]
# x: successes in each arm [vector]

validnx <- function(n, x) {
  # check validity of n[] patients in each treatment arm, x[] successes
  k <- length(n)                   #  number of treatment arms
  if( length(x) != k || min(n) < 0 || min(x) < 0 ||
      sum(x > n) != 0) return(FALSE)
  return(TRUE)
}
#validnx( c(2,3), c(3,4,5))

# m: number of patients to allocate (scalar)
# k: number of arms (scalar)
#  number of possible outcomes
# choose(m + k - 1, k - 1)  
#cme <- function(m, k){
  # possible outcomes in each arm
 # outcomes <- list()
#  for(i in 1:k){
#    outcomes[[i]] <- 0:m
#  }
#  combns <- expand.grid(outcomes)
#  combns <- combns[rowSums(combns)==m,]
#  rownames(combns) <- NULL
#  return(combns)
#}
#cme(m=4, k=3)
cme <- function(n, m, k){
  # Complete multinomial enumeration.  On successive calls, 
  #   return n[] containing all outcomes of m items into k categories
  #   Begin and end with a vector of all zero's
  if(sum(n==0) == k)           #  is this the initial vector of all zeros?
    return(c(m,  rep(0, k - 1))) # ... then assign the first n[] 
  for(j in 2 : k) {            #  Otherwise update: For each digit...
    n[j] <- n[j] + 1           #    ...increase its value
    s <- sum(n[j : k])         #  sum of this and all to the right
    n[1] <- m - s              #  n[1] is whatever is left over
    if(s <= m)return(n)        #  if n[] is valid, return this value
    n[j] <- 0                  #  Otherwise, reset and look at the next digit
  }                            #  We are done when j runs off the end
  rep(0, k)                    #  Return all zero's to signal completion
}
######################################################################
#  format the output from pda()
buildp <- function(E) {              
  # build table of design outcomes and their respective probabilities
  d <- dim(E)                        #  sparse  (N+1) x ... x (N+1)  array
  k <- length(d)                     #  number of treatment arms
  N <- d[1] - 1                      #  number of patients
  n <- rep(0, k)                     #  initialize loop on treatment alloc's
  tab <- NULL                        #  table to build
  repeat {                           #  loop through all possible allocations
    n <- cme(n, N, k)                #  generate next treatment allocation
    if(max(n) == 0)break             #  are we done generating allocations?
    prob <- E[rbind(n + 1)]          #  probability of this allocation
    tab <- rbind(tab, c(n, prob))    #  build table
  }
  tab <- data.frame(tab)             #  make it pretty
  colnames(tab) <- c(LETTERS[1 : k], "prob")
  tab
}

######################################################################

summarize <- function(tab) {          #  summarizes output from buildp
  # summary statistics from the table of treatment allocations
  # values: sum = sum of probabilities;    exa = expected pt's in each treatment
  summarize <- NULL
  k <- dim(tab)[2] - 1                #  number of treatment arms
  summarize$sum <- sum(tab[ , k + 1]) #  sum of probabilities (hopefully = 1!)
  summarize$expected <- Ea(tab)       #  expected allocation to each arm
  summarize$arms <- k                 #  number of treatment arms
  summarize$N <- sum(tab[1, 1:k])     #  number of patients
  summarize$design <- tab             #  original design matrix
  return(summarize)                   #  and all its components
}




