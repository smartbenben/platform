rm(list = ls(all = TRUE))
Packages <- c("foreach", "parallel", "doMC")
lapply(Packages, library, character.only = TRUE)
library(prodlim)
cores <- detectCores()
registerDoMC(cores)

setwd("/Users/ww369/Dropbox/platform trials/R Code")

source("bayesian winners.R")
source("utility functions.R")

# Nmax: the maximum sample size for the trial
# p[vector]: response rate for each treatment arm

# tau: soft thereshold for burn-in
# tau=NULL, equal allocation
# tau=0, Adaptive randomization (AR) starts from the first pt
# larger tau leads to more burn in

# efficacy boundary: 1-delta_e*(N/Nmax)^gam
# futility boundary: delta_f*(N/Nmax)^gam
# to turn off efficacy stop: set delta_e=0
# to turn off futility stop: set delta_f=0

Nmax = 10; p = c(0.1, 0.1); tau = 0
delta_e = 0; delta_f = 0; gam = 1

recursive_ptw <- function(Nmax, p, tau = NULL, delta_e = 0,
                    delta_f = 0, gam=0) {
  # number of arms (scalar)
  R <- length(p)
  EN <- pse <- psf <- NULL
  plus1 <- diag(nrow = R, ncol = R)
  for (N in 1:Nmax) {
    ## efficacy boundary 
    theta_e <- 1-delta_e*(N/Nmax)^gam
    ## futility boundary 
    theta_f <- delta_f*(N/Nmax)^gam
    
    ## generate all possible n and x, given N and R
    mat <- gen_nx(N = N, R = R)
    new <- NULL
    # i: one combination of (x, n)
    for(i in seq_len(nrow(mat))) {
        new_v <- mat[i, ]
        ## allocation vector n
        nvec <- new_v[1:R]
        ## one of the outcome vector x, given n
        xvec <- new_v[-c(1:R)]
        ## probability of being the best treatment
        p_best <- RBayeswin(n = nvec, x = xvec, 
                            priora = rep(1, R),
                            priorb = rep(1, R))
        nst <- which(p_best<theta_e & p_best>theta_f)
        winner <- loser <- rep(0, R)
        
        ## equal allocation for the first pt
        if(N==1){
          p_recur <- (1/R)*prod(dbinom(xvec, nvec, p))
        }
        if(N>1){
          ## trace the origins of (x, n) is the most recent past
          ori <- tracer(new_v = new_v, old = old[ , 1:(2*R)], R = R)
          idx <- row.match(as.data.frame(ori), old[ , 1:(2*R)], 0)
          ori <- old[idx, ]
          ## Cxn is the probability of observing x, n
          old_recur <- unlist(ori$Cxn)
          old_alloc <- ori[, grep("to", names(old))]
          p_recur <- 0
          ## (x, n) has nrow(ori) origins
          for(j in 1:nrow(ori)){
            # which arm the Nth pt is allocated to
            index <- which((nvec-ori[j, 1:R])==1)
            # is the Nth pt a response or not
            response <- unlist(xvec-ori[j,(R+1):(2*R)])[index]
            p_recur <- p_recur + dbinom(response, 1, p[index])*old_recur[j]*old_alloc[j,index]
            ## conditional probs of being a winner or loser
            winner <- winner + ifelse(p_best>theta_e, 1, 0)*dbinom(response, 1, p[index])*old_recur[j]*old_alloc[j,index]
            loser <- loser + ifelse(p_best<theta_f, 1, 0)*dbinom(response, 1, p[index])*old_recur[j]*old_alloc[j,index]
          }
        }
        
        p_alloc <- rep(0, R)
        ## if there is only one arm left, continue until Nmax
        if(length(nst)==1){
          p_alloc[nst] <- 1
        }
        if(length(nst)>1){
          num <- (p_best[nst])^((0.5*N/Nmax)^tau)
          den <- sum(p_best[nst]^((0.5*N/Nmax)^tau))
          ## the updated prob of allocation based on curent data
          p_alloc[nst] <- num/den
        }
        
        new <- rbind(new, c(mat[i, ], p_alloc, p_recur, winner, loser))
    }
    rownames(new) <- NULL
    new <- as.data.frame(new)
    name_n <- paste("n", LETTERS[1:R], sep="")
    name_x <- paste("x", LETTERS[1:R], sep="")
    name_pi <- paste("to", LETTERS[1:R], sep="")
    names(new) <- c(name_n, name_x, name_pi, "Cxn",
                    paste("winner", LETTERS[1:R], sep=""),
                    paste("loser", LETTERS[1:R], sep=""))
    
    
    ## without early stopping, should sum up to 1
    print(c(N, sum(new$Cxn)))
    cxn <- new$Cxn/sum(new$Cxn)
    EN <- rbind(EN, c(N, colSums(new[,(1:R)]*cxn)))
    
    win <- new[,grep("winner", names(new))]
    los <- new[,grep("loser", names(new))]
    
    pse <- rbind(pse, c(N, colSums(win)))
    psf <- rbind(psf, c(N, colSums(los)))
    
    old <- new
  }
  print(old)
  return(list("EN" = EN, "pse"=pse, "psf"=psf))
}

design <- recursive_ptw(Nmax =10, p = c(0.1, 0.4), tau = 1,
                  delta_e = 0.2, delta_f = 0.2, gam=0.5)
design$EN
round(design$pse,3)
round(design$psf,3)
round(colSums(design$pse),3)
round(colSums(design$psf),3)


