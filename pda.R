# Closed-form probability distribution of treatment allocation

# N: number of patients (scalar)
# p: success probabilities (vector)
# priora, priorb: priors for beta(alpha, beta)

pda <- function(N, p, nmin, early.stop=FALSE, 
                priora=NULL, priorb=NULL) {
  #  number of treatment arms
  k <- length(p)   
  print(noquote(c("Number of design points", choose(N + k - 1, k - 1))))
  #  set for first patient
  EN <- firstEN(p, nmin, early.stop,
                priora, priorb)  
  p_alloc <- c(1, Ea(buildp(EN)))
  #  expected allocation
  EaN <- c(1, Ea(buildp(EN)))   
  #  we are done if N = 1
  if(N == 1) 
    return(summary(buildp(EN))) 
  #  for each additional patient
  for(i in 2 : N) {    
    #  tree set for next patient
    EN <- nextEN(EN, p, nmin, early.stop,
                 priora, priorb) 
    #  build table of expected treatment arm allocation
    alloc <- buildp(EN)
    
    p_new <- rowSums(apply(alloc, 1, function(x){
      x[length(p)+1]*ERN1(n=unlist(x[1:length(p)]), p, 
                          nmin, early.stop, priora, priorb)
    }))
    p_alloc <- rbind(p_alloc, c(i,p_new))
    EaN <- rbind(EaN,c(i, Ea(alloc))) 
    
    #  progress report 
    print(noquote(c("N: ", i)))  
    #  print it out before completion of loop
    flush.console()                 
  }    
  
  #  output, summarized
  pda <- summarize(buildp(EN))      
  colnames(EaN) <- c("N", LETTERS[1 : k])
  #  expected treatment allocations
  pda$EaN <- EaN   
  pda$p_alloc <- p_alloc
  return(pda)
}


