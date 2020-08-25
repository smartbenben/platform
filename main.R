rm(list = ls())
setwd("/Users/ww369/Dropbox/platform trials/R Code")

source("bayesian winners.R")
source("utility functions.R")
source("next.alloc.R")
source("pda.R")

N=20
p=c(.5, .7)

out1 <- pda(N=20, theta_f=0.20, 
            p=c(.5, .7), 
            early.stop = FALSE,
            priora=rep(1,2),
            priorb=rep(1,2)) 

out2 <- pda(N=20, theta_f=0.20,
            p=c(.5, .7), 
            early.stop = TRUE,
            priora=rep(1,2),
            priorb=rep(1,2)) 

#pdf(file = "rate.pdf")
plot(xlim = c(1, N + .5), ylim = c(0,N+1), 
     type = "n", x=0, y=0,
     xlab = "Sample Size", 
     ylab="",
     cex.lab = 1.25)
#lines(1 : N, out1$EaN  %*% c(0,p), 
#      type = "l", col = "black", lwd = 4)
#lines(1 : maxn, ptw.out$EaN %*% c(0,p), type = "l", col = "green", lwd = 4)
for(j in 1 : length(p)){
  lines(out1$EaN[, 1], out1$EaN[, j+1], 
        type = "l", col = "red", lwd = 2)
  lines(out2$EaN[, 1], out2$EaN[, j+1], 
        type = "l", col = "black", lwd = 2)
  #lines(ptw.out$EaN[, 1], ptw.out$EaN[, j+1], type = "l", col = "blue",
  #      lwd = 2)
}
#text(labels = LETTERS[1:2], x = rep(N + 1,2),
#     y = c(6.3, 8))
#dev.off()
out1$expected

out2$expected




