rm(list = ls())
setwd("/Users/ww369/Dropbox/platform trials/R Code")

source("bayesian winners.R")
source("utility functions.R")
source("next.alloc.R")
source("pda.R")

N=30
p=c(.4, .7)

out1 <- pda(N=30, 
            p=c(.4, .7), 
            early.stop = FALSE,
            priora=rep(1,2),
            priorb=rep(1,2)) 

out2 <- pda(N=30, 
            p=c(.4, .7), 
            early.stop = TRUE,
            priora=rep(1,2),
            priorb=rep(1,2)) 

#pdf(file = "rate.pdf")
par(mfrow=c(1,2), mar=c(4,4,1,1))
plot(xlim = c(1, N + 1), ylim = c(0,N+1), 
     type = "n", x=0, y=0,
     xlab = "Sample Size", 
     ylab="Expected Allocation",
     las=1, cex.lab = 1.25)
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
text(labels = LETTERS[1:2], x = rep(N + 1,2),
    y = out1$expected)
legend("topleft", lty=1, lwd=2,
       col=c("red","black"),
       c("without early stopping", 
         "with early stopping"), bty="n")

plot(xlim = c(1, N + 1), ylim = c(0,1), 
     type = "n", x=0, y=0,
     xlab = "Sample Size", 
     ylab="Allocation Probability",
     las=1, cex.lab = 1.25)
for(j in 1 : length(p)){
  lines(out1$p_alloc[,1], out1$p_alloc[,j+1], 
        type = "l", col = "red", lwd = 2)
  lines(out1$p_alloc[,1], out2$p_alloc[,j+1], 
        type = "l", col = "black", lwd = 2)
  #lines(ptw.out$EaN[, 1], ptw.out$EaN[, j+1], type = "l", col = "blue",
  #      lwd = 2)
}
abline(h=c(0.1, 0.9), col="grey", lty=2)
text(labels = LETTERS[1:2], x = rep(N + 1,2),
     y = out1$p_alloc[N,])
#dev.off()
out1$expected

out2$expected

##what is alpha, power
out1$p_alloc[N,]
out2$p_alloc[N,]
