het <- function(x) 2*x*(1-x)
grDevices::pdf('/Users/Robin//Desktop/ConservationGenetics/Images/het.pdf', width=8, height=6)
curve(het, from=0, to=1, add=F, xlab="Allele frequency", ylab="Expected heterozygosity", lwd=3)
abline(v=.5, col='red', lty='dashed', lwd=3)
dev.off()


fixation <- function(p, N) -4*N*(p/(1-p))*log(p)

p <- seq(0, 1, .01)
N <- 1:10000

grDevices::pdf('/Users/Robin//Desktop/ConservationGenetics/Images/Fix.pdf', width=8, height=6)
plot(seq(0, 1, .01), fixation(p, 100), type='l', xlab='Starting frequency', 
     ylab='Generations generations', lwd=2, ylim=c(0, 5000))
lines(seq(0, 1, .01), fixation(p, 200), lwd=2)
lines(seq(0, 1, .01), fixation(p, 500), lwd=2)
lines(seq(0, 1, .01), fixation(p, 1000), lwd=2)
lines(seq(0, 1, .01), fixation(p, 2000), lwd=2)
abline(v=0, lty='dashed')
dev.off()
