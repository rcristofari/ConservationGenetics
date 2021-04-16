BottleNeckSim <- function(fqs=c(.1, .2, .5), nInd=c(1000, 20, 1000), nGen=c(40, 60, 100), rep=10){
  
require(ggplot2)
require(ggpubr)

bottle <- function(f, nInd, nGen){
  nChrom <- c(rep(nInd[1]*2, nGen[1]),
  rep(nInd[2]*2, nGen[2] - nGen[1]),
  rep(nInd[3]*2, nGen[3] - nGen[2]))
  
  freqs <- numeric(nGen[3])
  freqs[1] <- f
  for(g in 2:nGen[3]) freqs[g]  <- rbinom(1, nChrom[g], freqs[g-1])/nChrom[g]
  return(freqs)
}


drift_block <- function(f, nInd, nGen, rep){
  out <- data.frame()
  nChrom <- c(rep(nInd[1]*2, nGen[1]),
              rep(nInd[2]*2, nGen[2] - nGen[1]),
              rep(nInd[3]*2, nGen[3] - nGen[2]))
  
  for(r in 1:rep){
    G <- 1:nGen[3]
    D <- bottle(f, nInd, nGen)
    N <- rep(sprintf("N=%d", nInd), nGen[3])
    R <- rep(sprintf("f=%3s, replicate_%d", f, r), nGen[3])
    N <- nChrom/2
    out <- rbind(out, data.frame(G, D, N, R, N))
  }
  return(out)
}

drift_series <- function(fqs, nInd, nGen, rep){
  out <- data.frame()
  for(this_f in fqs){
    d <- drift_block(f=this_f, nInd=nInd, nGen=nGen, rep=rep)
    d$freq <- rep(sprintf("f=%3s", this_f), nrow(d))
    out <- rbind(out, d)  
  }
  return(out)
}



fixation_points <- function(d) {
  
  zeroes <- suppressWarnings(na.omit(data.frame(G=as.numeric(by(d, d$R, function(x) ifelse(min(x$G[x$D == 0]) != Inf, 
                                                                                           min(x$G[x$D == 0]), NA), simplify=T)),
                                                f=rep(0, nrow(d)),
                                                freq = substr(names(by(d, d$R, function(x) ifelse(min(x$G[x$D == 0]) != Inf, 
                                                                                                  min(x$G[x$D == 0]), NA), simplify=T)), 1, 5))))
  
  ones <- suppressWarnings(na.omit(data.frame(G=as.numeric(by(d, d$R, function(x) ifelse(min(x$G[x$D == 1]) != Inf, 
                                                                                         min(x$G[x$D == 1]), NA), simplify=T)),
                                              f=rep(1, nrow(d)),
                                              freq = substr(names(by(d, d$R, function(x) ifelse(min(x$G[x$D == 1]) != Inf, 
                                                                                                min(x$G[x$D == 1]), NA), simplify=T)), 1, 5))))
  
  fix <- rbind(zeroes, ones)
  
  return(fix)}


d <- drift_series(fqs=fqs, nInd=nInd, nGen=nGen, rep=rep)
fix <- fixation_points(d)
startpoints <- d[d$G==min(d$G),]
endpoints <- d[d$G==max(d$G),]

pHetStart <- mean(2*startpoints$D*(1-startpoints$D))
pHetEnd <- mean(2*endpoints$D*(1-endpoints$D))
HetLabel <- sprintf("Heterozygosity: \nstart=%#.3f; end=%#.3f", pHetStart, pHetEnd)

g1 <- ggplot(d) +
  theme_bw() +
  theme(legend.position = "None", axis.title.x=element_blank()) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_line(aes(x=G, y=D, group=R, col=freq), lwd=.25) +
  geom_point(data=fix, aes(x=G, y=f, col=freq), size=2) +
  geom_point(data=endpoints, aes(x=G, y=D, col=freq), size=2, shape=18) +
  geom_vline(xintercept=nGen[1:2], col='red', lty='dashed') +
  ylab("Allele frequency") +
  annotate("text", x=0, y=1, label=HetLabel, hjust=0, vjust=1, size=5)

g2 <- ggplot(d[1:max(d$G),]) +
  theme_bw() +
  geom_step(aes(x=G, y=N), lwd=1, col='red', direction='vh') +
  xlab("Generation") + ylab("Population Size")

g <- ggarrange(g1, g2, ncol=1, heights=c(2, 1))  

return(g) }
