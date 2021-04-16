DriftSim <- function(fqs=c(.1, .2, .5), nInd=1000, nGen=100, rep=50){

require(ggplot2)

drift <- function(f, nInd, nGen){
  nChrom <- 2*nInd
  freqs <- numeric(nGen)
  freqs[1] <- f
  for(g in 2:nGen) freqs[g]  <- rbinom(1, nChrom, freqs[g-1])/nChrom
  return(freqs)
}


drift_block <- function(f, nInd, nGen, rep){
  out <- data.frame()
  
  for(r in 1:rep){
    G <- 1:nGen
    D <- drift(f, nInd, nGen)
    N <- rep(sprintf("N=%d", nInd), nGen)
    R <- rep(sprintf("f=%3s, replicate_%d", f, r), nGen)
    out <- rbind(out, data.frame(G, D, N, R))
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

g <- ggplot(d) +
  theme_bw() +
  theme(legend.position = "None") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_line(aes(x=G, y=D, group=R, col=freq), lwd=.25) +
  geom_point(data=fix, aes(x=G, y=f, col=freq), size=2) +
  geom_point(data=endpoints, aes(x=G, y=D, col=freq), size=2, shape=18) +
  xlab("Generation") + ylab("Allele frequency") +
  annotate("text", x=0, y=1, label=sprintf("N = %d", nInd), hjust=0, vjust=1, size=6.5) +
  annotate("text", x=0, y=.9, label=HetLabel, hjust=0, vjust=1, size=6.5)

return(g)}