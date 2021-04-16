SizeChangeSim <- function(fqs=c(.1, .2, .5), nInd=c(10, 20, 30, 40, 50), nGen=c(5, 10, 15, 20, 25), rep=1, doPlot='full', fixedSites=F){
  
  # fqs: a vector of allele frequencies of aribtrary length. Each will have "rep" replicates.
  # nInd: a vector with the number of individuals at each period
  # nGen: a vector with the population size boundaries. Its length is the same as nInd (the first boundary is implicitly 0)
  
  require(ggplot2)
  require(ggpubr)
  
  # Test the input:
  if(class(fqs) != "numeric") stop("The allele frequencies must be supplied as a numeric vector")
  if(class(fqs) != "numeric") stop("The population sizes must be supplied as a numeric vector")
  if(class(fqs) != "numeric") stop("The time points must be supplied as a numeric vector")
  if(length(nInd) != length(nGen)) stop("The nInd and nGen vectors don't match")
  if(!(doPlot %in% c("full", "spectrum"))) stop("doPlot must be 'full' or 'spectrum'")
  
  bottle <- function(f, nChrom){
    freqs <- numeric(length(nChrom))
    freqs[1] <- f
    for(g in 2:length(nChrom)) freqs[g]  <- rbinom(1, nChrom[g], freqs[g-1])/nChrom[g]
    return(freqs)
  }
  
  
  drift_block <- function(f, nInd, nGen, rep){
    out <- data.frame()
    nChrom <- rep(nInd[1]*2, nGen[1])
    for(g in 2:length(nInd)) nChrom <- c(nChrom, rep(nInd[g]*2, nGen[g] - nGen[g-1]))
    
    for(r in 1:rep){
      G <- 0:(length(nChrom)-1)
      D <- bottle(f, nChrom)
      R <- rep(sprintf("f=%3s, replicate_%d", f, r), length(nChrom))
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
  
  if(doPlot =="full") {
  g1 <- ggplot(d) +
    theme_bw() +
    theme(legend.position = "None", axis.title.x=element_blank()) +
    coord_cartesian(ylim=c(0, 1)) +
    geom_line(aes(x=G, y=D, group=R, col=freq), lwd=.25) +
    geom_point(data=fix, aes(x=G, y=f, col=freq), size=2) +
    geom_point(data=endpoints, aes(x=G, y=D, col=freq), size=2, shape=18) +
    geom_vline(xintercept=nGen, col='red', lty='dashed') +
    ylab("Allele frequency") +
    annotate("text", x=0, y=1, label=HetLabel, hjust=0, vjust=1, size=5)
  
  g2 <- ggplot(d[1:max(d$G),]) +
    theme_bw() +
    geom_step(aes(x=G, y=N), lwd=1, col='red', direction='hv') +
    xlab("Generation") + ylab("Population Size")
  
  g <- ggarrange(g1, g2, ncol=1, heights=c(2, 1))  
  
  } else {
  
  sdata <- data.frame(fqs=c(fqs, endpoints$D), class=as.factor(c(rep("A_start", length(fqs)), rep("B_end", nrow(endpoints)))))
  if(fixedSites==F) sdata <- sdata[which(sdata$fqs > 0 & sdata$fqs < 1),]
  
   g1 <- ggplot(sdata) +
     theme_bw() +
     geom_histogram(aes(x=fqs), fill='lightgray', col='black', bins=50) +
     facet_wrap(.~class, ncol=1, scales = "free_y") +
     theme(strip.background = element_blank(), strip.text.x = element_blank()) +
     xlab("Allele frequency") + ylab ("") +
     theme(panel.grid = element_blank())
   g2 <- ggplot(d[1:max(d$G),]) +
     theme_bw() +
     geom_step(aes(x=G, y=N), lwd=1, col='red', direction='hv') +
     xlab("Generation") + ylab("Population Size")
   
   g <- ggarrange(g1, g2, ncol=1, heights=c(2, 1))  
   
  
  }
  
  
  return(g) 
  }

rawspectrum <- rexp(10000, 9)
spectrum <- ifelse(rawspectrum>=1, 1-1e-6, rawspectrum)

# Bottleneck
SizeChangeSim(fqs=spectrum, nInd=c(1000, 500, 100, 10, 20, 40, 80, 100, 200, 400, 800), nGen=c(20, 25, 30, 35, 36, 37, 38, 39, 40, 45, 100), rep=1, doPlot='spectrum', fixedSites = T)

# Exponential growth
SizeChangeSim(fqs=spectrum, nInd=c(1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000), nGen=c(20, 40, 60, 80, 100, 120, 140, 160), rep=1, doPlot='spectrum', fixedSites = T)

# Gradual collapse
SizeChangeSim(fqs=spectrum, nInd=c(128000, 64000, 32000, 16000, 8000, 4000, 2000, 1000), nGen=c(20, 40, 60, 80, 100, 120, 140, 160), rep=1, doPlot='spectrum', fixedSites = T)
