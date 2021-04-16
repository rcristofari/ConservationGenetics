splitSim <- function(fqs=c(.1, .2, .5), nInd=c(1000, 50, 500), nGen=c(50, 100), m=c(.1, .1)){
  
  # fqs: a vector of allele frequencies of aribtrary length. Each will have "rep" replicates.
  # nInd: a vector with the number of individuals before the split, in afte rthe split in pop1 and in pop2 
  # nGen: a vector with the time of the split and the end of the simulation
  # m: the migration rates frop pop2 to pop1 and from pop1 to pop2, in proportion of the pop that comes from the other pop at each generation.
  
  require(ggplot2)
  require(ggpubr)
  
  nChrom_0 <- rep(nInd[1]*2, nGen[1])
  nChrom_1 <- rep(nInd[2]*2, nGen[2]-nGen[1])
  nChrom_2 <- rep(nInd[3]*2, nGen[2]-nGen[1])
  
  
  drift <- function(f, nChrom){
    freqs <- numeric(length(nChrom))
    freqs[1] <- f
    for(g in 2:length(nChrom)) freqs[g]  <- rbinom(1, nChrom[g], freqs[g-1])/nChrom[g]
    return(freqs)
  }
  
  
  drift_block <- function(f, nChrom, l){
    out <- data.frame(
      G=0:(length(nChrom)-1),
      D=drift(f, nChrom),
      R=rep(sprintf("locus_%04d", l), length(nChrom)))
    return(out)
  }
  
  drift_series <- function(fqs, nChrom){
    out <- data.frame()
    l<-1
    for(this_f in fqs){
      d <- drift_block(f=this_f, nChrom=nChrom, l=l)
      d$freq <- rep(sprintf("f=%3s", this_f), nrow(d))
      out <- rbind(out, d)
      l <- l+1
    }
    return(out)
  }
  
  # Drift during the ancestral period:
  ancestralDrift <- drift_series(fqs=fqs, nChrom=nChrom_0)
  fqsAtSplit <- ancestralDrift[ancestralDrift$G==max(ancestralDrift$G),"D"]
  
  
  Migr_drift <- function(f, nChrom_1, nChrom_2, m, l){
    freqs_1 <- numeric(length(nChrom_1))
    freqs_2 <- numeric(length(nChrom_2))
    
    freqs_1[1] <- rbinom(1, nChrom_1[1], f)/nChrom_1[1]
    freqs_2[1] <- rbinom(1, nChrom_2[1], f)/nChrom_2[1]
    
    for(g in 2:length(nChrom_1)){
      freqs_1[g]  <- rbinom(1, nChrom_1[g], ((1-m[1])*freqs_1[g-1] + m[1]*freqs_2[g-1]))/nChrom_1[g]
      freqs_2[g]  <- rbinom(1, nChrom_2[g], ((1-m[2])*freqs_2[g-1] + m[2]*freqs_1[g-1]))/nChrom_2[g]
    }
    return(data.frame(G=(0:(length(nChrom_1)-1) + nGen[1]),
                      D1=freqs_1, 
                      D2=freqs_2,
                      R=rep(sprintf("locus_%04d", l), length(nChrom_1))))
  }
  
  Migr_series <- function(fqs, nChrom_1, nChrom_2, m){
    out <- data.frame()
    l<-1
    for(this_f in fqs){
      d <- Migr_drift(f=this_f, nChrom_1, nChrom_2, m, l)
      d$freq <- rep(sprintf("f=%3s", this_f), nrow(d))
      out <- rbind(out, d)
      l <- l+1
    }
    return(out)
  }
  

  # Drift in the resulting populations:
  popDrift <- Migr_series(fqs=fqsAtSplit, nChrom_1, nChrom_2, m)
  
  pop1Drift <- popDrift[,-3]
  pop2Drift <- popDrift[,-2]
  names(pop1Drift)[2] <- "D"
  names(pop2Drift)[2] <- "D"
  
  gdata <- rbind(pop1Drift, ancestralDrift, pop2Drift)
  gdata$group <- c(rep("A_POP1", nrow(pop1Drift)), rep("B_Ancestral", nrow(ancestralDrift)), rep("C_POP2", nrow(pop2Drift)))
                   

  # Fst through time:
  Fst <- numeric(nGen[2]-nGen[1])
  i <- 1
  c1 <- nInd[2]/(nInd[2]+nInd[3])
  c2 <- nInd[3]/(nInd[2]+nInd[3])
  
  for(f in nGen[1]:(nGen[2]-1)){
    this_gen <- popDrift[popDrift$G == f,]
    fst_D1D2 <- function(p1, p2) {
      p <- c1*p1 + c2*p2
      fst <- (p*(1-p) - (c1*p1*(1-p1) + c2*p2*(1-p2))) / p*(1-p)
      return(fst)
    }
    fsts <- ifelse(!is.na(fst_D1D2(this_gen$D1, this_gen$D2)), fst_D1D2(this_gen$D1, this_gen$D2), 0)
    Fst[i] <- mean(fsts)
    i <- i+1
    }
  fstdata <- data.frame(G=0:(nGen[2]-1), Fst=c(rep(NA, nGen[1]), Fst))
  arrows <- data.frame(x0=((nGen[1]+nGen[2])/2)*.95, x1=((nGen[1]+nGen[2])/2)*1.05, ymin=0.2, ymax=0.8,
                        group=factor("B_Ancestral",levels=c("A_POP1","B_Ancestral","C_POP2")))
  
  HetStart <- mean(2*fqs*(1-fqs))
  HetPop1 <- mean(2*pop1Drift[pop1Drift$G==max(pop1Drift$G),"D"]*(1-pop1Drift[pop1Drift$G==max(pop1Drift$G),"D"]))
  HetPop2 <- mean(2*pop2Drift[pop2Drift$G==max(pop2Drift$G),"D"]*(1-pop2Drift[pop2Drift$G==max(pop2Drift$G),"D"]))
  HetData <- data.frame(het = sprintf("H = %#.3f", c(HetStart, HetPop1, HetPop2)),
                        group=as.factor(c("B_Ancestral", "A_POP1", "C_POP2")), 
                        y=1,
                        x=0) 
  
  
  
  g <- ggplot(gdata) +
    theme_bw() +
    theme(legend.position = "None", axis.title.x=element_blank()) +
    coord_cartesian(ylim=c(0, 1)) +
    geom_line(aes(x=G, y=D, group=R, col=freq), lwd=.25) +
    facet_wrap(.~group, ncol=1) +
    geom_vline(xintercept=nGen[1], col='red', lty='dashed', lwd=1) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    ylab("Allele frequency") + 
    geom_segment(data=arrows, aes(x=x0, xend=x0, y=ymin, yend=ymax), lwd=10*m[1], arrow=arrow(length=unit(0.4,"cm"))) +
    geom_segment(data=arrows, aes(x=x1, xend=x1, y=ymax, yend=ymin), lwd=10*m[2], arrow=arrow(length=unit(0.4,"cm"))) +
    geom_text(data=HetData, aes(x=x, y=y, label=het), hjust=0, vjust=1, size=4)
  
  
  gfst <- ggplot(fstdata) +
    theme_bw() +
    theme(legend.position = "None", axis.title.x=element_blank()) +
    geom_line(aes(x=G, y=Fst), col='red') +
    geom_hline(yintercept=0, lty='dashed') +
    geom_vline(xintercept=nGen[1], col='red', lty='dashed', lwd=1) +
    ylab("Fst") + xlab("")
    
  p <- ggarrange(g, gfst, ncol=1, heights=c(3, 1))  
  
  return(p) }
