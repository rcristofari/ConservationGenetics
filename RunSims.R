

DriftSim(fqs=c(.1, .2, .5), nInd=1000, nGen=100, rep=50)


#ggsave("~/Desktop/ConservationGenetics/Images/drift_5.pdf", width=12, height=6)


BottleNeckSim(fqs=c(.1, .2, .5), nInd=c(10000, 50, 1000), nGen=c(40, 60, 200), rep=50)


##############################################################################################################


#ggsave("~/Desktop/ConservationGenetics/Images/drift_5.pdf", width=12, height=6)

# function to show allele freqs before and after bottleneck / low period.

#Cheetah
g <- BottleNeckSim(fqs=c(.1, .2, .5), nInd=c(142500, 37700, 100), nGen=c(35000, 39000, 43000), rep=50)
BottleNeckSim(fqs=c(.1, .2, .5), nInd=c(142500, 37700, 100), nGen=c(35000, 39000, 43000), rep=50)
ggsave("~/Desktop/ConservationGenetics/cheetah.pdf", width=12, height=8)

# White tailed deer:
BottleNeckSim(fqs=c(.1, .2, .5), nInd=c(1000000, 5, 100000), nGen=c(50, 51, 100), rep=50)

# With a more realistic distribution:
rawspectrum <- rexp(1000, 9)
spectrum <- ifelse(rawspectrum>=1, 1-1e-6, rawspectrum)
BottleNeckSim(fqs=spectrum, nInd=c(100000, 5, 10000), nGen=c(10, 11, 100), rep=1)


BottleNeckSim(fqs=spectrum, nInd=c(100000, 5, 1000), nGen=c(10, 15, 100), rep=1)





SizeChangeSim(fqs=c(.1, .2, .5), nInd=c(1000, 500, 100, 10, 20, 40, 80, 100, 200, 400, 800), nGen=c(20, 25, 30, 35, 36, 37, 38, 39, 40, 45, 100), rep=100)


rawspectrum <- rexp(100, 9)
spectrum <- ifelse(rawspectrum>=1, 1-1e-6, rawspectrum)

splitSim(fqs=spectrum, nInd=c(1000, 20, 1000), nGen=c(15, 200), m=c(0.025, 0))
#ggsave("~/Desktop/ConservationGenetics/Images/Fst_1.pdf", width=12, height=8)
