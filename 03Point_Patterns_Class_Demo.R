## ----message=FALSE-----------------------------------------------------------
library(tidyverse)
library(spatstat)


## ----------------------------------------------------------------------------
n <- 100
patternA <- data.frame(x=rnorm(n),y=rnorm(n),id="A")
patternB <- data.frame(x=runif(n),y=runif(n),id="B")
patternC <- data.frame(expand.grid(x = seq(0, 1, length.out = sqrt(n)),
                                   y = seq(0, 1, length.out = sqrt(n))),
                       id="C")

patternD <- data.frame(expand.grid(x = rep(seq(0, 1, length.out = n/20),2),
                                   y = rep(seq(0, 1, length.out = n/20),2)),
                       id="D")
patternD$x <- jitter(patternD$x)
patternD$y <- jitter(patternD$y)

simDat <- bind_rows(patternA,patternB,patternC,patternD)

simDat <- simDat %>% group_by(id) %>% 
  mutate(x = scales::rescale(x), 
         y=scales::rescale(y))

ggplot(simDat,aes(x=x,y=y)) + 
  geom_point() + coord_fixed() + 
  facet_wrap(~id)


patternA <- simDat %>% filter(id=="A")
patternA <- ppp(x = patternA$x, y=patternA$y,
                xrange=c(0,1), yrange=c(0,1),
                unitname="km")
summary(patternA)

patternAEnv <- envelope(patternA, "Lest", verbose=FALSE)
plot(patternAEnv)

patternB <- simDat %>% filter(id=="B")
patternB <- ppp(x = patternB$x, y=patternB$y,
                xrange=c(0,1), yrange=c(0,1),
                unitname="km")
summary(patternB)

patternBEnv <- envelope(patternB, "Lest", verbose=FALSE)
plot(patternBEnv)

patternC <- simDat %>% filter(id=="C")
patternC <- ppp(x = patternC$x, y=patternC$y,
                xrange=c(0,1), yrange=c(0,1),
                unitname="km")
summary(patternC)

patternCEnv <- envelope(patternC, "Lest", verbose=FALSE)
plot(patternCEnv)

patternD <- simDat %>% filter(id=="D")
patternD <- ppp(x = patternD$x, y=patternD$y,
                xrange=c(0,1), yrange=c(0,1),
                unitname="km")
summary(patternD)

patternDEnv <- envelope(patternD, "Lest", verbose=FALSE)
plot(patternDEnv)



## ----------------------------------------------------------------------------
data(sporophores) # run ?sporophores to see the help file. It's cool.
summary(sporophores)
plot(sporophores, chars=c(16,1,2), cex=0.6, leg.args=list(cex=1.1))
points(0,0,pch=16, cex=2)
text(15,8,"Tree", cex=0.75)
levels(marks(sporophores))

## ----------------------------------------------------------------------------
# or directly in envelope (note i and j arguments)
cross_heb_pub <- envelope(sporophores, "Lcross", i="L pubescens", j = "Hebloma spp", verbose=FALSE)
plot(cross_heb_pub)

cross_pub_laccata <- envelope(sporophores, "Lcross", i="L pubescens", j = "L laccata", verbose=FALSE)
plot(cross_pub_laccata)

cross_heb_laccata <- envelope(sporophores, "Lcross", i="Hebloma spp", j = "L laccata", verbose=FALSE)
plot(cross_heb_laccata)

# ----echo=FALSE, eval=FALSE--------------------------------------------------
# Three levels, six comparisons. can do pairwise or use alltypes
# L laccata, L pubescens, Hebloma spp
sporophoresL <- alltypes(X = sporophores,fun = "Lcross",envelope=TRUE)
# so ugly
plot(sporophoresL)

