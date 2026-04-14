## ----echo=FALSE, include=FALSE------------------------------------------------
set.seed(1984)
knitr::purl("03Point_Patterns.qmd", output = "03Point_Patterns.R", documentation = 1)


## ----message=FALSE------------------------------------------------------------
library(tidyverse)
library(spatstat)


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
data(japanesepines)
data(redwood)
data(bei)

summary(japanesepines)
summary(redwood)
summary(bei)

plot(japanesepines)
plot(redwood)
plot(bei)


## -----------------------------------------------------------------------------
tibble(x=japanesepines$x * japanesepines$window$units$multiplier, 
           y=japanesepines$y * japanesepines$window$units$multiplier) %>%
  ggplot(mapping = aes(x=x,y=y)) +
  geom_point() +
  labs(x=paste("x (",
               japanesepines$window$units$plural,
               ")",sep=""),
       y=paste("y (",
               japanesepines$window$units$plural,
               ")",sep="")) +
  coord_fixed() + theme_minimal()


## -----------------------------------------------------------------------------
x <- rnorm(n)
hist(x,freq=FALSE)
lines(density(x),col="red")


## -----------------------------------------------------------------------------
ggplot() + 
  geom_histogram(mapping = aes(x=x,after_stat(density)),fill="grey") + 
  geom_density(mapping = aes(x=x))


## -----------------------------------------------------------------------------
jp.den <- density(japanesepines)
summary(jp.den)


## -----------------------------------------------------------------------------
persp(jp.den,theta = 30, phi = 30)


## -----------------------------------------------------------------------------
plot(jp.den)
contour(jp.den,add=TRUE)
points(japanesepines,pch=20)


## -----------------------------------------------------------------------------
plot(density(redwood)) # points(redwood,pch=20)
plot(density(bei)) # points(bei,pch=20)


## -----------------------------------------------------------------------------
n <- 100
japanesepinesK <- envelope(japanesepines, fun = Kest, nsim = n, verbose=FALSE)
beiK <- envelope(bei, fun = Kest, nsim = n, verbose=FALSE)
redwoodK <- envelope(redwood, fun = Kest, nsim = n, verbose=FALSE)


## -----------------------------------------------------------------------------
plot(japanesepinesK)


## -----------------------------------------------------------------------------
plot(redwoodK)


## -----------------------------------------------------------------------------
plot(beiK)


## -----------------------------------------------------------------------------
plot(envelope(redwood, fun = Kest, nsim = 10,
              verbose=FALSE, rmax=0.5),main="")


## -----------------------------------------------------------------------------
japanesepinesK2 <- Kest(japanesepines, verbose=FALSE, 
                       correction=c("none","isotropic","border"))
plot(japanesepinesK2)


## -----------------------------------------------------------------------------
redwoodL <- envelope(redwood, fun = Lest, nsim = n, verbose=FALSE)
plot(redwoodL)


## -----------------------------------------------------------------------------
ggplot(redwoodK, mapping = aes(x=r, ymin = lo-pi*r^2, ymax=hi-pi*r^2)) +
  geom_ribbon(fill="grey",alpha=0.5) + 
  geom_line(mapping = aes(y=theo-pi*r^2),col="red", linetype="dashed") +
  geom_line(mapping = aes(y=obs-pi*r^2)) +
  labs(y=expression(K(r) - pi~r^2), x = "r") +
  scale_x_continuous(expand = c(0,0))


## -----------------------------------------------------------------------------
data(longleaf)
summary(longleaf)
plot(longleaf)


## -----------------------------------------------------------------------------
hist(longleaf$marks)


## -----------------------------------------------------------------------------
bigTrees <- subset(longleaf, marks > 50)
plot(envelope(bigTrees, fun = Lest, nsim = n, verbose=FALSE),
     main = "Big Trees")


## -----------------------------------------------------------------------------
longleaf2 <- cut(longleaf, breaks=c(0,30,Inf), 
                 labels=c("Sapling","Adult"))
summary(longleaf2)
plot(longleaf2)


## -----------------------------------------------------------------------------
adults <- subset(longleaf2, marks == "Adult",drop=TRUE)
plot(envelope(adults, fun = Lest, verbose=FALSE), main="Adult")

saplings <- subset(longleaf2, marks == "Sapling",drop=TRUE)
plot(envelope(saplings, fun = Lest, verbose=FALSE), main="Sapling")


## -----------------------------------------------------------------------------
adults.den <- density(adults)
plot(adults.den)
points(saplings, pch=20)


## -----------------------------------------------------------------------------
longleaf2L <- envelope(longleaf2, "Lcross", 
                       i ="Sapling",j = "Adult",
                       verbose=FALSE)
plot(longleaf2L)


## ----eval=FALSE---------------------------------------------------------------
# foo <- clickppp(n = 50) # 50 points
# # KDE
# plot(density(foo))
# points(foo,pch=20)
# # L
# env <- envelope(foo, Lest, nsim = 100)
# plot(env)


## ----eval=FALSE---------------------------------------------------------------
# demo(data)


## ----eval=F-------------------------------------------------------------------
# patternD <- simDat %>% filter(id=="D")
# patternD <- ppp(x = patternD$x, y=patternD$y,
#                 xrange=c(0,1), yrange=c(0,1),
#                 unitname="km")
# summary(patternD)
# plot(patternD)


## -----------------------------------------------------------------------------
data(sporophores) # run ?sporophores to see the help file. It's cool.
summary(sporophores)
plot(sporophores, chars=c(16,1,2), cex=0.6, leg.args=list(cex=1.1))
points(0,0,pch=16, cex=2)
text(15,8,"Tree", cex=0.75)


## -----------------------------------------------------------------------------
# using subset
heb_pub <- subset(sporophores, 
                  marks %in% c("L pubescens","Hebloma spp"), drop=TRUE)
cross_heb_pub <- envelope(heb_pub, "Lcross", verbose=FALSE)
# or directly in envelope (note i and j arguments)
cross_heb_pub <- envelope(sporophores, "Lcross", i="L pubescens", j = "Hebloma spp", verbose=FALSE)


## ----echo=FALSE, eval=FALSE---------------------------------------------------
# # Three levels, six comparisons. can do pairwise or use alltypes
# # L laccata, L pubescens, Hebloma spp
# sporophoresL <- alltypes(X = sporophores,fun = "Lcross",envelope=TRUE)
# # so ugly
# plot(sporophoresL)
# 
# # vs going pairwise e.g.,
# # L pubescens, Hebloma spp
# heb_pub <- subset(sporophores, marks %in% c("L pubescens","Hebloma spp"), drop=TRUE)
# cross_heb_pub <- envelope(heb_pub, "Lcross", verbose=FALSE)
# plot(cross_heb_pub)
# # can use i and j args to envelope as well
# cross_heb_pub <- envelope(sporophores, "Lcross", i="L pubescens", j = "Hebloma spp", verbose=FALSE)
# plot(cross_heb_pub)


## -----------------------------------------------------------------------------
n <- 12
area <- 100 # in m2
x <- runif(n = n, min = 0,max = 10)
y <- runif(n = n, min = 0,max = 10)
ggplot() + 
  geom_point(aes(x=x,y=y),size=3) + 
  coord_fixed()


## -----------------------------------------------------------------------------
Dmat <- as.matrix(dist(cbind(x,y)))
Dmat


## -----------------------------------------------------------------------------
Dvec <- c(Dmat[upper.tri(Dmat)],Dmat[lower.tri(Dmat)])
Dvec


## -----------------------------------------------------------------------------
r <- seq(0,2.5,by=0.1)


## -----------------------------------------------------------------------------
Kr <- numeric()
for(i in 1:length(r)){
  Kr[i] = (area/(n * (n-1))) * sum(Dvec <= r[i])
}
p1 <- ggplot() +
  geom_line(aes(x=r,y=pi*r^2),color="red", linetype="dashed") +
  geom_line(aes(x=r,y=Kr),color="blue",linewidth=1) +
  labs(y="K(r)",x="r")
p1


## -----------------------------------------------------------------------------
xy <- as.ppp(cbind(x,y),W=c(0,10,0,10))
xy
xyK <- Kest(xy,r=r,correction = "none")

# match 
p1 + geom_line(data=xyK,mapping = aes(x=r,y=un),linetype="dashed",color="white")


