## ----echo=FALSE---------------------------------------------------------------
set.seed(184)


## ----message=FALSE------------------------------------------------------------
library(gstat)
library(fields)
library(ggplot2)
library(ggrepel)


## -----------------------------------------------------------------------------
dat <- data.frame(x = c(2,3,9,6,5), y = c(2,7,9,5,3), z = c(3,4,2,4,6))


## -----------------------------------------------------------------------------
dat2 <- rbind(dat, data.frame(x=5, y=5, z="?"))
dat2$lbls <- paste("(", dat2$x, ",", dat2$y, ") z=", dat2$z, sep="")
dat2$cols <- factor(c(rep("black",5), "purple"))

ggplot(data=dat2, mapping=aes(x=x, y=y, col=cols, label=lbls)) + 
  geom_point() +
  geom_label_repel() +
  lims(x=c(0,10), y=c(0,10)) +
  scale_color_identity() +
  coord_equal() +
  labs(subtitle = "Kriging: X marks the spot")


## -----------------------------------------------------------------------------
spher.covar <- function(dist, nugget, sill, range) {
   sigma <- sill + nugget
   covar <- matrix(0, nrow=nrow(dist), ncol=ncol(dist))
   for (i in 1:nrow(covar)) {
      for (j in 1:ncol(covar)) {
         if (dist[i,j] > range) 
            covar[i,j] <- 0
         else if (dist[i,j] == 0)
            covar[i,j] <- sigma
         else covar[i,j] <- sigma - (nugget + (sigma - nugget) * 
                                        (((3*dist[i,j])/(2*range)) - ((dist[i,j]^3)/(2*range^3))))
      }
   }
   return(covar)
}

exp.covar <- function(dist, nugget, sill, range) {
   sigma <- sill + nugget
   covar <- matrix(0, nrow=nrow(dist), ncol=ncol(dist))
   for (i in 1:nrow(covar)) {
      for (j in 1:ncol(covar)) {
         if (dist[i,j] == 0) {
            covar[i,j] <- sigma
         } else {
            covar[i,j] <- sill * exp(-dist[i,j]/range)
         }
      }
   }
   return(covar)
}

gauss.covar <- function(dist, nugget, sill, range) {
   sigma <- sill + nugget
   covar <- matrix(0, nrow=nrow(dist), ncol=ncol(dist))
   for (i in 1:nrow(covar)) {
      for (j in 1:ncol(covar)) {
         if (dist[i,j] == 0) {
            covar[i,j] <- sigma
         } else {
            covar[i,j] <- sill * exp(-dist[i,j]^2/range^2)
         }
      }
   }
   return(covar)
}


## -----------------------------------------------------------------------------
dmat <- round(rdist(cbind(dat$x, dat$y)), digits=4)               # distances among observed points
dvec <- round(rdist(cbind(dat$x, dat$y), cbind(5,5)), digits=4)   # distances to prediction location


## -----------------------------------------------------------------------------
vmod   <- "Gau"
nugget <- 2.5
sill   <- 7.5
range  <- 10

K <- gauss.covar(dmat, nugget, sill, range)  # covariance matrix among observed points
k <- gauss.covar(dvec, nugget, sill, range)  # covariance vector to prediction location


## -----------------------------------------------------------------------------
wts.sk <- solve(K) %*% k
est.sk  <- sum(wts.sk * (dat$z - mean(dat$z))) + mean(dat$z)
var.sk  <- (sill + nugget) - t(k) %*% solve(K) %*% k


## -----------------------------------------------------------------------------
K <- cbind(K, rep(1, nrow(K)))
K <- rbind(K, c(rep(1, ncol(K)-1), 0))
k <- c(k, 1)

wts.ok <- solve(K) %*% k
est.ok  <- sum(wts.ok[-length(wts.ok)] * dat$z)  # drop the Lagrange multiplier from the sum
var.ok  <- (sill + nugget) - t(k) %*% solve(K) %*% k


## -----------------------------------------------------------------------------
dat.vgram <- vgm(model=vmod, psill=sill, nugget=nugget, range=range)

est2.sk <- krige(z~1, ~x+y, data=dat, model=dat.vgram, 
                 newdata=data.frame(x=5,y=5), beta=mean(dat$z))
est2.ok <- krige(z~1, ~x+y, data=dat, model=dat.vgram, 
                 newdata=data.frame(x=5,y=5))

res.df <- data.frame(
  type = c('manual SK', 'gstat SK', 'manual OK', 'gstat OK'),
  est  = c(round(est.sk, 3), round(est2.sk$var1.pred, 3), 
           round(est.ok, 3), round(est2.ok$var1.pred, 3)),
  var  = c(round(var.sk, 3), round(est2.sk$var1.var, 3), 
           round(var.ok, 3), round(est2.ok$var1.var, 3))
)
res.df


## -----------------------------------------------------------------------------
dat3 <- dat2
dat3$lbls <- paste("(", dat3$x, ",", dat3$y, ") z=", dat3$z,
                   ", wt=", c(round(wts.sk, 2)[,1], ""), sep="")
dat3$lbls[6] <- paste("(5,5) z=", round(est.sk,3),
                      ", var=", round(sqrt(var.sk), 3), sep="")

ggplot(data=dat3, mapping=aes(x=x, y=y, col=cols, label=lbls)) + 
  geom_point() +
  geom_label_repel() +
  lims(x=c(0,10), y=c(0,10)) +
  scale_color_identity() +
  coord_equal() +
  labs(subtitle = paste("Simple Kriging; Model=", vmod, 
                        ", nug=", nugget, ", psill=", sill, 
                        ", range=", range, sep=""))


## -----------------------------------------------------------------------------
dat3$lbls[6] <- paste("(5,5) z=", round(est.ok,3),
                      ", var=", round(sqrt(var.ok), 3), sep="")

ggplot(data=dat3, mapping=aes(x=x, y=y, col=cols, label=lbls)) + 
  geom_point() +
  geom_label_repel() +
  lims(x=c(0,10), y=c(0,10)) +
  scale_color_identity() +
  coord_equal() +
  labs(subtitle = paste("Ordinary Kriging; Model=", vmod, 
                        ", nug=", nugget, ", psill=", sill, 
                        ", range=", range, sep=""))

