## ----echo=FALSE---------------------------------------------------------------
set.seed(2613)


## ----message=FALSE------------------------------------------------------------
library(tidyverse)
library(sf)
library(gstat)
library(PNWColors)


## ----echo=FALSE---------------------------------------------------------------
anem7 <- pnw_palette(name="Anemone", n=7, type="discrete")


## -----------------------------------------------------------------------------
n   <- 50
x   <- runif(n, min=-4, max=4)
eps <- rnorm(n, sd=15)
y   <- 30 + 2*x^3 - 5*x + eps
dat <- data.frame(x=x, y=y)


## -----------------------------------------------------------------------------
ggplot(dat, aes(x=x, y=y)) +
  geom_point(size=2.5, alpha=0.8) +
  theme_minimal()


## -----------------------------------------------------------------------------
test_ov_idx <- sample(1:nrow(dat), size=round(nrow(dat) * 0.33))
train_ov    <- dat[-test_ov_idx, ]
test_ov     <- dat[test_ov_idx,  ]

degrees <- 1:10
rsq_in  <- numeric(length(degrees))
rsq_out <- numeric(length(degrees))

for(d in degrees){
  fit_d      <- lm(y ~ poly(x, d), data=train_ov)
  rsq_in[d]  <- cor(train_ov$y, predict(fit_d))^2
  rsq_out[d] <- cor(test_ov$y,  predict(fit_d, newdata=test_ov))^2
}

bv <- data.frame(
  degree = rep(degrees, 2),
  rsq    = c(rsq_in, rsq_out),
  set    = rep(c("Training", "Test"), each=length(degrees))
)

ggplot(bv, aes(x=degree, y=rsq, color=set)) +
  geom_line(linewidth=1.1) +
  geom_point(size=2.5) +
  scale_x_continuous(breaks=degrees) +
  scale_color_manual(values=c(anem7[7], anem7[1])) +
  labs(x="Polynomial degree", y=bquote(R^2), color=NULL) +
  theme_minimal()


## -----------------------------------------------------------------------------
test_idx            <- sample(1:nrow(dat), size=round(nrow(dat) * 0.25))
dat$split           <- "train"
dat$split[test_idx] <- "test"

ggplot(dat, aes(x=x, y=y, color=split)) +
  geom_point(size=3) +
  scale_color_manual(values=c(anem7[1], anem7[7])) +
  theme_minimal()


## -----------------------------------------------------------------------------
train     <- dat[dat$split == "train", ]
test      <- dat[dat$split == "test",  ]

fit       <- lm(y ~ poly(x, 3), data=train)
test$yhat <- predict(fit, newdata=test)

rsq  <- cor(test$y, test$yhat)^2
rmse <- sqrt(mean((test$yhat - test$y)^2))
mae  <- mean(abs(test$yhat - test$y))

c("rsq"=rsq,"rmse"=rmse,"mae"=mae)


## -----------------------------------------------------------------------------
ggplot(test, aes(x=y, y=yhat)) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_point(size=3, color=anem7[7]) +
  coord_fixed(xlim=range(test$y, test$yhat),
              ylim=range(test$y, test$yhat)) +
  labs(x="Observed", y="Predicted") +
  theme_minimal()


## -----------------------------------------------------------------------------
rmseNull <- sqrt(mean((mean(train$y) - test$y)^2))
maeNull  <- mean(abs(mean(train$y) - test$y))

rmseNull
1 - (rmse / rmseNull)


## -----------------------------------------------------------------------------
k    <- 10
dat2 <- dat[sample(nrow(dat)), c("x","y")]  # shuffle, drop the split column
dat2$fold <- cut(seq(1, nrow(dat2)), breaks=k, labels=FALSE)

results <- data.frame(fold=1:k, rsq=NA, rmse=NA, mae=NA)

for(i in 1:k){
  train_i <- dat2[dat2$fold != i, ]
  test_i  <- dat2[dat2$fold == i, ]
  
  fit_i  <- lm(y ~ poly(x, 3), data=train_i)
  yhat_i <- predict(fit_i, newdata=test_i)
  
  results$rsq[i]  <- cor(test_i$y, yhat_i)^2
  results$rmse[i] <- sqrt(mean((yhat_i - test_i$y)^2))
  results$mae[i]  <- mean(abs(yhat_i - test_i$y))
}

results


## -----------------------------------------------------------------------------
colMeans(results[, c("rsq","rmse","mae")])


## -----------------------------------------------------------------------------
results %>%
  pivot_longer(cols=c(rsq, rmse, mae),
               names_to="metric", values_to="value") %>%
  ggplot(aes(x=value)) +
  geom_histogram(bins=6, fill=anem7[5], color=anem7[7]) +
  labs(x=NULL, y=NULL) +
  facet_wrap(~metric, scales="free_x") +
  theme_minimal()


## ----message=FALSE------------------------------------------------------------
data(meuse.all)
meuse_sf <- st_as_sf(meuse.all, coords=c("x","y"), crs=28992)

test_m_idx                 <- sample(1:nrow(meuse_sf), size=round(nrow(meuse_sf) * 0.25))
meuse_sf$split             <- "train"
meuse_sf$split[test_m_idx] <- "test"

ggplot(meuse_sf) +
  geom_sf(aes(color=split, shape=split), size=3) +
  scale_color_manual(values=c(anem7[1], anem7[7])) +
  theme_minimal() +
  labs(title="Random 75/25 holdout",
       subtitle="Test points scattered among training points")

