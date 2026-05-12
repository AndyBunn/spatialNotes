## -----------------------------------------------------------------------------
x <- mtcars$wt
y <- mtcars$mpg

beta1 <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
beta0 <- mean(y) - beta1 * mean(x)
c(intercept = beta0, slope = beta1)


## -----------------------------------------------------------------------------
# Same thing with lm()
coef(lm(mpg ~ wt, data = mtcars))


## -----------------------------------------------------------------------------
X <- as.matrix(mtcars$wt)
X <- cbind(1, X)   # prepend a column of ones for the intercept
y <- as.matrix(mtcars$mpg)


## -----------------------------------------------------------------------------
beta <- solve(t(X) %*% X) %*% t(X) %*% y
beta


## -----------------------------------------------------------------------------
# And again, verify with lm()
coef(lm(mpg ~ wt, data = mtcars))


## -----------------------------------------------------------------------------
# OLS (from above)
beta_ols <- solve(t(X) %*% X) %*% t(X) %*% y

# GLS with Sigma = identity matrix -- should match OLS exactly
Sigma <- diag(nrow(X))
Sigma_inv <- solve(Sigma)
beta_gls <- solve(t(X) %*% Sigma_inv %*% X) %*% t(X) %*% Sigma_inv %*% y

cbind(OLS = beta_ols, GLS = beta_gls)

