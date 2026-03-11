# Aside: Understanding Methods and Generics in R {- .aside-chapter}

## You’ve Come a Long Way in R

In your journey through the curriculum, you’ve learned a lot of R. You can wrangle data, fit models, make plots, write functions, etc. But you might also be starting to wonder:  **What’s actually going on under the hood?**

For example, how does `summary()` know what to do when you give it the output from `lm`? Or why does `plot()` sometimes give you a histogram and other times a scatterplot?

A good way to start answering those questions is to learn about **methods** and **generic functions**.

## Why Use Generics?

It’s nice when a function does the *right thing* depending on the kind of object you give it. You don’t want to memorize different function names for every kind of model—you just want to call `summary()`, `plot()`, or `print()` and get something useful. That’s what **generic functions** are for.

In R, many common functions like `summary()`, `plot()`, and `print()` are *generic*—they automatically call the version that matches the class of your object.

Let's run a linear model and use that as an example.


``` r
fit <- lm(mpg ~ wt, data = mtcars)
```

Ok. We have a new object called `fit` that is a `list` with a bunch of stuff in it that authors of the `lm` function think will be useful to you. Look at its structure to see what's in there.


``` r
str(fit)
```

```
## List of 12
##  $ coefficients : Named num [1:2] 37.29 -5.34
##   ..- attr(*, "names")= chr [1:2] "(Intercept)" "wt"
##  $ residuals    : Named num [1:32] -2.28 -0.92 -2.09 1.3 -0.2 ...
##   ..- attr(*, "names")= chr [1:32] "Mazda RX4" "Mazda RX4 Wag" "Datsun 710" "Hornet 4 Drive" ...
##  $ effects      : Named num [1:32] -113.65 -29.116 -1.661 1.631 0.111 ...
##   ..- attr(*, "names")= chr [1:32] "(Intercept)" "wt" "" "" ...
##  $ rank         : int 2
##  $ fitted.values: Named num [1:32] 23.3 21.9 24.9 20.1 18.9 ...
##   ..- attr(*, "names")= chr [1:32] "Mazda RX4" "Mazda RX4 Wag" "Datsun 710" "Hornet 4 Drive" ...
##  $ assign       : int [1:2] 0 1
##  $ qr           :List of 5
##   ..$ qr   : num [1:32, 1:2] -5.657 0.177 0.177 0.177 0.177 ...
##   .. ..- attr(*, "dimnames")=List of 2
##   .. .. ..$ : chr [1:32] "Mazda RX4" "Mazda RX4 Wag" "Datsun 710" "Hornet 4 Drive" ...
##   .. .. ..$ : chr [1:2] "(Intercept)" "wt"
##   .. ..- attr(*, "assign")= int [1:2] 0 1
##   ..$ qraux: num [1:2] 1.18 1.05
##   ..$ pivot: int [1:2] 1 2
##   ..$ tol  : num 1e-07
##   ..$ rank : int 2
##   ..- attr(*, "class")= chr "qr"
##  $ df.residual  : int 30
##  $ xlevels      : Named list()
##  $ call         : language lm(formula = mpg ~ wt, data = mtcars)
##  $ terms        :Classes 'terms', 'formula'  language mpg ~ wt
##   .. ..- attr(*, "variables")= language list(mpg, wt)
##   .. ..- attr(*, "factors")= int [1:2, 1] 0 1
##   .. .. ..- attr(*, "dimnames")=List of 2
##   .. .. .. ..$ : chr [1:2] "mpg" "wt"
##   .. .. .. ..$ : chr "wt"
##   .. ..- attr(*, "term.labels")= chr "wt"
##   .. ..- attr(*, "order")= int 1
##   .. ..- attr(*, "intercept")= int 1
##   .. ..- attr(*, "response")= int 1
##   .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
##   .. ..- attr(*, "predvars")= language list(mpg, wt)
##   .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
##   .. .. ..- attr(*, "names")= chr [1:2] "mpg" "wt"
##  $ model        :'data.frame':	32 obs. of  2 variables:
##   ..$ mpg: num [1:32] 21 21 22.8 21.4 18.7 18.1 14.3 24.4 22.8 19.2 ...
##   ..$ wt : num [1:32] 2.62 2.88 2.32 3.21 3.44 ...
##   ..- attr(*, "terms")=Classes 'terms', 'formula'  language mpg ~ wt
##   .. .. ..- attr(*, "variables")= language list(mpg, wt)
##   .. .. ..- attr(*, "factors")= int [1:2, 1] 0 1
##   .. .. .. ..- attr(*, "dimnames")=List of 2
##   .. .. .. .. ..$ : chr [1:2] "mpg" "wt"
##   .. .. .. .. ..$ : chr "wt"
##   .. .. ..- attr(*, "term.labels")= chr "wt"
##   .. .. ..- attr(*, "order")= int 1
##   .. .. ..- attr(*, "intercept")= int 1
##   .. .. ..- attr(*, "response")= int 1
##   .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
##   .. .. ..- attr(*, "predvars")= language list(mpg, wt)
##   .. .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
##   .. .. .. ..- attr(*, "names")= chr [1:2] "mpg" "wt"
##  - attr(*, "class")= chr "lm"
```

That list has 12 things in it including the coefficients, residuals, fitted values and a lot more. But getting the information you want in terms of interpreting the linear model isn't obvious. For instance, the object `fit` doesn’t directly contain the $R^2$ or p-values. Those are computed when you call `summary(fit)`.


``` r
summary(fit)
```

```
## 
## Call:
## lm(formula = mpg ~ wt, data = mtcars)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -4.5432 -2.3647 -0.1252  1.4096  6.8727 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  37.2851     1.8776  19.858  < 2e-16 ***
## wt           -5.3445     0.5591  -9.559 1.29e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 3.046 on 30 degrees of freedom
## Multiple R-squared:  0.7528,	Adjusted R-squared:  0.7446 
## F-statistic: 91.38 on 1 and 30 DF,  p-value: 1.294e-10
```

Ok. Check that out. We got a really nice formatted output that helps us interpret the linear model using the `summary` function. But note that the output of `summary(fit)` is really different than getting a summary of a numeric vector like the `mpg` in the `mtcars` data that went into the model.


``` r
summary(mtcars$mpg)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   10.40   15.43   19.20   20.09   22.80   33.90
```

### Same Function, Different Behavior — What Gives?

We used the same function (`summary`) and got two very different results.  
Answer? **Generics and methods.**

When we used `summary(fit)`, we're not calling some one-size-fits-all summary function. We're calling a **specific method** written for objects of class `"lm"`.

## How This Works

Here’s what R does when you call `summary(fit)`:

1. It checks the class of the object (`class(fit)`).
2. It looks for a function named `summary.classname`—in this case, `summary.lm`.
3. If it finds that method, it runs it.
4. If not, it falls back to `summary.default()`.

You can look under the hood:


``` r
class(fit)
```

```
## [1] "lm"
```

``` r
methods(summary)
```

```
##  [1] summary.aov                         summary.aovlist*                   
##  [3] summary.aspell*                     summary.check_packages_in_dir*     
##  [5] summary.connection                  summary.data.frame                 
##  [7] summary.Date                        summary.default                    
##  [9] summary.difftime                    summary.ecdf*                      
## [11] summary.factor                      summary.glm                        
## [13] summary.infl*                       summary.lm                         
## [15] summary.loess*                      summary.manova                     
## [17] summary.matrix                      summary.mlm*                       
## [19] summary.nls*                        summary.packageStatus*             
## [21] summary.POSIXct                     summary.POSIXlt                    
## [23] summary.ppr*                        summary.prcomp*                    
## [25] summary.princomp*                   summary.proc_time                  
## [27] summary.rlang_error*                summary.rlang_message*             
## [29] summary.rlang_trace*                summary.rlang_warning*             
## [31] summary.rlang:::list_of_conditions* summary.srcfile                    
## [33] summary.srcref                      summary.stepfun                    
## [35] summary.stl*                        summary.table                      
## [37] summary.tukeysmooth*                summary.warnings                   
## see '?methods' for accessing help and source code
```

``` r
summary.lm
```

```
## function (object, correlation = FALSE, symbolic.cor = FALSE, 
##     ...) 
## {
##     z <- object
##     p <- z$rank
##     rdf <- z$df.residual
##     if (p == 0) {
##         r <- z$residuals
##         n <- length(r)
##         w <- z$weights
##         if (is.null(w)) {
##             rss <- sum(r^2)
##         }
##         else {
##             rss <- sum(w * r^2)
##             r <- sqrt(w) * r
##         }
##         resvar <- rss/rdf
##         ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
##         class(ans) <- "summary.lm"
##         ans$aliased <- is.na(coef(object))
##         ans$residuals <- r
##         ans$df <- c(0L, n, length(ans$aliased))
##         ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL, 
##             c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
##         ans$sigma <- sqrt(resvar)
##         ans$r.squared <- ans$adj.r.squared <- 0
##         ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
##         if (correlation) 
##             ans$correlation <- ans$cov.unscaled
##         return(ans)
##     }
##     if (is.null(z$terms)) 
##         stop("invalid 'lm' object:  no 'terms' component")
##     if (!inherits(object, "lm")) 
##         warning("calling summary.lm(<fake-lm-object>) ...")
##     Qr <- qr.lm(object)
##     n <- NROW(Qr$qr)
##     if (is.na(z$df.residual) || n - p != z$df.residual) 
##         warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
##     r <- z$residuals
##     f <- z$fitted.values
##     if (!is.null(z$offset)) {
##         f <- f - z$offset
##     }
##     w <- z$weights
##     if (is.null(w)) {
##         mss <- if (attr(z$terms, "intercept")) 
##             sum((f - mean(f))^2)
##         else sum(f^2)
##         rss <- sum(r^2)
##     }
##     else {
##         mss <- if (attr(z$terms, "intercept")) {
##             m <- sum(w * f/sum(w))
##             sum(w * (f - m)^2)
##         }
##         else sum(w * f^2)
##         rss <- sum(w * r^2)
##         r <- sqrt(w) * r
##     }
##     resvar <- rss/rdf
##     if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) * 
##         1e-30) 
##         warning("essentially perfect fit: summary may be unreliable")
##     p1 <- 1L:p
##     R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
##     se <- sqrt(diag(R) * resvar)
##     est <- z$coefficients[Qr$pivot[p1]]
##     tval <- est/se
##     ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
##     ans$residuals <- r
##     ans$coefficients <- cbind(Estimate = est, `Std. Error` = se, 
##         `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf, 
##             lower.tail = FALSE))
##     ans$aliased <- is.na(z$coefficients)
##     ans$sigma <- sqrt(resvar)
##     ans$df <- c(p, rdf, NCOL(Qr$qr))
##     if (p != attr(z$terms, "intercept")) {
##         df.int <- if (attr(z$terms, "intercept")) 
##             1L
##         else 0L
##         ans$r.squared <- mss/(mss + rss)
##         ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
##             df.int)/rdf)
##         ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
##             numdf = p - df.int, dendf = rdf)
##     }
##     else ans$r.squared <- ans$adj.r.squared <- 0
##     ans$cov.unscaled <- R
##     dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
##         1)]
##     if (correlation) {
##         ans$correlation <- (R * resvar)/outer(se, se)
##         dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
##         ans$symbolic.cor <- symbolic.cor
##     }
##     if (!is.null(z$na.action)) 
##         ans$na.action <- z$na.action
##     class(ans) <- "summary.lm"
##     ans
## }
## <bytecode: 0x120863240>
## <environment: namespace:stats>
```

``` r
getS3method("summary", "lm")
```

```
## function (object, correlation = FALSE, symbolic.cor = FALSE, 
##     ...) 
## {
##     z <- object
##     p <- z$rank
##     rdf <- z$df.residual
##     if (p == 0) {
##         r <- z$residuals
##         n <- length(r)
##         w <- z$weights
##         if (is.null(w)) {
##             rss <- sum(r^2)
##         }
##         else {
##             rss <- sum(w * r^2)
##             r <- sqrt(w) * r
##         }
##         resvar <- rss/rdf
##         ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
##         class(ans) <- "summary.lm"
##         ans$aliased <- is.na(coef(object))
##         ans$residuals <- r
##         ans$df <- c(0L, n, length(ans$aliased))
##         ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL, 
##             c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
##         ans$sigma <- sqrt(resvar)
##         ans$r.squared <- ans$adj.r.squared <- 0
##         ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
##         if (correlation) 
##             ans$correlation <- ans$cov.unscaled
##         return(ans)
##     }
##     if (is.null(z$terms)) 
##         stop("invalid 'lm' object:  no 'terms' component")
##     if (!inherits(object, "lm")) 
##         warning("calling summary.lm(<fake-lm-object>) ...")
##     Qr <- qr.lm(object)
##     n <- NROW(Qr$qr)
##     if (is.na(z$df.residual) || n - p != z$df.residual) 
##         warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
##     r <- z$residuals
##     f <- z$fitted.values
##     if (!is.null(z$offset)) {
##         f <- f - z$offset
##     }
##     w <- z$weights
##     if (is.null(w)) {
##         mss <- if (attr(z$terms, "intercept")) 
##             sum((f - mean(f))^2)
##         else sum(f^2)
##         rss <- sum(r^2)
##     }
##     else {
##         mss <- if (attr(z$terms, "intercept")) {
##             m <- sum(w * f/sum(w))
##             sum(w * (f - m)^2)
##         }
##         else sum(w * f^2)
##         rss <- sum(w * r^2)
##         r <- sqrt(w) * r
##     }
##     resvar <- rss/rdf
##     if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) * 
##         1e-30) 
##         warning("essentially perfect fit: summary may be unreliable")
##     p1 <- 1L:p
##     R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
##     se <- sqrt(diag(R) * resvar)
##     est <- z$coefficients[Qr$pivot[p1]]
##     tval <- est/se
##     ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
##     ans$residuals <- r
##     ans$coefficients <- cbind(Estimate = est, `Std. Error` = se, 
##         `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf, 
##             lower.tail = FALSE))
##     ans$aliased <- is.na(z$coefficients)
##     ans$sigma <- sqrt(resvar)
##     ans$df <- c(p, rdf, NCOL(Qr$qr))
##     if (p != attr(z$terms, "intercept")) {
##         df.int <- if (attr(z$terms, "intercept")) 
##             1L
##         else 0L
##         ans$r.squared <- mss/(mss + rss)
##         ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
##             df.int)/rdf)
##         ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
##             numdf = p - df.int, dendf = rdf)
##     }
##     else ans$r.squared <- ans$adj.r.squared <- 0
##     ans$cov.unscaled <- R
##     dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
##         1)]
##     if (correlation) {
##         ans$correlation <- (R * resvar)/outer(se, se)
##         dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
##         ans$symbolic.cor <- symbolic.cor
##     }
##     if (!is.null(z$na.action)) 
##         ans$na.action <- z$na.action
##     class(ans) <- "summary.lm"
##     ans
## }
## <bytecode: 0x120863240>
## <environment: namespace:stats>
```

## What’s the Difference Between a Generic and a Method?

- A **generic function** is the function you call—like `summary()` or `plot()`.  
  It doesn’t actually do the work itself. Instead, it decides **which method to run** based on the class of the object.

- A **method** is the specific function written for a particular class.  
  For example:
  - `summary.lm()` handles linear models.
  - `summary.data.frame()` handles data frames.
  - `summary.default()` is the fallback if no class-specific method is found.

You can look at the generic directly:


``` r
summary
```

```
## function (object, ...) 
## UseMethod("summary")
## <bytecode: 0x12086cb80>
## <environment: namespace:base>
```

This call to `UseMethod("summary")` tells R to dispatch to the right method depending on the class of the object.

---

## Common Generic Functions

Here are some other generic functions you're likely to run into:

- `print()` – display a human-readable version  
- `summary()` – provide a richer summary  
- `plot()` – show plots  
- `predict()` – make predictions  
- `residuals()` – return residuals  
- `fitted()` – return fitted values  
- `coef()` – extract coefficients  
- `update()` – refit a model with changes  
- `AIC()` / `BIC()` – compare model fit  
- `anova()` – compare nested models  
- `terms()` / `formula()` – extract parts of model structure  

## How to Explore and Understand Methods

Want to know what’s going on behind a function call? Here are a few tools:


``` r
methods(summary)
```

```
##  [1] summary.aov                         summary.aovlist*                   
##  [3] summary.aspell*                     summary.check_packages_in_dir*     
##  [5] summary.connection                  summary.data.frame                 
##  [7] summary.Date                        summary.default                    
##  [9] summary.difftime                    summary.ecdf*                      
## [11] summary.factor                      summary.glm                        
## [13] summary.infl*                       summary.lm                         
## [15] summary.loess*                      summary.manova                     
## [17] summary.matrix                      summary.mlm*                       
## [19] summary.nls*                        summary.packageStatus*             
## [21] summary.POSIXct                     summary.POSIXlt                    
## [23] summary.ppr*                        summary.prcomp*                    
## [25] summary.princomp*                   summary.proc_time                  
## [27] summary.rlang_error*                summary.rlang_message*             
## [29] summary.rlang_trace*                summary.rlang_warning*             
## [31] summary.rlang:::list_of_conditions* summary.srcfile                    
## [33] summary.srcref                      summary.stepfun                    
## [35] summary.stl*                        summary.table                      
## [37] summary.tukeysmooth*                summary.warnings                   
## see '?methods' for accessing help and source code
```

``` r
methods(plot)
```

```
##  [1] plot.acf*           plot.data.frame*    plot.decomposed.ts*
##  [4] plot.default        plot.dendrogram*    plot.density*      
##  [7] plot.ecdf           plot.factor*        plot.formula*      
## [10] plot.function       plot.hclust*        plot.histogram*    
## [13] plot.HoltWinters*   plot.isoreg*        plot.lm*           
## [16] plot.medpolish*     plot.mlm*           plot.ppr*          
## [19] plot.prcomp*        plot.princomp*      plot.profile*      
## [22] plot.profile.nls*   plot.R6*            plot.raster*       
## [25] plot.spec*          plot.stepfun        plot.stl*          
## [28] plot.table*         plot.ts             plot.tskernel*     
## [31] plot.TukeyHSD*     
## see '?methods' for accessing help and source code
```

``` r
methods(class = "lm")
```

```
##  [1] add1           alias          anova          case.names     coerce        
##  [6] confint        cooks.distance deviance       dfbeta         dfbetas       
## [11] drop1          dummy.coef     effects        extractAIC     family        
## [16] formula        hatvalues      influence      initialize     kappa         
## [21] labels         logLik         model.frame    model.matrix   nobs          
## [26] plot           predict        print          proj           qr            
## [31] residuals      rstandard      rstudent       show           simulate      
## [36] slotsFromS3    summary        variable.names vcov          
## see '?methods' for accessing help and source code
```

``` r
getS3method("summary", "lm")
```

```
## function (object, correlation = FALSE, symbolic.cor = FALSE, 
##     ...) 
## {
##     z <- object
##     p <- z$rank
##     rdf <- z$df.residual
##     if (p == 0) {
##         r <- z$residuals
##         n <- length(r)
##         w <- z$weights
##         if (is.null(w)) {
##             rss <- sum(r^2)
##         }
##         else {
##             rss <- sum(w * r^2)
##             r <- sqrt(w) * r
##         }
##         resvar <- rss/rdf
##         ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
##         class(ans) <- "summary.lm"
##         ans$aliased <- is.na(coef(object))
##         ans$residuals <- r
##         ans$df <- c(0L, n, length(ans$aliased))
##         ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL, 
##             c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
##         ans$sigma <- sqrt(resvar)
##         ans$r.squared <- ans$adj.r.squared <- 0
##         ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
##         if (correlation) 
##             ans$correlation <- ans$cov.unscaled
##         return(ans)
##     }
##     if (is.null(z$terms)) 
##         stop("invalid 'lm' object:  no 'terms' component")
##     if (!inherits(object, "lm")) 
##         warning("calling summary.lm(<fake-lm-object>) ...")
##     Qr <- qr.lm(object)
##     n <- NROW(Qr$qr)
##     if (is.na(z$df.residual) || n - p != z$df.residual) 
##         warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
##     r <- z$residuals
##     f <- z$fitted.values
##     if (!is.null(z$offset)) {
##         f <- f - z$offset
##     }
##     w <- z$weights
##     if (is.null(w)) {
##         mss <- if (attr(z$terms, "intercept")) 
##             sum((f - mean(f))^2)
##         else sum(f^2)
##         rss <- sum(r^2)
##     }
##     else {
##         mss <- if (attr(z$terms, "intercept")) {
##             m <- sum(w * f/sum(w))
##             sum(w * (f - m)^2)
##         }
##         else sum(w * f^2)
##         rss <- sum(w * r^2)
##         r <- sqrt(w) * r
##     }
##     resvar <- rss/rdf
##     if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) * 
##         1e-30) 
##         warning("essentially perfect fit: summary may be unreliable")
##     p1 <- 1L:p
##     R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
##     se <- sqrt(diag(R) * resvar)
##     est <- z$coefficients[Qr$pivot[p1]]
##     tval <- est/se
##     ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
##     ans$residuals <- r
##     ans$coefficients <- cbind(Estimate = est, `Std. Error` = se, 
##         `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf, 
##             lower.tail = FALSE))
##     ans$aliased <- is.na(z$coefficients)
##     ans$sigma <- sqrt(resvar)
##     ans$df <- c(p, rdf, NCOL(Qr$qr))
##     if (p != attr(z$terms, "intercept")) {
##         df.int <- if (attr(z$terms, "intercept")) 
##             1L
##         else 0L
##         ans$r.squared <- mss/(mss + rss)
##         ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
##             df.int)/rdf)
##         ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
##             numdf = p - df.int, dendf = rdf)
##     }
##     else ans$r.squared <- ans$adj.r.squared <- 0
##     ans$cov.unscaled <- R
##     dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
##         1)]
##     if (correlation) {
##         ans$correlation <- (R * resvar)/outer(se, se)
##         dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
##         ans$symbolic.cor <- symbolic.cor
##     }
##     if (!is.null(z$na.action)) 
##         ans$na.action <- z$na.action
##     class(ans) <- "summary.lm"
##     ans
## }
## <bytecode: 0x120863240>
## <environment: namespace:stats>
```

``` r
class(fit)
```

```
## [1] "lm"
```

## Why This Matters

Understanding how generics and methods work helps you read other people’s code, write better functions, and make sense of R’s behavior when things feel a little magic.
