## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
library(tidyverse)
set.seed(42)


## -----------------------------------------------------------------------------
n <- 80
x <- rnorm(n, mean = 50, sd = 10)
y <- 2 + 0.6 * x + rnorm(n, sd = 6)

ggplot(data.frame(x, y), aes(x = x, y = y)) +
  geom_point(alpha = 0.5) + theme_minimal()



## -----------------------------------------------------------------------------
cor(x,y)
cor(y,x)


## -----------------------------------------------------------------------------
# Regression of y on x
fit_yx <- lm(y ~ x)
# Regression of x on y
fit_xy <- lm(x ~ y)


## -----------------------------------------------------------------------------
a <- coef(fit_xy)[1]
b <- coef(fit_xy)[2]

# Rearranged: y = -a/b + (1/b) * x
intercept_rearranged <- -a / b
slope_rearranged <- 1 / b


## ----echo=FALSE---------------------------------------------------------------
ggplot(data.frame(x, y), aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = coef(fit_yx)[1], slope = coef(fit_yx)[2],
              color = "steelblue", linewidth = 1, aes(linetype = "lm(y ~ x)")) +
  geom_abline(intercept = intercept_rearranged, slope = slope_rearranged,
              color = "firebrick", linewidth = 1, aes(linetype = "lm(x ~ y) rearranged")) +
  scale_linetype_manual(name = "Model",
                        values = c("lm(y ~ x)" = "solid",
                                   "lm(x ~ y) rearranged" = "dashed"),
                        guide = guide_legend()) +
  labs(title = "Two regressions, one scatter plot",
       subtitle = "Same data. Different lines.") +
  theme_minimal()


## -----------------------------------------------------------------------------
r <- cor(x, y)
slope_formula <- r * (sd(y) / sd(x))
slope_lm <- coef(fit_yx)[2]

round(c(from_formula = slope_formula, from_lm = slope_lm), 6)


## -----------------------------------------------------------------------------
# Small correlation, large n
r_small <- 0.11
n_large <- 1000
t_stat <- r_small * sqrt(n_large - 2) / sqrt(1 - r_small^2)
p_val <- 2 * pt(-abs(t_stat), df = n_large - 2)
p_val

