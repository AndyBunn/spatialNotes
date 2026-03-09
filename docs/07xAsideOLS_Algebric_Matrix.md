# Aside: OLS via Algebra and Matrices {- .aside-chapter}

## Introduction

Ordinary Least Squares (OLS) regression is a method for estimating the parameters in a linear regression model. This document explains OLS regression first using algebra and then using matrix notation, with a practical implementation in R.

## Algebraic Approach to OLS

In simple linear regression, we model the relationship between a predictor variable $X$ and a response variable $y$ using the equation:

$$ y_i = \beta_0 + \beta_1 X_i + \epsilon_i $$

where:
- $y_i$ is the response variable for observation $i$,
- $X_i$ is the predictor variable for observation $i$,
- $\beta_0$ is the intercept,
- $\beta_1$ is the slope,
- $\epsilon_i$ is the error term.

### Goal: Minimize the Sum of Squared Errors

The OLS method aims to find the values of $\beta_0$ and $\beta_1$ that minimize the sum of squared errors (SSE):

$$SSE = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2$$

where $\hat{y}_i = \beta_0 + \beta_1 X_i$ is the predicted value of $y_i$.

### Deriving the Coefficients

To minimize the SSE, we take the partial derivatives of the SSE with respect to $\beta_0$ and $\beta_1$ and set them to zero:

$$\frac{\partial SSE}{\partial \beta_0} = -2 \sum_{i=1}^{n} (y_i - \beta_0 - \beta_1 X_i) = 0$$
$$\frac{\partial SSE}{\partial \beta_1} = -2 \sum_{i=1}^{n} (y_i - \beta_0 - \beta_1 X_i) X_i = 0$$

Solving these equations gives the OLS estimates:

$$\hat{\beta_1} = \frac{\sum_{i=1}^{n} (X_i - \bar{X})(y_i - \bar{y})}{\sum_{i=1}^{n} (X_i - \bar{X})^2}$$
$$\hat{\beta_0} = \bar{y} - \hat{\beta_1} \bar{X}$$

where $\bar{X}$ and $\bar{y}$ are the means of $X$ and $y$, respectively.

## Matrix Approach to OLS

The algebraic method works well for simple linear regression, but for multiple regression, the matrix approach is more convenient. In matrix notation, the linear model is:

$$\mathbf{y} = \mathbf{X} \boldsymbol{\beta} + \boldsymbol{\epsilon}$$

where:

* $\mathbf{y}$ is an $n \times 1$ vector of the response variable,
* $\mathbf{X}$ is an $n \times p$ matrix of predictors (with a column of ones for the intercept),
* $\boldsymbol{\beta}$ is a $p \times 1$ vector of coefficients,
* $\boldsymbol{\epsilon}$ is an $n \times 1$ vector of errors.

### Goal: Minimize the Sum of Squared Errors

The OLS method minimizes the sum of squared errors, which in matrix form is:

$$SSE = (\mathbf{y} - \mathbf{X} \boldsymbol{\beta})^T (\mathbf{y} - \mathbf{X} \boldsymbol{\beta})$$

### Deriving the Coefficients

To find the coefficients that minimize the SSE, we set the derivative with respect to $\boldsymbol{\beta}$ to zero:

$$\frac{\partial SSE}{\partial \boldsymbol{\beta}} = -2 \mathbf{X}^T (\mathbf{y} - \mathbf{X} \boldsymbol{\beta}) = 0$$

Solving for $\boldsymbol{\beta}$ gives the OLS estimate:

$$\hat{\boldsymbol{\beta}} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}$$

## Practical Implementation in R

### Step 1: Prepare the Data

We'll use a small dataset with one predictor (X) and one response variable (y).


``` r
# Create a simple dataset
set.seed(123)  # For reproducibility
X <- matrix(rnorm(10), ncol=1)  # 10 random predictors
y <- 3 + 2 * X + rnorm(10)  # Response variable with some noise

# Display the data
data.frame(X = X, y = y)
```

```
##              X         y
## 1  -0.56047565  3.103131
## 2  -0.23017749  2.899459
## 3   1.55870831  6.518188
## 4   0.07050839  3.251699
## 5   0.12928774  2.702734
## 6   1.71506499  8.217043
## 7   0.46091621  4.419683
## 8  -1.26506123 -1.496740
## 9  -0.68685285  2.327650
## 10 -0.44566197  1.635885
```

### Step 2: Add an Intercept Column to the Predictor Matrix

We need to add a column of ones to the predictor matrix to account for the intercept term in the regression model.


``` r
# Add intercept column to the predictor matrix
X_with_intercept <- cbind(1, X)
# Display the new predictor matrix
X_with_intercept
```

```
##       [,1]        [,2]
##  [1,]    1 -0.56047565
##  [2,]    1 -0.23017749
##  [3,]    1  1.55870831
##  [4,]    1  0.07050839
##  [5,]    1  0.12928774
##  [6,]    1  1.71506499
##  [7,]    1  0.46091621
##  [8,]    1 -1.26506123
##  [9,]    1 -0.68685285
## [10,]    1 -0.44566197
```

### Step 3: Compute the OLS Coefficients

The OLS coefficients can be computed using the formula $\hat{\boldsymbol{\beta}} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}$.


``` r
# Compute OLS coefficients
beta <- solve(t(X_with_intercept) %*% X_with_intercept) %*% t(X_with_intercept) %*% y

# Print the coefficients
beta
```

```
##          [,1]
## [1,] 3.161708
## [2,] 2.628661
```

### Step 4: Make Predictions

We can use the computed coefficients to make predictions for the response variable.


``` r
# Make predictions
y_pred <- X_with_intercept %*% beta

# Print predictions
y_pred
```

```
##             [,1]
##  [1,]  1.6884072
##  [2,]  2.5566491
##  [3,]  7.2590236
##  [4,]  3.3470504
##  [5,]  3.5015614
##  [6,]  7.6700323
##  [7,]  4.3733002
##  [8,] -0.1637095
##  [9,]  1.3562044
## [10,]  1.9902135
```

### Step 5: Verify with `lm` Function (Optional)

To verify our matrix-based implementation, we can compare the results with the `lm` function in R.


``` r
# Verify with lm function
model <- lm(y ~ X)
summary(model)
```

```
## 
## Call:
## lm(formula = y ~ X)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.33303 -0.64421 -0.02448  0.49596  1.41472 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   3.1617     0.2852  11.086 3.91e-06 ***
## X             2.6287     0.3141   8.368 3.15e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.8988 on 8 degrees of freedom
## Multiple R-squared:  0.8975,	Adjusted R-squared:  0.8847 
## F-statistic: 70.03 on 1 and 8 DF,  p-value: 3.153e-05
```

## Conclusion

In this document, we explained OLS regression using both algebra and matrix notation. We demonstrated how to derive the OLS estimates and implemented the matrix-based method in R. This approach helps build a strong understanding of the underlying mechanics of linear regression.
