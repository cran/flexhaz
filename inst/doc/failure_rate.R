## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(flexhaz)

## -----------------------------------------------------------------------------
# Exponential: constant hazard h(t) = lambda
exp_dist <- dfr_exponential(lambda = 0.5)

# Weibull: power-law hazard h(t) = (k/sigma)(t/sigma)^(k-1)
weib_dist <- dfr_weibull(shape = 2, scale = 3)

# Gompertz: exponentially increasing hazard h(t) = a*exp(b*t)
gomp_dist <- dfr_gompertz(a = 0.01, b = 0.1)

# Log-logistic: non-monotonic hazard (increases then decreases)
ll_dist <- dfr_loglogistic(alpha = 10, beta = 2)

print(exp_dist)
is_dfr_dist(exp_dist)

## -----------------------------------------------------------------------------
# Custom: linear increasing hazard h(t) = a + b*t
linear_dist <- dfr_dist(
  rate = function(t, par, ...) {
    a <- par[1]
    b <- par[2]
    a + b * t
  },
  par = c(a = 0.1, b = 0.01)
)

## -----------------------------------------------------------------------------
h <- hazard(exp_dist)
H <- cum_haz(exp_dist)

# Evaluate at specific times
h(1)      # Hazard at t=1
h(c(1, 2, 3))  # Vectorized

H(2)      # Cumulative hazard at t=2

## -----------------------------------------------------------------------------
S <- surv(exp_dist)
F <- cdf(exp_dist)

# Verify S(t) + F(t) = 1
t <- 2
c(survival = S(t), cdf = F(t), sum = S(t) + F(t))

## -----------------------------------------------------------------------------
# Use stats::density which flexhaz implements as density.dfr_dist
pdf_fn <- density(exp_dist)

# For exponential: f(t) = lambda * exp(-lambda * t)
t <- 1
lambda <- 0.5
c(computed = pdf_fn(t), analytical = lambda * exp(-lambda * t))

## -----------------------------------------------------------------------------
Q <- inv_cdf(exp_dist)

# Median of exponential is log(2)/lambda
median_computed <- Q(0.5)
median_analytical <- log(2) / 0.5
c(computed = median_computed, analytical = median_analytical)

## -----------------------------------------------------------------------------
samp <- sampler(exp_dist)

set.seed(42)
samples <- samp(1000)

# Compare sample mean to theoretical mean (1/lambda = 2)
c(sample_mean = mean(samples), theoretical = 1/0.5)

## -----------------------------------------------------------------------------
h <- hazard(exp_dist)

# Use default lambda = 0.5
h(1)

# Override with lambda = 2
h(1, par = c(2))

## -----------------------------------------------------------------------------
# Simulate exact failure times from exponential(lambda=1)
set.seed(123)
true_lambda <- 1
n <- 50
times <- rexp(n, rate = true_lambda)

# Create data frame with standard survival format
# delta = 1 means exact observation, delta = 0 means censored
df_exact <- data.frame(t = times, delta = rep(1, n))
head(df_exact)

## -----------------------------------------------------------------------------
dist <- dfr_dist(
  rate = function(t, par, ...) rep(par[1], length(t)),
  par = NULL  # No default - must be supplied
)

ll <- loglik(dist)

# Evaluate at different parameter values
ll(df_exact, par = c(0.5))  # lambda = 0.5
ll(df_exact, par = c(1.0))  # lambda = 1.0 (true value)
ll(df_exact, par = c(2.0))  # lambda = 2.0

## -----------------------------------------------------------------------------
s <- score(dist)
s(df_exact, par = c(1.0))  # Should be close to 0 at MLE

## -----------------------------------------------------------------------------
H_ll <- hess_loglik(dist)
hess <- H_ll(df_exact, par = c(1.0))
hess  # Should be negative (concave at maximum)

## ----warning=FALSE------------------------------------------------------------
solver <- fit(dist)

# Find MLE starting from initial guess
result <- solver(df_exact, par = c(0.5), method = "BFGS")

# Extract fitted parameters (the fisher_mle class from likelihood.model uses coef())
coef(result)

# Compare to analytical MLE: lambda_hat = n / sum(t)
analytical_mle <- n / sum(times)
c(fitted = coef(result), analytical = analytical_mle, true = true_lambda)

## ----warning=FALSE------------------------------------------------------------
# Some observations are censored (patient still alive at study end)
df_mixed <- data.frame(
  t = c(1, 2, 3, 4, 5, 6, 7, 8),
  delta = c(1, 1, 1, 0, 0, 1, 1, 0)  # 0 = censored
)

ll <- loglik(dist)

# Fit with censored data
solver <- fit(dist)
result <- solver(df_mixed, par = c(0.5), method = "BFGS")
coef(result)

## ----warning=FALSE------------------------------------------------------------
# Create Weibull DFR
weibull <- dfr_dist(
  rate = function(t, par, ...) {
    k <- par[1]
    sigma <- par[2]
    (k / sigma) * (t / sigma)^(k - 1)
  }
)

# Simulate Weibull data (shape=2, scale=3)
set.seed(456)
true_shape <- 2
true_scale <- 3
n <- 100

# Use inverse CDF sampling
u <- runif(n)
weibull_times <- true_scale * (-log(u))^(1/true_shape)

df_weibull <- data.frame(t = weibull_times, delta = rep(1, n))

# Fit
solver <- fit(weibull)
result <- solver(df_weibull, par = c(1.5, 2.5), method = "BFGS")

c(fitted_shape = coef(result)[1], true_shape = true_shape)
c(fitted_scale = coef(result)[2], true_scale = true_scale)

## ----fig.alt="Bathtub hazard curve showing three phases: high infant mortality at t=0 that decreases, then a constant useful life period, followed by increasing wear-out hazard."----
# h(t) = a * exp(-b*t) + c + d * t^k
# Three components: infant mortality + baseline + wear-out
bathtub <- dfr_dist(
  rate = function(t, par, ...) {
    a <- par[1]  # infant mortality magnitude
    b <- par[2]  # infant mortality decay rate
    c <- par[3]  # baseline hazard (useful life)
    d <- par[4]  # wear-out coefficient
    k <- par[5]  # wear-out exponent
    a * exp(-b * t) + c + d * t^k
  },
  par = c(a = 1, b = 2, c = 0.02, d = 0.001, k = 2)
)

# Plot the hazard function
t_seq <- seq(0.01, 15, length.out = 200)
h <- hazard(bathtub)
plot(t_seq, sapply(t_seq, h), type = "l",
     xlab = "Time", ylab = "Hazard rate",
     main = "Bathtub hazard curve")

## -----------------------------------------------------------------------------
# Proportional hazards with covariate x
# h(t, x) = h0(t) * exp(beta * x)
# where h0(t) = Weibull baseline

ph_model <- dfr_dist(
  rate = function(t, par, x = 0, ...) {
    k <- par[1]
    sigma <- par[2]
    beta <- par[3]
    baseline <- (k / sigma) * (t / sigma)^(k - 1)
    baseline * exp(beta * x)
  },
  par = c(shape = 2, scale = 3, beta = 0.5)
)

h <- hazard(ph_model)

# Hazard for different covariate values
h(2, x = 0)  # Baseline
h(2, x = 1)  # Higher risk group

## -----------------------------------------------------------------------------
# Support is (0, Inf) for all DFR distributions
support <- sup(exp_dist)
print(support)

# Parameters
params(exp_dist)

## ----warning=FALSE, fig.alt="Cox-Snell residuals Q-Q plot showing model fit assessment."----
# Fit a model and check residuals
set.seed(99)
test_times <- rexp(80, rate = 0.3)
test_df <- data.frame(t = test_times, delta = 1)

# Fit exponential
exp_fitted <- dfr_exponential()
solver <- fit(exp_fitted)
fit_result <- solver(test_df, par = c(0.5))
lambda_hat <- coef(fit_result)

# Create fitted distribution with estimated parameters
exp_final <- dfr_exponential(lambda = lambda_hat)

# Q-Q plot of Cox-Snell residuals
qqplot_residuals(exp_final, test_df)

## -----------------------------------------------------------------------------
mart_resid <- residuals(exp_final, test_df, type = "martingale")
summary(mart_resid)

# Large positive values: failed earlier than expected
# Large negative values: survived longer than expected

