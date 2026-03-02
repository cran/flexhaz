## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----eval=FALSE---------------------------------------------------------------
# # From r-universe
# install.packages("flexhaz", repos = "https://queelius.r-universe.dev")

## ----setup--------------------------------------------------------------------
library(flexhaz)

## -----------------------------------------------------------------------------
# Exponential with failure rate lambda = 0.1 (MTTF = 10 time units)
exp_dist <- dfr_exponential(lambda = 0.1)
print(exp_dist)

## -----------------------------------------------------------------------------
# Get closures
S <- surv(exp_dist)
h <- hazard(exp_dist)
f <- density(exp_dist)

# Evaluate at specific times
S(10)   # P(survive past t=10) = exp(-0.1 * 10) ~ 0.37
h(10)   # Hazard at t=10 = 0.1 (constant for exponential)
f(10)   # PDF at t=10

## -----------------------------------------------------------------------------
samp <- sampler(exp_dist)
set.seed(42)
failure_times <- samp(20)

# Sample mean should be close to MTTF = 1/lambda = 10
mean(failure_times)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 50
df <- data.frame(
  t = rexp(n, rate = 0.08),   # True lambda = 0.08
  delta = rep(1, n)            # All exact observations
)
head(df)

## ----warning=FALSE------------------------------------------------------------
# Template with no parameters (to be estimated)
exp_template <- dfr_exponential()

solver <- fit(exp_template)
result <- solver(df, par = c(0.1))  # Initial guess
coef(result)

# Compare to analytical MLE: lambda_hat = n / sum(t)
n / sum(df$t)

## ----fig.alt="Cox-Snell residuals Q-Q plot showing points close to diagonal line, indicating good fit."----
fitted_dist <- dfr_exponential(lambda = coef(result))
qqplot_residuals(fitted_dist, df)

## ----warning=FALSE------------------------------------------------------------
df_censored <- data.frame(
  t = c(5, 8, 12, 15, 20, 25, 30, 30, 30, 30),
  delta = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)  # Last 5 right-censored at t=30
)

solver <- fit(dfr_exponential())
result <- solver(df_censored, par = c(0.1))
coef(result)

## ----warning=FALSE------------------------------------------------------------
df_left <- data.frame(
  t = c(5, 10, 15, 20),
  delta = c(-1, -1, 1, 0)   # left-censored, left-censored, exact, right-censored
)

solver <- fit(dfr_exponential())
result <- solver(df_left, par = c(0.1))
coef(result)

## ----warning=FALSE------------------------------------------------------------
clinical_data <- data.frame(
  time = c(5, 8, 12, 15, 20),
  status = c(1, 1, 1, 0, 0)
)

dist <- dfr_exponential()
dist$ob_col <- "time"
dist$delta_col <- "status"

solver <- fit(dist)
result <- solver(clinical_data, par = c(0.1))
coef(result)

## ----eval=FALSE---------------------------------------------------------------
# dist <- dfr_dist(
#   rate = function(t, par, ...) rep(par[1], length(t)),
#   cum_haz_rate = function(t, par, ...) par[1] * t,
#   ob_col = "time",
#   delta_col = "status"
# )

## ----fig.alt="Two hazard curves: exponential (flat line) and Weibull shape=2 (increasing curve)."----
# Weibull with increasing failure rate (wear-out)
weib_dist <- dfr_weibull(shape = 2, scale = 20)

# Compare hazard functions
plot(weib_dist, what = "hazard", xlim = c(0, 30), col = "blue",
     main = "Hazard Comparison")
plot(dfr_exponential(lambda = 0.05), what = "hazard", add = TRUE, col = "red")
legend("topleft", c("Weibull (shape=2)", "Exponential"),
       col = c("blue", "red"), lwd = 2)

## ----define-dist--------------------------------------------------------------
bathtub <- dfr_dist(
  rate = function(t, par, ...) {
    a <- par[[1]]  # infant mortality magnitude
    b <- par[[2]]  # infant mortality decay
    c <- par[[3]]  # baseline hazard
    d <- par[[4]]  # wear-out coefficient
    k <- par[[5]]  # wear-out exponent
    a * exp(-b * t) + c + d * t^k
  },
  par = c(a = 0.5, b = 0.5, c = 0.01, d = 5e-5, k = 2.5)
)

## ----plot-hazard, fig.alt="Bathtub hazard curve showing high infant mortality at left, flat useful life in the middle, and rising wear-out on the right."----
h <- hazard(bathtub)
t_seq <- seq(0.01, 30, length.out = 300)
plot(t_seq, sapply(t_seq, h), type = "l", lwd = 2, col = "darkred",
     xlab = "Time", ylab = "Hazard rate h(t)",
     main = "Bathtub Hazard Function")

## ----plot-surv, fig.alt="Survival curve derived from the bathtub hazard, showing rapid early decline that stabilizes then drops again."----
plot(bathtub, what = "survival", xlim = c(0, 30),
     col = "darkblue", lwd = 2,
     main = "Survival Function S(t)")

## ----simulate-fit, warning=FALSE----------------------------------------------
# Generate failure times, right-censored at t = 25
set.seed(42)
samp <- sampler(bathtub)
raw_times <- samp(80)
censor_time <- 25

df <- data.frame(
  t = pmin(raw_times, censor_time),
  delta = as.integer(raw_times <= censor_time)
)

cat("Observed failures:", sum(df$delta), "out of", nrow(df), "\n")
cat("Right-censored:", sum(1 - df$delta), "\n")

# Fit the model
solver <- fit(bathtub)
result <- solver(df,
  par = c(0.3, 0.3, 0.02, 1e-4, 2),
  method = "Nelder-Mead"
)

## ----inspect------------------------------------------------------------------
cat("Estimated parameters:\n")
print(coef(result))

cat("\nTrue parameters:\n")
print(c(a = 0.5, b = 0.5, c = 0.01, d = 5e-5, k = 2.5))

## ----qq-plot, fig.alt="Cox-Snell residuals Q-Q plot for the bathtub model."----
fitted_bathtub <- dfr_dist(
  rate = function(t, par, ...) {
    par[[1]] * exp(-par[[2]] * t) + par[[3]] + par[[4]] * t^par[[5]]
  },
  par = coef(result)
)
qqplot_residuals(fitted_bathtub, df)

