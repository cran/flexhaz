## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(flexhaz)

## -----------------------------------------------------------------------------
# Simulated capacitor failure data (hours)
set.seed(42)
n_tested <- 50
max_test_time <- 1000

# Generate from Weibull (unknown to analyst)
true_shape <- 2.3
true_scale <- 800
raw_times <- rweibull(n_tested, shape = true_shape, scale = true_scale)

# Apply censoring at test end
capacitor_data <- data.frame(
  hours = pmin(raw_times, max_test_time),
  failed = as.integer(raw_times <= max_test_time)
)

# Summary
cat("Failures:", sum(capacitor_data$failed), "\n")
cat("Censored:", sum(1 - capacitor_data$failed), "\n")
head(capacitor_data)

## ----warning=FALSE------------------------------------------------------------
# Prepare data in flexhaz format
df <- data.frame(t = capacitor_data$hours, delta = capacitor_data$failed)

# Fit exponential
exp_solver <- fit(dfr_exponential())
exp_result <- exp_solver(df, par = c(0.001))
exp_lambda <- coef(exp_result)
cat("Exponential lambda:", round(exp_lambda, 6), "\n")

# Fit Weibull
weib_solver <- fit(dfr_weibull())
weib_result <- weib_solver(df, par = c(2, 600))
weib_params <- coef(weib_result)
cat("Weibull shape:", round(weib_params[1], 3), "\n")
cat("Weibull scale:", round(weib_params[2], 1), "\n")

## -----------------------------------------------------------------------------
# Compute log-likelihood at fitted parameters
ll_exp <- loglik(dfr_exponential())
ll_weib <- loglik(dfr_weibull())

aic_exp <- -2 * ll_exp(df, exp_lambda) + 2 * 1   # 1 parameter
aic_weib <- -2 * ll_weib(df, weib_params) + 2 * 2  # 2 parameters

cat("AIC (Exponential):", round(aic_exp, 2), "\n")
cat("AIC (Weibull):", round(aic_weib, 2), "\n")
cat("Winner:", ifelse(aic_weib < aic_exp, "Weibull", "Exponential"), "\n")

## ----fig.alt="Hazard comparison showing exponential (flat line) vs Weibull (increasing curve) for capacitors."----
# Create fitted distributions
exp_fit <- dfr_exponential(lambda = exp_lambda)
weib_fit <- dfr_weibull(shape = weib_params[1], scale = weib_params[2])

# Compare hazard functions
plot(weib_fit, what = "hazard", xlim = c(0, 1200),
     col = "blue", lwd = 2, main = "Capacitor Failure Rate")
plot(exp_fit, what = "hazard", add = TRUE, col = "red", lwd = 2, lty = 2)
abline(h = exp_lambda, col = "red", lwd = 1, lty = 3)
legend("topleft", c("Weibull", paste0("Exponential (\u03bb=", round(exp_lambda, 4), ")")),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))

## -----------------------------------------------------------------------------
# Use fitted Weibull distribution
Q <- inv_cdf(weib_fit)

# Calculate B-lives
B10 <- Q(0.10)   # 10% failure quantile
B50 <- Q(0.50)   # Median
B90 <- Q(0.90)   # 90% failure quantile

cat("B10 life:", round(B10, 1), "hours\n")
cat("B50 life:", round(B50, 1), "hours\n")
cat("B90 life:", round(B90, 1), "hours\n")

## ----fig.alt="Survival curve with horizontal lines showing B10, B50, and B90 life points."----
plot(weib_fit, what = "survival", xlim = c(0, 1200),
     main = "Capacitor Reliability Curve",
     col = "darkblue", lwd = 2)

# Add B-life reference lines
abline(h = c(0.90, 0.50, 0.10), lty = 2, col = "gray50")
abline(v = c(B10, B50, B90), lty = 2, col = "gray50")

# Labels
text(B10, 0.95, paste0("B10=", round(B10)), pos = 4, cex = 0.8)
text(B50, 0.55, paste0("B50=", round(B50)), pos = 4, cex = 0.8)
text(B90, 0.15, paste0("B90=", round(B90)), pos = 4, cex = 0.8)

## -----------------------------------------------------------------------------
# Warranty period (convert years to hours: 2 years ≈ 17520 hours)
warranty_hours <- 2 * 365 * 24

# But our capacitor test was at accelerated conditions
# Assume acceleration factor of 20x
field_warranty_equivalent <- warranty_hours / 20

# Fraction failing during warranty
S <- surv(weib_fit)
survival_at_warranty <- S(field_warranty_equivalent)
failure_rate <- 1 - survival_at_warranty

cat("Warranty period (equivalent hours):", round(field_warranty_equivalent, 1), "\n")
cat("Expected survival at warranty end:", round(survival_at_warranty * 100, 1), "%\n")
cat("Expected failure rate:", round(failure_rate * 100, 2), "%\n")
cat("Claims per 1000 units:", round(failure_rate * 1000, 1), "\n")

## ----fig.alt="Gompertz hazard curve showing exponentially increasing failure rate with proposed maintenance interval."----
# System with aging characteristics
system <- dfr_gompertz(a = 0.001, b = 0.002)

# Target: keep hazard below 0.05 (5% per time unit)
h <- hazard(system)

# Find time when hazard reaches threshold
# h(t) = a * exp(b*t) = 0.05
# t = log(0.05/a) / b
threshold <- 0.05
a <- 0.001
b <- 0.002
maintenance_time <- log(threshold / a) / b

cat("Hazard threshold:", threshold, "\n")
cat("Maintenance interval:", round(maintenance_time, 1), "time units\n")

# Visualize
plot(system, what = "hazard", xlim = c(0, 2500),
     main = "Maintenance Planning", col = "darkgreen", lwd = 2)
abline(h = threshold, col = "red", lty = 2, lwd = 2)
abline(v = maintenance_time, col = "blue", lty = 2, lwd = 2)
legend("topleft", c("Hazard", "Threshold", "Maintenance"),
       col = c("darkgreen", "red", "blue"), lty = c(1, 2, 2), lwd = 2)

## ----fig.alt="Hazard decomposition showing electronic and mechanical failure modes combining to bathtub shape."----
# Electronic component: constant failure rate (random defects)
# Mechanical component: Weibull wear-out

dfr_competing <- function(lambda_elec = NULL, k_mech = NULL, sigma_mech = NULL) {
  par <- if (!is.null(lambda_elec) && !is.null(k_mech) && !is.null(sigma_mech)) {
    c(lambda_elec, k_mech, sigma_mech)
  } else NULL

  dfr_dist(
    rate = function(t, par, ...) {
      lambda <- par[[1]]
      k <- par[[2]]
      sigma <- par[[3]]
      # Combined hazard = electronic + mechanical
      lambda + (k / sigma) * (t / sigma)^(k - 1)
    },
    par = par
  )
}

# Product with both failure modes
product <- dfr_competing(lambda_elec = 0.0005, k_mech = 2.5, sigma_mech = 10000)

# Decompose contributions
t_seq <- seq(1, 15000, length.out = 200)
h_total <- sapply(t_seq, function(ti) hazard(product)(ti))
h_elec <- rep(0.0005, length(t_seq))
h_mech <- (2.5 / 10000) * (t_seq / 10000)^(2.5 - 1)

# Plot decomposition
plot(t_seq, h_total, type = "l", col = "black", lwd = 3,
     xlab = "Time (hours)", ylab = "Hazard rate",
     main = "Competing Failure Modes")
lines(t_seq, h_elec, col = "blue", lwd = 2, lty = 2)
lines(t_seq, h_mech, col = "red", lwd = 2, lty = 2)
legend("topleft", c("Total", "Electronic", "Mechanical"),
       col = c("black", "blue", "red"), lwd = c(3, 2, 2), lty = c(1, 2, 2))

## ----fig.alt="Cox-Snell residuals Q-Q plot for Weibull model showing good fit."----
# Check Weibull fit for capacitor data
qqplot_residuals(weib_fit, df)

