## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(flexhaz)

## ----simulate-data------------------------------------------------------------
set.seed(42)
true_k <- 2
true_sigma <- 3
n <- 100

u <- runif(n)
times <- true_sigma * (-log(u))^(1/true_k)
df <- data.frame(t = times, delta = rep(1, n))

cat("Sample size:", n, "\n")
cat("True parameters: k =", true_k, ", sigma =", true_sigma, "\n")

## ----minimal-approach---------------------------------------------------------
weibull_minimal <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (k / sigma) * (t / sigma)^(k - 1)
    }
)

ll <- loglik(weibull_minimal)
s <- score(weibull_minimal)
H <- hess_loglik(weibull_minimal)

test_par <- c(1.8, 2.8)

cat("Log-likelihood:", ll(df, par = test_par), "\n")
cat("Score (numerical):", s(df, par = test_par), "\n")
cat("Hessian (numerical):\n")
print(round(H(df, par = test_par), 4))

## ----analytical-approach------------------------------------------------------
weibull_full <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (k / sigma) * (t / sigma)^(k - 1)
    },
    cum_haz_rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (t / sigma)^k
    },
    score_fn = function(df, par, ...) {
        k <- par[1]
        sigma <- par[2]
        t <- df$t
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))

        n_events <- sum(delta == 1)
        t_ratio <- t / sigma
        log_t_ratio <- log(t_ratio)

        dk <- n_events / k + sum(delta * log_t_ratio) -
            sum(t_ratio^k * log_t_ratio)
        dsigma <- -n_events * k / sigma +
            (k / sigma) * sum(t_ratio^k)

        c(dk, dsigma)
    },
    hess_fn = function(df, par, ...) {
        k <- par[1]
        sigma <- par[2]
        t <- df$t
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))

        n_events <- sum(delta == 1)
        t_ratio <- t / sigma
        log_t_ratio <- log(t_ratio)
        t_ratio_k <- t_ratio^k

        d2k <- -n_events / k^2 - sum(t_ratio_k * log_t_ratio^2)
        d2k_sigma <- -n_events / sigma +
            (1 / sigma) * sum(t_ratio_k) +
            (k / sigma) * sum(t_ratio_k * log_t_ratio)
        d2sigma <- n_events * k / sigma^2 -
            k * (k + 1) / sigma^2 * sum(t_ratio_k)

        matrix(c(d2k, d2k_sigma, d2k_sigma, d2sigma), nrow = 2)
    }
)

s_full <- score(weibull_full)
H_full <- hess_loglik(weibull_full)

cat("Score (analytical):", s_full(df, par = test_par), "\n")
cat("Hessian (analytical):\n")
print(round(H_full(df, par = test_par), 4))

## ----score-only---------------------------------------------------------------
weibull_score_only <- dfr_dist(
    rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (k / sigma) * (t / sigma)^(k - 1)
    },
    cum_haz_rate = function(t, par, ...) {
        k <- par[1]
        sigma <- par[2]
        (t / sigma)^k
    },
    score_fn = function(df, par, ...) {
        k <- par[1]
        sigma <- par[2]
        t <- df$t
        delta <- if ("delta" %in% names(df)) df$delta else rep(1, nrow(df))

        n_events <- sum(delta == 1)
        t_ratio <- t / sigma
        log_t_ratio <- log(t_ratio)

        dk <- n_events / k + sum(delta * log_t_ratio) -
            sum(t_ratio^k * log_t_ratio)
        dsigma <- -n_events * k / sigma +
            (k / sigma) * sum(t_ratio^k)

        c(dk, dsigma)
    }
    # No hess_fn -> numDeriv::hessian fallback
)

H_mixed <- hess_loglik(weibull_score_only)
cat("Hessian (numerical fallback):\n")
print(round(H_mixed(df, par = test_par), 4))

## ----constructors, warning=FALSE----------------------------------------------
# These include score_fn and hess_fn (where available)
weib <- dfr_weibull(shape = 2, scale = 3)
exp_d <- dfr_exponential(lambda = 0.5)

# Fit the model
solver <- fit(weib)
result <- solver(df, par = c(1.5, 2.5))

cat("Fitted parameters:\n")
print(coef(result))
cat("\nTrue parameters: k =", true_k, ", sigma =", true_sigma, "\n")
cat("\n95% Confidence intervals:\n")
print(confint(result))

## ----standard-errors----------------------------------------------------------
# Hessian at MLE
mle_par <- coef(result)
hess_mle <- H_full(df, par = mle_par)

# Observed Fisher information = -Hessian
obs_info <- -hess_mle

# Standard errors = sqrt(diag(inverse Fisher info))
se <- sqrt(diag(solve(obs_info)))

cat("Standard errors:\n")
cat("  SE(k) =", round(se[1], 4), "\n")
cat("  SE(sigma) =", round(se[2], 4), "\n")

