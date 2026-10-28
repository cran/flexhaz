get_params <- function(par, default = NULL) {
    if (is.null(par) || all(is.na(par))) {
        return(default)
    }
    if (!is.null(default) && length(par) != length(default)) {
        stop(sprintf(
            "Parameter length mismatch: got %d, expected %d",
            length(par), length(default)
        ))
    }
    if (!is.null(default)) {
        par[is.na(par)] <- default[is.na(par)]
    }
    par
}

require_params <- function(par, default) {
    par <- get_params(par, default)
    if (is.null(par)) {
        stop("Parameters required: provide via 'par' argument or in distribution object")
    }
    par
}

get_delta <- function(df, delta_col = "delta") {
    delta <- if (delta_col %in% names(df)) df[[delta_col]] else rep(1, nrow(df))
    if (!(is.numeric(delta) || is.logical(delta)) ||
        !is.null(dim(delta)) || length(delta) != nrow(df) ||
        anyNA(delta) || any(!delta %in% c(-1, 0, 1))) {
        stop("event indicators must be non-missing numeric or logical values ",
             "coded 1 (exact), 0 (right-censored), or -1 (left-censored)",
             call. = FALSE)
    }
    delta
}

# Validate once at each inference entry point, including analytical derivatives
# that bypass loglik(). A missing status column still means all events observed.
get_observations <- function(df, ob_col = "t", delta_col = "delta") {
    if (!is.data.frame(df)) {
        stop("observations must be a data frame", call. = FALSE)
    }
    t <- df[[ob_col]]
    if (!is.numeric(t) || !is.null(dim(t)) || length(t) != nrow(df) ||
        any(!is.finite(t)) || any(t < 0)) {
        stop("observation times in '", ob_col,
             "' must be finite, non-negative numeric values", call. = FALSE)
    }
    list(t = t, delta = get_delta(df, delta_col))
}
