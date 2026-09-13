#' Generalised logistic function
#'
#' @description
#' Calculate a 4- or 5-parameter logistic (sigmoidal) curve. The 4-parameter
#' symmetric form is fit by [analyse_kinetics()] with `method = "sigmoidal"`
#' and `shape = "symmetric"` (*default*), and by [stats::nls()] via the
#' self-starting wrapper [SSlogistic()].
#'
#' @param t A numeric vector of the predictor variable (time).
#' @param A A numeric parameter for the starting asymptote of the response
#'   variable.
#' @param B A numeric parameter for the ending asymptote of the response
#'   variable.
#' @param xmid A numeric parameter for the time at the *inflection point* (the
#'   steepest point) of the curve, in units of the predictor variable `t`.
#' @param slope A numeric parameter for the response rate `dx/dt` at the
#'   inflection `xmid`.
#' @param asym A numeric parameter for the asymmetry index of the curve; the
#'   fraction of the amplitude `(y(xmid) - A) / (B - A)` at which the
#'   inflection `xmid` occurs, in `(0, 1)`. `asym = 0.5` is symmetric and
#'   equivalent to the 4-parameter form. If `NULL` (*default*), a symmetric
#'   4-parameter model is used.
#'
#' @details
#' The 5-parameter Richards form is exported for advanced use directly with
#' [stats::nls()] but is not used by [analyse_kinetics()] due to convergence
#' instability. For asymmetric responses, prefer [gompertz()] /
#' [gompertz_left()], which are more stable.
#'
#' ## Model equations
#'
#' Both forms are re-parameterised from the Richards generalised logistic
#' model so `xmid` is the time at inflection and `slope` is the response rate
#' `dx/dt` at the inflection.
#'
#' - 4-parameter (symmetric):
#'   `A + (B - A) / (1 + exp(-4 * slope * (t - xmid) / (B - A)))`
#' - 5-parameter (asymmetric):
#'   `A + (B - A) / (1 + exp(-k * (t - xmid)))^(1 / v)` with
#'   `v = -log(2) / log(asym)` and `k = 2 * slope * v / ((B - A) * asym)`.
#'
#' The inflection is at `t = xmid` with `dx/dt = slope` and
#' `y(xmid) = A + (B - A) * asym` for any `asym` in `(0, 1)`:
#'
#' - `asym = 0.5` (`v = 1`) collapses to the 4-parameter form.
#' - `asym -> 0` gives an early-acceleration curve (inflection near `A`).
#' - `asym -> 1` gives a late-acceleration curve (inflection near `B`).
#' - `asym = 0.368` (`1/e`) approximates a right-inflection [gompertz()] curve.
#' - `asym = 0.632` (`1 - 1/e`) approximates a left-inflection
#'   [gompertz_left()] curve.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSlogistic()], [gompertz()],
#'   [gompertz_left()], [sigmoidal_drift()], [monoexponential()]
#'
#' @examples
#' ## create an asymmetric logistic curve with random noise
#' set.seed(15)
#' t <- 1:60
#' x <- logistic(t, A = 10, B = 100, xmid = 30, slope = 4, asym = 0.3) +
#'     rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## 5-parameter fit with the self-starting wrapper
#' model <- nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data)
#' summary(model)
#'
#' y <- predict(model, data)
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         ggplot2::ggplot(data, ggplot2::aes(t, x)) +
#'             theme_mnirs() +
#'             ggplot2::geom_point() +
#'             ggplot2::geom_line(ggplot2::aes(y = y))
#'     }
#' }
#'
#' @export
logistic <- function(t, A, B, xmid, slope, asym = NULL) {
    if (is.null(asym)) {
        ## 4-parameter symmetric
        y <- A + (B - A) / (1 + exp(-4 * slope * (t - xmid) / (B - A)))
    } else {
        ## 5-parameter Richards re-parameterised: asym is the
        ## inflection-height fraction (y(xmid) - A) / (B - A).
        ## clamp to keep log(asym) finite during nls iteration
        asym <- min(max(asym, 1e-6), 1 - 1e-6)
        v <- -log(2) / log(asym)
        k <- 2 * slope * v / ((B - A) * asym)
        y <- A + (B - A) / (1 + exp(-k * (t - xmid)))^(1 / v)
    }
    return(y)
}


#' Gompertz growth functions
#'
#' @description
#' Calculate 4-parameter Gompertz (asymmetric sigmoidal) curves. Model
#' families fit by [analyse_kinetics()] with `method = "sigmoidal"` and
#' `shape = "gompertz"` or `"gompertz_left"`, and by [stats::nls()] via the
#' self-starting wrappers [SSgompertz()] and [SSgompertz_left()].
#'
#' @inheritParams logistic
#'
#' @details
#' `gompertz()` (right-Gompertz) is asymmetric with the inflection point
#' `xmid` closer to the starting asymptote `A`: early acceleration away from
#' `A`, and a slow approach to the ending asymptote `B`. Appropriate for
#' fast-onset, slow-tail responses.
#'
#' `gompertz_left()` (left-Gompertz) has the inflection point closer to the
#' ending asymptote `B`: slow departure from `A`, and late acceleration toward
#' `B`. Appropriate for slow-onset, fast-tail responses.
#'
#' ## Model equations
#'
#' Both forms are re-parameterised so `xmid` is the time at inflection and
#' `slope` is the response rate `dx/dt` at the inflection, with
#' `k = slope * e / (B - A)`.
#'
#' - `gompertz()`: `A + (B - A) * exp(-exp(-k * (t - xmid)))`. Inflection
#'   height fixed at `A + (B - A) / e`; 36.8% of the amplitude.
#' - `gompertz_left()`: `A + (B - A) * (1 - exp(-exp(k * (t - xmid))))`.
#'   Inflection height fixed at `A + (B - A) * (1 - 1/e)`; 63.2% of the
#'   amplitude.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSgompertz()], [SSgompertz_left()],
#'   [logistic()], [sigmoidal_drift()]
#'
#' @examples
#' ## create a Gompertz curve with random noise
#' set.seed(15)
#' t <- 1:60
#' x <- gompertz(t, A = 10, B = 100, xmid = 30, slope = 4) +
#'     rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## fit with the self-starting wrapper
#' model <- nls(x ~ SSgompertz(t, A, B, xmid, slope), data = data)
#' summary(model)
#'
#' y <- predict(model, data)
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         ggplot2::ggplot(data, ggplot2::aes(t, x)) +
#'             theme_mnirs() +
#'             ggplot2::geom_point() +
#'             ggplot2::geom_line(ggplot2::aes(y = y))
#'     }
#' }
#'
#' @export
gompertz <- function(t, A, B, xmid, slope) {
    k <- slope * exp(1) / (B - A)
    y <- A + (B - A) * exp(-exp(-k * (t - xmid)))
    return(y)
}


#' @rdname gompertz
#' @export
gompertz_left <- function(t, A, B, xmid, slope) {
    k <- slope * exp(1) / (B - A)
    y <- A + (B - A) * (1 - exp(-exp(k * (t - xmid))))
    return(y)
}


#' Sigmoid curve with gradient
#'
#' [sigmoid_core()] evaluates a 4-parameter sigmoid of the given `shape`
#' and its partial derivatives on the canonical parameters, shared by the
#' `selfStart` model functions of [SSlogistic()], [SSgompertz()],
#' [SSgompertz_left()], and [SSsigmoidal_drift()]. Every shape is a
#' function `W(u)` of `u = k * (t - xmid)` with rate `k = c * slope /
#' (B - A)` (`c = 4` symmetric, `e` Gompertz), so with `P = dW/du` the
#' partials share one form. [sigmoid_model()] attaches the gradient over
#' the parameters written as bare symbols in `mCall` (see
#' [free_params()]), so [stats::nls()] skips [stats::numericDeriv()].
#'
#' @param mCall A matched call to the model function.
#' @inheritParams sigmoidal_drift
#'
#' @returns [sigmoid_core()]: a list of the curve `val`, the partial
#'   derivatives by parameter name, and the rate `k`. [sigmoid_model()]: a
#'   numeric vector of predicted values with a `"gradient"` attribute when
#'   any parameter is free.
#'
#' @keywords internal
sigmoid_core <- function(t, A, B, xmid, slope, shape) {
    D <- B - A
    cc <- if (shape == "symmetric") 4 else exp(1)
    k <- cc * slope / D
    u <- k * (t - xmid)
    ## shape response W on [0, 1] and its rate P = dW/du
    if (shape == "symmetric") {
        W <- 1 / (1 + exp(-u))
        P <- W * (1 - W)
    } else if (shape == "gompertz") {
        W <- exp(-exp(-u))
        P <- W * exp(-u)
    } else {
        H <- exp(-exp(u))
        W <- 1 - H
        P <- H * exp(u)
    }
    uP <- u * P
    return(list(
        val = A + D * W,
        A = 1 - W + uP,
        B = W - uP,
        xmid = -D * k * P,
        slope = cc * (t - xmid) * P,
        k = k
    ))
}


#' @rdname sigmoid_core
#' @keywords internal
sigmoid_model <- function(mCall, t, A, B, xmid, slope, shape) {
    g <- sigmoid_core(t, A, B, xmid, slope, shape)
    val <- g$val
    free <- free_params(mCall, c("A", "B", "xmid", "slope"))
    if (length(free) > 0L) {
        attr(val, "gradient") <- do.call(cbind, g[free])
    }
    return(val)
}


## `selfStart` model fns: the exported curves plus the analytic gradient.
## the 5-parameter logistic has no gradient (numericDeriv)
logistic_model <- function(t, A, B, xmid, slope, asym = NULL) {
    if (!is.null(asym)) {
        return(logistic(t, A, B, xmid, slope, asym))
    }
    return(sigmoid_model(match.call(), t, A, B, xmid, slope, "symmetric"))
}

gompertz_model <- function(t, A, B, xmid, slope) {
    return(sigmoid_model(match.call(), t, A, B, xmid, slope, "gompertz"))
}

gompertz_left_model <- function(t, A, B, xmid, slope) {
    return(sigmoid_model(match.call(), t, A, B, xmid, slope, "gompertz_left"))
}


#' Initiate self-starting logistic model
#'
#' [logistic_init()]: Returns initial values for the parameters in a
#' `selfStart` model.
#'
#' @param mCall A matched call to the function `model`.
#' @param data A data frame with predictor `t` and the response variable.
#' @param LHS The left-hand side expression of the model formula.
#' @param ... Additional arguments, including `fixed`, a named list of
#'   user-fixed parameter values from [init_fixed()] used to seed the
#'   remaining free estimates.
#'
#' @returns [logistic_init()]: Initial starting estimates for parameters in
#'   the model called by [SSlogistic()].
#'
#' @keywords internal
logistic_init <- function(mCall, data, LHS, ...) {
    tx <- sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)
    has_asym <- "asym" %in% names(mCall)

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()

    ## asymptotes from first and last ceiling(n/5) values
    ab <- init_asymptotes(x, n)
    A_init <- fixed$A %||% ab$A
    B_init <- fixed$B %||% ab$B

    ## linearisation for 4-param: log((B - y) / (y - A)) ~ t
    lo <- min(A_init, B_init)
    hi <- max(A_init, B_init)
    eps <- (hi - lo) * 1e-3
    x_clip <- pmin(pmax(x, lo + eps), hi - eps)
    xf <- log((B_init - x_clip) / (x_clip - A_init))

    xmid_init <- NA_real_
    slope_init <- NA_real_
    finite_idx <- is.finite(xf)
    if (sum(finite_idx) >= 3L) {
        b <- slope(
            xf[finite_idx],
            t[finite_idx],
            intercept = TRUE,
            bypass_checks = TRUE,
            min_obs = 2L
        )
        a <- attr(b, "intercept")
        if (is.finite(b) && b != 0) {
            xmid_init <- -a / b
            slope_init <- -b * (B_init - A_init) / 4
        }
    }

    ## fallbacks for degenerate data
    t_range <- diff(range(t))
    if (!is.finite(xmid_init) || xmid_init < min(t) || xmid_init > max(t)) {
        xmid_init <- t[which.min(abs(x - (A_init + B_init) / 2))]
    }
    if (!is.finite(slope_init) || slope_init == 0) {
        slope_init <- if (t_range > 0) {
            (B_init - A_init) / t_range
        } else {
            sign(B_init - A_init)
        }
    }

    if (!has_asym) {
        return(c(
            A = A_init,
            B = B_init,
            xmid = fixed$xmid %||% xmid_init,
            slope = fixed$slope %||% slope_init
        ))
    }

    ## 5-param: empirical inflection from smoothed derivative
    infl <- init_inflection(x, t, A_init, B_init)
    asym_emp <- (x[infl$idx] - A_init) / (B_init - A_init)
    asym_init <- min(max(asym_emp, 0.1), 0.9)

    return(c(
        A = A_init,
        B = B_init,
        xmid = fixed$xmid %||% infl$xmid,
        slope = fixed$slope %||% infl$slope,
        asym = fixed$asym %||% asym_init
    ))
}


#' Initiate self-starting Gompertz model
#'
#' [gompertz_init()]: Returns initial values for the parameters in a
#' `selfStart` model. Used by both [SSgompertz()] and [SSgompertz_left()];
#' the symmetric logistic linearisation does not apply to Gompertz forms, so
#' initialisation is derivative-based via [init_inflection()].
#'
#' @inheritParams logistic_init
#'
#' @returns [gompertz_init()]: Initial starting estimates for parameters
#'   in the model called by [SSgompertz()] or [SSgompertz_left()].
#'
#' @keywords internal
gompertz_init <- function(mCall, data, LHS, ...) {
    tx <- sortedXyData(mCall[["t"]], LHS, data)
    x <- tx[["y"]]
    t <- tx[["x"]]
    n <- length(x)

    ## user-fixed parameter values seed the remaining free estimates
    fixed <- list(...)$fixed %||% list()

    ab <- init_asymptotes(x, n)
    A_init <- fixed$A %||% ab$A
    B_init <- fixed$B %||% ab$B
    infl <- init_inflection(x, t, A_init, B_init)

    return(c(
        A = A_init,
        B = B_init,
        xmid = fixed$xmid %||% infl$xmid,
        slope = fixed$slope %||% infl$slope
    ))
}


#' Estimate baseline and asymptote from the first/last quintile of `x`
#'
#' Shared helper used by self-start initialisers for logistic / Gompertz
#' model families.
#'
#' @param x A numeric vector of the response variable (sorted by `t`).
#' @param n An integer length of `x`.
#'
#' @returns A list with elements `A` (starting asymptote estimate) and `B`
#'   (ending asymptote estimate).
#'
#' @keywords internal
init_asymptotes <- function(x, n = length(x)) {
    n_asymp <- max(1L, ceiling(n / 5))
    A_init <- mean(x[seq_len(n_asymp)])
    B_init <- mean(x[seq(n - n_asymp + 1L, n)])
    return(list(A = A_init, B = B_init))
}


#' Estimate inflection point from a smoothed first derivative
#'
#' Shared helper that locates the empirical inflection (peak of
#' `|dx/dt|` after smoothing) and returns the corresponding `xmid` and
#' `slope` initial values. Falls back to the half-response point and a
#' mean-rate slope when the derivative is degenerate.
#'
#' @param x A numeric vector of the response variable (sorted by `t`).
#' @param t A numeric vector of the predictor variable.
#' @param A_init Estimated starting asymptote.
#' @param B_init Estimated ending asymptote.
#'
#' @returns A list with elements `idx` (integer index into `x`), `xmid`
#'   (numeric `t` value at the inflection), and `slope` (numeric `dx/dt`
#'   at the inflection).
#'
#' @keywords internal
init_inflection <- function(x, t, A_init, B_init) {
    dx_dt <- diff(x) / diff(t)
    win <- max(3L, 2L * (length(dx_dt) %/% 20L) + 1L)
    dx_smooth <- as.numeric(stats::filter(dx_dt, rep(1 / win, win), sides = 2L))
    dx_smooth[!is.finite(dx_smooth)] <- 0

    i_infl <- which.max(abs(dx_smooth))
    xmid_init <- t[i_infl]
    slope_init <- dx_smooth[i_infl]

    ## fallback: half-response point with mean-rate slope
    t_range <- diff(range(t))
    if (!is.finite(xmid_init) || xmid_init < min(t) || xmid_init > max(t)) {
        i_infl <- which.min(abs(x - (A_init + B_init) / 2))
        xmid_init <- t[i_infl]
    }
    if (!is.finite(slope_init) || slope_init == 0) {
        slope_init <- if (t_range > 0) {
            (B_init - A_init) / t_range
        } else {
            sign(B_init - A_init)
        }
    }

    return(list(idx = i_infl, xmid = xmid_init, slope = slope_init))
}


#' Self-starting logistic model
#'
#' @description
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [logistic()], for use with [stats::nls()]. Supports both the 4-parameter
#' symmetric (A, B, xmid, slope) and 5-parameter asymmetric (A, B, xmid,
#' slope, asym) forms; arity is inferred from the formula passed to
#' [stats::nls()].
#'
#' @usage
#' SSlogistic(t, A, B, xmid, slope, asym)
#'
#' @inheritParams logistic
#'
#' @details
#' ## Model formulas
#'
#' - 4-parameter: `x ~ SSlogistic(t, A, B, xmid, slope)`
#' - 5-parameter: `x ~ SSlogistic(t, A, B, xmid, slope, asym)`
#'
#' The 4-parameter form is used by [analyse_kinetics()] with
#' `method = "sigmoidal"` and `shape = "symmetric"`. The 5-parameter
#' asymmetric form is retained for advanced/experimental use only;
#' [analyse_kinetics()] instead dispatches to [SSgompertz()] /
#' [SSgompertz_left()] for asymmetric shapes, which are more stable.
#' [stats::nls()] reads the free parameters from the formula right-hand side,
#' so omitting `asym` incurs no degrees-of-freedom penalty.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its name
#' in the formula, e.g. `x ~ SSlogistic(t, A = 0, B, xmid, slope)` fixes the
#' starting asymptote at `A = 0`. Fixed parameters are excluded from
#' estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [logistic()], [analyse_kinetics()], [stats::nls()],
#'   [stats::selfStart()], [stats::SSfpl()], [SSgompertz()]
#'
#' @examples
#' ## create an asymmetric logistic curve with random noise
#' set.seed(15)
#' t <- 1:60
#' x <- logistic(t, A = 10, B = 100, xmid = 30, slope = 4, asym = 0.3) +
#'     rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## 4-parameter fit
#' model4 <- nls(x ~ SSlogistic(t, A, B, xmid, slope), data = data)
#' summary(model4)
#'
#' ## 5-parameter fit on the same data
#' model5 <- nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = data)
#' summary(model5)
#'
#' ## fix the starting asymptote `A` at a known value
#' model_fixed <- nls(x ~ SSlogistic(t, A = 10, B, xmid, slope), data = data)
#' summary(model_fixed)
#'
#' y4 <- predict(model4, data)
#' y5 <- predict(model5, data)
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         ggplot2::ggplot(data, ggplot2::aes(t, x)) +
#'             theme_mnirs() +
#'             ggplot2::geom_point() +
#'             ggplot2::geom_line(ggplot2::aes(y = y5, colour = "5-param")) +
#'             ggplot2::geom_line(ggplot2::aes(y = y4, colour = "4-param"))
#'     }
#' }
#'
#' @export
SSlogistic <- selfStart(
    model = logistic_model,
    initial = init_fixed(logistic_init, c("A", "B", "xmid", "slope", "asym")),
    parameters = c("A", "B", "xmid", "slope", "asym")
)


#' Self-starting Gompertz models
#'
#' @description
#' Creates initial coefficient estimates for `selfStart` wrappers around
#' [gompertz()] and [gompertz_left()], for use with [stats::nls()]. Both
#' wrappers use the same 4-parameter (A, B, xmid, slope) interface.
#'
#' @usage
#' SSgompertz(t, A, B, xmid, slope)
#'
#' SSgompertz_left(t, A, B, xmid, slope)
#'
#' @inheritParams logistic
#'
#' @details
#' ## Model formulas
#'
#' - Right-Gompertz: `x ~ SSgompertz(t, A, B, xmid, slope)`
#' - Left-Gompertz: `x ~ SSgompertz_left(t, A, B, xmid, slope)`
#'
#' Used by [analyse_kinetics()] with `method = "sigmoidal"` and
#' `shape = "gompertz"` or `"gompertz_left"`. Starting estimates locate the
#' inflection from a smoothed first derivative. `SSgompertz()` masks
#' [stats::SSgompertz()].
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its name
#' in the formula, e.g. `x ~ SSgompertz(t, A = 0, B, xmid, slope)` fixes the
#' starting asymptote at `A = 0`. Fixed parameters are excluded from
#' estimation and are not returned by [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [gompertz()], [gompertz_left()], [analyse_kinetics()],
#'   [SSlogistic()], [stats::nls()], [stats::selfStart()],
#'   [stats::SSgompertz()]
#'
#' @examples
#' ## create a Gompertz curve with random noise
#' set.seed(15)
#' t <- 1:60
#' x <- gompertz(t, A = 10, B = 100, xmid = 30, slope = 4) +
#'     rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' model <- nls(x ~ SSgompertz(t, A, B, xmid, slope), data = data)
#' summary(model)
#'
#' ## fix the starting asymptote `A` at a known value
#' model_fixed <- nls(x ~ SSgompertz(t, A = 10, B, xmid, slope), data = data)
#' summary(model_fixed)
#'
#' ## left-Gompertz
#' set.seed(16)
#' x2 <- gompertz_left(t, A = 10, B = 100, xmid = 30, slope = 4) +
#'     rnorm(length(t), 0, 2)
#' data2 <- data.frame(t, x = x2)
#'
#' model_left <- nls(x ~ SSgompertz_left(t, A, B, xmid, slope), data = data2)
#' summary(model_left)
#'
#' @export
SSgompertz <- selfStart(
    model = gompertz_model,
    initial = init_fixed(gompertz_init, c("A", "B", "xmid", "slope")),
    parameters = c("A", "B", "xmid", "slope")
)


#' @rdname SSgompertz
#' @export
SSgompertz_left <- selfStart(
    model = gompertz_left_model,
    initial = init_fixed(gompertz_init, c("A", "B", "xmid", "slope")),
    parameters = c("A", "B", "xmid", "slope")
)


#' Analyse logistic kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "sigmoidal")`. Fits a 4-parameter sigmoidal
#' curve to each `nirs_channel` within a single *"mnirs"* data frame via
#' [fit_sigmoidal()] with one of three shapes: `"symmetric"`, `"gompertz"`,
#' or `"gompertz_left"`. See [analyse_kinetics()] for user-facing
#' documentation.
#'
#' @param shape Character; the 4-parameter sigmoidal shape to fit. One of
#'   `"symmetric"` (*default*; calls [SSlogistic()]), `"gompertz"`
#'   (early-inflection; calls [SSgompertz()]), or `"gompertz_left"`
#'   (late-inflection; calls [SSgompertz_left()]).
#' @param fix An *optional* named list of model parameters (`A`, `B`,
#'   `xmid`, `slope`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Fixed parameters are excluded from estimation
#'   and reported at their fixed values. Applied to every channel, or
#'   per-channel as a list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`.
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#' @inheritParams analyse_monoexponential
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `A`, `B`, `xmid`, `slope`, `xmid_fitted`.
#'   Per-channel metadata are attached as attributes:
#'   - `"model"`: an [nls][stats::nls] model object, or `NULL` for channels
#'     where fitting failed.
#'   - `"fitted_data"`: a named list of per-channel data frames with
#'     columns `window_idx` and `fitted`.
#'   - `"diagnostics"`: a `data.frame` with one row per `nirs_channel`
#'     containing model fit diagnostics.
#'   - `"channel_args"`: a `data.frame` with one row per `nirs_channel`
#'     recording the resolved arguments used.
#'
#' @seealso [analyse_kinetics()], [logistic()], [SSlogistic()],
#'   [gompertz()], [gompertz_left()], [SSgompertz()], [SSgompertz_left()]
#'
#' @keywords internal
analyse_logistic <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    shape = c("symmetric", "gompertz", "gompertz_left"),
    fix = NULL,
    control = NULL,
    start_time = NULL,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    env = rlang::caller_env()
) {
    ## validation ==================================================
    args <- list(...)
    ## interval label; falls back to the `data` argument name when unsupplied
    interval_name <- args$interval_name %||% deparse(substitute(data))

    ## shared prologue: validate data, resolve channels/time, broadcast and
    ## validate per-channel args
    setup <- setup_kinetics_worker(
        data,
        enquo(nirs_channels),
        enquo(time_channel),
        # fmt: skip
        arg_list = mget(c(
            "shape", "fix", "control", "start_time", "direction", "end_window"
        )),
        choices = list(
            shape = c("symmetric", "gompertz", "gompertz_left"),
            direction = c("auto", "positive", "negative")
        ),
        fix_params = c("A", "B", "xmid", "slope"),
        verbose = verbose,
        env = env
    )

    return(analyse_kinetics_channels(
        data,
        setup$nirs_channels,
        setup$time_channel,
        setup$per_channel,
        fit_sigmoidal,
        verbose,
        interval_name,
        extra_args = args,
        method = "sigmoidal",
        env = env
    ))
}


#' Fit a sigmoidal model to one channel
#'
#' Channel fitter of [analyse_logistic()] (see
#' [analyse_kinetics_channels()]), also the fallback of
#' [fit_sigmoidal_drift()]. Self-starting [SSlogistic()], [SSgompertz()],
#' or [SSgompertz_left()] per the channel `shape` via [stats::nls()], with
#' the requested `direction` enforced on `B - A` and the sign of `slope`
#' ([enforce_direction()]).
#'
#' @inheritParams fit_monoexponential
#'
#' @returns The `coefs`/`model`/`fitted_data`/`diag` list of
#'   [build_fit_results()], or [build_na_results()] when the fit fails.
#'
#' @keywords internal
fit_sigmoidal <- function(x, t, valid, .a, ctx) {
    x_fit <- x[valid$idx]
    t_fit <- t[valid$idx]
    ## NA scaffold (method columns only) for convergence failure
    na_cols <- kinetics_coef_cols$sigmoidal
    ## the self-start fn of the channel shape
    fn <- as.name(paste0(
        "SS",
        if (.a$shape == "symmetric") "logistic" else .a$shape
    ))
    ## columns carry the channel names so the model predicts on them
    params <- c("A", "B", "xmid", "slope")
    nm <- fit_names(ctx$nirs, ctx$time_channel, params)
    fit_data <- list2DF(setNames(list(x_fit, t_fit), nm))

    ## build nls formula with any fixed params as constants
    model <- tryCatch(
        embed_fit_call(nls(
            build_ss_formula(fn, params, .a$fix, nm[[1L]], nm[[2L]]),
            fit_data,
            control = fit_control(.a$control)
        )),
        error = \(e) {
            warn_fit_failed(fn, e, ctx$nirs, ctx$interval_name, env = ctx$env)
        }
    )
    if (is.null(model)) {
        return(build_na_results(na_cols))
    }

    coefs <- full_coefs(model, params, .a$fix)

    ## enforce direction: bounded refit on D = B - A and slope sign.
    ## data-scaled slope floor: slope pinned here is a degenerate
    ## flat fit, not a genuine response
    want <- if (.a$direction == "positive") 1 else -1
    slope_eps <- diff(range(x_fit)) / diff(range(t_fit)) * 1e-6
    slope_free <- !"slope" %in% names(.a$fix)
    enforced <- enforce_direction(
        model,
        coefs,
        fit_data,
        direction = .a$direction,
        amp_fn = fn,
        lower = if (slope_free) {
            c(slope = if (want > 0) slope_eps else -Inf)
        },
        upper = if (slope_free) {
            c(slope = if (want > 0) Inf else -slope_eps)
        },
        fix = .a$fix,
        control = .a$control,
        .nirs = ctx$nirs,
        interval_name = ctx$interval_name,
        env = ctx$env
    )
    if (is.null(enforced)) {
        return(build_na_results(na_cols))
    }
    model <- enforced$model
    coefs <- enforced$coefs

    ## predict response at the inflection point xmid, which is already
    ## elapsed from start_time, matching the fit time base
    xmid_fitted <- as.numeric(
        stats::predict(model, setNames(data.frame(coefs[["xmid"]]), nm[[2L]]))
    )

    return(build_fit_results(
        list2DF(list(
            A = coefs[["A"]],
            B = coefs[["B"]],
            xmid = coefs[["xmid"]],
            slope = coefs[["slope"]],
            xmid_fitted = xmid_fitted
        )),
        model,
        x_fit,
        t_fit,
        valid,
        env = ctx$env
    ))
}
