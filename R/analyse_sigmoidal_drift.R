#' Sigmoidal-drift function
#'
#' @description
#' Calculate a two-phase curve: a *fast* sigmoidal primary response of the
#' given `shape` plus a *slow* linear secondary drift beginning near the
#' ending asymptote. Model family fit by [analyse_kinetics()] with
#' `method = "sigmoidal_drift"`, and by [stats::nls()] via the self-starting
#' wrapper [SSsigmoidal_drift()].
#'
#' @param slope_B A numeric parameter for the linear drift rate `dx/dt` of the
#'   secondary phase at the ending asymptote `B`, in response units per unit
#'   of the predictor variable `t`.
#' @param drift_fraction A numeric fraction of the primary amplitude `B - A`
#'   in `(0.5, 1)` at which the linear drift begins, where the sigmoid
#'   reaches `A + drift_fraction * (B - A)`.
#' @param shape Character; the 4-parameter sigmoidal shape. One of
#'   `"symmetric"` (*default*; [logistic()]), `"gompertz"` ([gompertz()]), or
#'   `"gompertz_left"` ([gompertz_left()]).
#' @inheritParams logistic
#'
#' @details
#' ## Model equation
#'
#' `S(t) + slope_B * pmax(t - onset, 0)`
#'
#' `S(t)` is the 4-parameter sigmoid of the given `shape` with asymptotes `A`
#' and `B`, inflection `xmid`, and inflection rate `slope` (see [logistic()]
#' and [gompertz()]). The drift is a hinge line anchored at zero at the onset,
#' so it is exactly zero up to the onset.
#'
#' The drift onset is not a free estimate: it is the analytic inverse of each
#' shape at the `drift_fraction` fraction `f` of its amplitude,
#' `onset = xmid + u / k`:
#'
#' - `shape = "symmetric"`: `k = 4 * slope / (B - A)`; `u = log(f / (1 - f))`.
#' - `shape = "gompertz"`: `k = slope * e / (B - A)`; `u = -log(-log(f))`.
#' - `shape = "gompertz_left"`: `k = slope * e / (B - A)`;
#'   `u = log(-log(1 - f))`.
#'
#' The `"gompertz"` form places its onset furthest past `xmid` (slow tail)
#' and `"gompertz_left"` nearest (fast tail).
#'
#' The excursion point `texc` is where the drift rate overtakes the decaying
#' primary rate, `|S'(t)| = |slope_B|`, floored at the drift onset.
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [analyse_kinetics()], [SSsigmoidal_drift()], [logistic()],
#'   [gompertz()], [gompertz_left()], [exponential_drift()]
#'
#' @examples
#' ## create a sigmoidal curve with late linear drift and random noise
#' set.seed(13)
#' t <- 1:120
#' x <- sigmoidal_drift(
#'     t, A = 10, B = 100, xmid = 40, slope = 4,
#'     slope_B = -0.4, drift_fraction = 0.95
#' ) + rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## the drift onset fraction is held constant in the formula
#' model <- nls(
#'     x ~ SSsigmoidal_drift(
#'         t, A, B, xmid, slope, slope_B, drift_fraction = 0.95
#'     ),
#'     data = data,
#'     algorithm = "port",
#'     control = nls.control(warnOnly = TRUE)
#' )
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
sigmoidal_drift <- function(
    t,
    A,
    B,
    xmid,
    slope,
    slope_B,
    drift_fraction,
    shape = c("symmetric", "gompertz", "gompertz_left")
) {
    shape <- match.arg(shape)
    ## primary sigmoid + hinge-linear secondary drift from the onset
    S <- switch(shape, symmetric = logistic, gompertz = gompertz, gompertz_left)
    onset <- sigdrift_onset(A, B, xmid, slope, drift_fraction, shape)
    return(S(t, A, B, xmid, slope) + slope_B * pmax(t - onset, 0))
}


#' Rate constant of a sigmoidal shape
#'
#' The rate `k` such that the sigmoid of the given `shape` is a function of
#' `u = k * (t - xmid)`: `4 * slope / (B - A)` for `"symmetric"`, else
#' `slope * e / (B - A)`. Positive for a consistent fit, where `slope` and
#' `B - A` share a sign.
#'
#' @inheritParams sigmoidal_drift
#'
#' @returns A numeric rate in units of `1 / t`.
#'
#' @keywords internal
sigdrift_rate <- function(A, B, slope, shape) {
    return(slope * (if (shape == "symmetric") 4 else exp(1)) / (B - A))
}


#' Drift onset time of the sigmoidal-drift model
#'
#' The time at which a sigmoid of the given `shape` reaches the
#' `drift_fraction` fraction of its amplitude, by the analytic inverse of
#' each shape (see [sigmoidal_drift()]). Vectorised over the numeric
#' parameters; `shape` is a single string.
#'
#' @inheritParams sigmoidal_drift
#'
#' @returns A numeric vector of onset times.
#'
#' @keywords internal
sigdrift_onset <- function(A, B, xmid, slope, drift_fraction, shape) {
    ## a fraction outside (0, 1) has no onset
    if (any(drift_fraction <= 0 | drift_fraction >= 1, na.rm = TRUE)) {
        stop("`drift_fraction` must be a fraction of the amplitude in (0, 1).")
    }
    f <- drift_fraction
    u <- switch(
        shape,
        symmetric = log(f / (1 - f)),
        gompertz = -log(-log(f)),
        gompertz_left = log(-log1p(-f))
    )
    return(xmid + u / sigdrift_rate(A, B, slope, shape))
}


#' Excursion point of the sigmoidal-drift model
#'
#' The time past the inflection at which the drift rate overtakes the
#' decaying sigmoid rate, `|S'(t)| = |slope_B|`, floored at the drift onset
#' (see [sigdrift_onset()]): the turning point of the curve when the
#' phases oppose, or where the linear trend takes over a monotonic
#' response. A drift at least as fast as the peak sigmoid rate `slope`
#' takes over from the onset. Scalar parameters only.
#'
#' @inheritParams sigmoidal_drift
#'
#' @details
#' With `ratio = |slope_B / slope|` and `u = k * (t - xmid)` (see
#' [sigdrift_rate()]), the sigmoid rate relative to its peak is
#' `4 * L * (1 - L)` with `L = 1 / (1 + exp(-u))` for `"symmetric"`, solved
#' as `u = 2 * atanh(sqrt(1 - r))`; `exp(1 - u - exp(-u))` for
#' `"gompertz"`; and `exp(1 + u - exp(u))` for `"gompertz_left"`. The
#' Gompertz forms have no closed inverse and are solved by
#' [stats::uniroot()] on a bracket containing the single post-inflection
#' root.
#'
#' @returns A numeric excursion time.
#'
#' @keywords internal
sigdrift_texc <- function(A, B, xmid, slope, slope_B, drift_fraction, shape) {
    onset <- sigdrift_onset(A, B, xmid, slope, drift_fraction, shape)
    ratio <- abs(slope_B / slope)
    if (!is.finite(ratio) || ratio >= 1) {
        return(onset)
    }
    u <- switch(
        shape,
        symmetric = 2 * atanh(sqrt(1 - ratio)),
        gompertz = stats::uniroot(
            \(u) u + exp(-u) - 1 + log(ratio),
            c(0, 1 - log(ratio))
        )$root,
        gompertz_left = stats::uniroot(
            \(u) exp(u) - u - 1 + log(ratio),
            c(0, log(2 * (1 - log(ratio))))
        )$root
    )
    return(max(onset, xmid + u / sigdrift_rate(A, B, slope, shape)))
}


#' Initiate self-starting sigmoidal-drift model
#'
#' [sigdrift_init()]: Returns initial values for the parameters in a
#' `selfStart` model. The `shape` written in the model call seeds the
#' matching sigmoid (`"symmetric"` when absent).
#'
#' @inheritParams logistic_init
#'
#' @returns [sigdrift_init()]: Initial starting estimates for parameters in
#'   the model called by [SSsigmoidal_drift()].
#'
#' @keywords internal
sigdrift_init <- function(mCall, data, LHS, ...) {
    fixed <- list(...)$fixed %||% list()
    tx <- sortedXyData(mCall[["t"]], LHS, data)
    shape <- eval(mCall[["shape"]], data) %||% "symmetric"
    return(sigdrift_start(tx[["y"]], tx[["x"]], fixed, shape))
}


#' Starting estimates for the sigmoidal-drift model
#'
#' Vector-level initialiser behind [sigdrift_init()], called directly by
#' the kinetics worker on the fit window. The sigmoid is seeded as for
#' [SSgompertz()] ([init_asymptotes()], [init_inflection()]), the drift
#' onset resolved from that seed, and the residual from the seeded sigmoid
#' past the onset regressed on time from the onset: the intercept corrects
#' the asymptote `B` and the slope is the drift. A second pass re-seeds the
#' sigmoid on the drift-corrected response, correcting an inflection
#' biased by the drift. Fewer than two points past the onset seed a zero
#' drift. User-fixed values are held.
#'
#' @param x A numeric vector of the response variable (sorted by `t`).
#' @param t A numeric vector of the predictor variable.
#' @param fixed A named list of user-fixed parameter values.
#' @inheritParams sigmoidal_drift
#'
#' @returns A named numeric vector of starting estimates in model order.
#'
#' @keywords internal
sigdrift_start <- function(x, t, fixed = list(), shape = "symmetric") {
    ## `[[` throughout: `$` would partial-match `slope` to `slope_B`
    p <- fixed[["drift_fraction"]] %||% 0.95
    ab <- init_asymptotes(x)
    A <- fixed[["A"]] %||% ab$A

    refine <- \(s) {
        ## inflection from the drift-corrected response
        xd <- x - s$slope_B * pmax(t - s$onset, 0)
        infl <- init_inflection(xd, t, A, s$B)
        xmid <- fixed[["xmid"]] %||% infl$xmid
        slope <- fixed[["slope"]] %||% infl$slope
        onset <- sigdrift_onset(A, s$B, xmid, slope, p, shape)

        ## the residual from the seeded sigmoid past the onset is the
        ## asymptote error plus the drift from the onset
        i <- t >= onset
        r <- x[i] - sigmoidal_drift(t[i], A, s$B, xmid, slope, 0, p, shape)
        b <- slope(r, t[i] - onset, intercept = TRUE, bypass_checks = TRUE)
        list(
            B = fixed[["B"]] %||% (s$B + (attr(b, "intercept") %||% 0)),
            xmid = xmid,
            slope = slope,
            slope_B = fixed[["slope_B"]] %||% (if (is.na(b)) 0 else b),
            onset = onset
        )
    }

    ## an onset at the record end makes the initial drift correction zero
    s <- refine(refine(list(
        B = fixed[["B"]] %||% ab$B,
        slope_B = fixed[["slope_B"]] %||% 0,
        onset = max(t)
    )))
    return(c(
        A = A,
        B = s$B,
        xmid = s$xmid,
        slope = s$slope,
        slope_B = s$slope_B,
        drift_fraction = p
    ))
}


#' Sigmoidal-drift model with gradient
#'
#' Model function of [SSsigmoidal_drift()]: [sigmoidal_drift()] plus the
#' partial derivatives for the parameters written as bare symbols in the
#' call (see [free_params()]), so [stats::nls()] skips
#' [stats::numericDeriv()]. The sigmoid partials come from
#' [sigmoid_core()]; the drift onset `xmid + u_f / k` moves with every
#' sigmoid parameter through the rate `k`, and the hinge derivatives are
#' one-sided at the onset.
#'
#' @inheritParams sigmoidal_drift
#'
#' @returns A numeric vector of predicted values with a `"gradient"`
#'   attribute when any parameter is free.
#'
#' @keywords internal
sigdrift_model <- function(
    t,
    A,
    B,
    xmid,
    slope,
    slope_B,
    drift_fraction,
    shape = "symmetric"
) {
    g <- sigmoid_core(t, A, B, xmid, slope, shape)
    f <- drift_fraction
    ## the shape inverse at the fraction, in units of u (see sigdrift_onset)
    u_f <- switch(
        shape,
        symmetric = log(f / (1 - f)),
        gompertz = -log(-log(f)),
        gompertz_left = log(-log1p(-f))
    )
    onset <- xmid + u_f / g$k
    h <- pmax(t - onset, 0)
    val <- g$val + slope_B * h
    free <- free_params(
        match.call(),
        c("A", "B", "xmid", "slope", "slope_B", "drift_fraction")
    )
    if (length(free) > 0L) {
        ## d/dp = dS/dp - slope_B * [t > onset] * d(onset)/dp
        on <- slope_B * (t > onset)
        kD <- g$k * (B - A)
        du_f <- switch(
            shape,
            symmetric = 1 / (f * (1 - f)),
            gompertz = -1 / (f * log(f)),
            gompertz_left = -1 / ((1 - f) * log1p(-f))
        )
        grad <- cbind(
            A = g$A + on * u_f / kD,
            B = g$B - on * u_f / kD,
            xmid = g$xmid - on,
            slope = g$slope + on * u_f / (g$k * slope),
            slope_B = h,
            drift_fraction = -on * du_f / g$k
        )
        attr(val, "gradient") <- grad[, free, drop = FALSE]
    }
    return(val)
}


#' Self-starting sigmoidal-drift model
#'
#' @description
#' Creates initial coefficient estimates for a `selfStart` wrapper around
#' [sigmoidal_drift()], for use with [stats::nls()]: a 4-parameter sigmoid
#' (A, B, xmid, slope) with a linear drift `slope_B` at its ending asymptote
#' from the onset fraction `drift_fraction`.
#'
#' @usage
#' SSsigmoidal_drift(t, A, B, xmid, slope, slope_B, drift_fraction, shape)
#'
#' @inheritParams sigmoidal_drift
#'
#' @details
#' ## Model formula
#'
#' `x ~ SSsigmoidal_drift(t, A, B, xmid, slope, slope_B,
#' drift_fraction = 0.95, shape = "gompertz")`
#'
#' `drift_fraction` should be written as a constant, and `shape` is a string
#' constant (`"symmetric"` when omitted); neither is estimated. The hinge at
#' the drift onset is not differentiable, so `algorithm = "port"` with
#' `control = nls.control(warnOnly = TRUE)` is recommended.
#'
#' Starting estimates seed the sigmoid as for [SSgompertz()], resolve the
#' drift onset from that seed, and regress the residual past the onset on
#' time to seed `slope_B` and correct the asymptote `B`.
#'
#' ## Fixing parameters
#'
#' Any parameter may be held constant by writing a value in place of its name
#' in the formula, e.g.
#' `x ~ SSsigmoidal_drift(t, A = 0, B, xmid, slope, slope_B,
#' drift_fraction = 0.95)` fixes the starting asymptote at `A = 0`. Fixed
#' parameters are excluded from estimation and are not returned by
#' [stats::coef()].
#'
#' @returns A numeric vector of predicted values the same length as the
#'   predictor variable `t`.
#'
#' @seealso [sigmoidal_drift()], [analyse_kinetics()], [stats::nls()],
#'   [stats::selfStart()], [SSlogistic()], [SSgompertz()],
#'   [SSexponential_drift()]
#'
#' @examples
#' ## create a Gompertz curve with late linear drift and random noise
#' set.seed(13)
#' t <- 1:120
#' x <- sigmoidal_drift(
#'     t, A = 10, B = 100, xmid = 40, slope = 4,
#'     slope_B = -0.4, drift_fraction = 0.95, shape = "gompertz"
#' ) + rnorm(length(t), 0, 2)
#' data <- data.frame(t, x)
#'
#' ## fit with the drift onset held at 95% of the amplitude
#' model <- nls(
#'     x ~ SSsigmoidal_drift(
#'         t, A, B, xmid, slope, slope_B,
#'         drift_fraction = 0.95, shape = "gompertz"
#'     ),
#'     data = data,
#'     algorithm = "port",
#'     control = nls.control(warnOnly = TRUE)
#' )
#' summary(model)
#'
#' @export
SSsigmoidal_drift <- selfStart(
    model = sigdrift_model,
    initial = init_fixed(
        sigdrift_init,
        c("A", "B", "xmid", "slope", "slope_B", "drift_fraction")
    ),
    parameters = c("A", "B", "xmid", "slope", "slope_B", "drift_fraction")
)


#' Analyse sigmoidal-drift kinetics across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "sigmoidal_drift")`. Fits a two-phase
#' sigmoidal + linear-drift curve to each `nirs_channel` within a single
#' *"mnirs"* data frame via [fit_sigmoidal_drift()], falling back to
#' [fit_sigmoidal()] where the drift is unsupported (see
#' `kinetics_fallbacks`). See [analyse_kinetics()] for user-facing
#' documentation.
#'
#' @param drift_fraction A numeric fraction of the amplitude in `(0.5, 1)` at
#'   which the drift onset is held (*default* `0.95`). Always held constant.
#'   Applied to every channel, or per-channel as a list keyed by channel
#'   name, e.g. `drift_fraction = list(smo2 = 0.9)`.
#' @param fix An *optional* named list of model parameters (`A`, `B`,
#'   `xmid`, `slope`, `slope_B`) to hold constant during fitting, e.g.
#'   `fix = list(A = 0)`. Applied to every channel, or per-channel as a
#'   list of lists keyed by channel name, e.g.
#'   `fix = list(smo2 = list(A = 0))`.
#' @inheritParams analyse_logistic
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#' @inheritParams analyse_monoexponential
#'
#' @returns A `data.frame` with one row per `nirs_channel` and columns
#'   `nirs_channels`, `model`, `A`, `B`, `xmid`, `slope`, `texc`,
#'   `slope_B`, `drift_fraction`, `xmid_fitted`, `texc_fitted`. `texc` is
#'   the excursion point where the drift rate overtakes the decaying
#'   sigmoid rate, never before the drift onset (see [sigdrift_texc()]).
#'   Per-channel metadata are attached as attributes:
#'   - `"model"`: an [nls][stats::nls] model object, or `NULL` for channels
#'     where fitting failed.
#'   - `"fitted_data"`: a named list of per-channel data frames with
#'     columns `window_idx` and `fitted`.
#'   - `"diagnostics"`: a `data.frame` with one row per `nirs_channel`
#'     containing model fit diagnostics.
#'   - `"channel_args"`: a `data.frame` with one row per `nirs_channel`
#'     recording the resolved arguments used.
#'   - `"warnings"`: a `data.frame` of conditions captured during fitting.
#'
#' @seealso [analyse_kinetics()], [sigmoidal_drift()],
#'   [SSsigmoidal_drift()]
#'
#' @keywords internal
analyse_sigmoidal_drift <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    shape = c("symmetric", "gompertz", "gompertz_left"),
    drift_fraction = 0.95,
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
            "shape", "drift_fraction", "fix", "control", "start_time",
            "direction", "end_window"
        )),
        choices = list(
            shape = c("symmetric", "gompertz", "gompertz_left"),
            direction = c("auto", "positive", "negative")
        ),
        fix_params = c("A", "B", "xmid", "slope", "slope_B"),
        verbose = verbose,
        env = env
    )

    ## an unsupported drift falls back to the sigmoidal (see
    ## `kinetics_fallbacks`); the undocumented `model_fallback = FALSE`
    ## keeps the raw fit
    return(analyse_kinetics_channels(
        data,
        setup$nirs_channels,
        setup$time_channel,
        resolve_drift_frac(setup$per_channel, env),
        fit_sigmoidal_drift,
        verbose,
        interval_name,
        extra_args = args,
        method = "sigmoidal_drift",
        fallback = !isFALSE(args$model_fallback),
        env = env
    ))
}


#' Fit a sigmoidal-drift model to one channel
#'
#' Channel fitter of [analyse_sigmoidal_drift()] (see
#' [analyse_kinetics_channels()]). Self-starting [SSsigmoidal_drift()] of
#' the channel `shape` via [stats::nls()] with `algorithm = "port"`,
#' seeded by [sigdrift_start()], with the requested `direction` enforced
#' on `B - A` and the sign of `slope` ([enforce_direction()]).
#'
#' @inheritParams fit_monoexponential
#'
#' @returns The `coefs`/`model`/`fitted_data`/`diag` list of
#'   [build_fit_results()], or [build_na_results()] when the fit fails.
#'
#' @keywords internal
fit_sigmoidal_drift <- function(x, t, valid, .a, ctx) {
    x_fit <- x[valid$idx]
    t_fit <- t[valid$idx]
    ## NA scaffold (method columns only) for convergence failure
    na_cols <- kinetics_coef_cols$sigmoidal_drift
    params <- c("A", "B", "xmid", "slope", "slope_B", "drift_fraction")
    fn <- quote(SSsigmoidal_drift)

    ## the drift onset fraction is always held constant; the shape rides
    ## in the formula as a string constant beside the fixed parameters
    .a$fix <- c(.a$fix, list(drift_fraction = .a$drift_fraction))
    fix_all <- c(.a$fix, list(shape = .a$shape))
    free <- setdiff(params, names(.a$fix))
    ## columns carry the channel names so the model predicts on them
    nm <- fit_names(ctx$nirs, ctx$time_channel, params)
    fit_data <- list2DF(setNames(list(x_fit, t_fit), nm))
    on_error <- \(e) {
        warn_fit_failed(fn, e, ctx$nirs, ctx$interval_name, env = ctx$env)
    }

    ## the hinge is non-smooth, so port often stops short of its
    ## certificate on usable coefficients, which are kept with a warning
    model <- if (nrow(fit_data) <= length(free)) {
        on_error(simpleError(sprintf(
            "%d observation%s for %d free parameters.",
            nrow(fit_data),
            if (nrow(fit_data) == 1L) "" else "s",
            length(free)
        )))
    } else {
        formula <- build_ss_formula(
            fn,
            c(params, "shape"),
            fix_all,
            nm[[1L]],
            nm[[2L]]
        )
        tryCatch(
            {
                start <- sigdrift_start(x_fit, t_fit, .a$fix, .a$shape)
                embed_fit_call(suppressWarnings(nls(
                    formula,
                    fit_data,
                    start = start[free],
                    algorithm = "port",
                    control = fit_control(
                        .a$control,
                        maxiter = 500L,
                        warnOnly = TRUE
                    )
                )))
            },
            error = on_error
        )
    }
    model <- accept_port_fit(model, on_error)
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
        fix = fix_all,
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

    ## xmid and texc are elapsed from start_time, matching the fit
    ## time base; the coefficients are in model-argument order
    cf <- c(as.list(coefs), shape = .a$shape)
    texc_val <- do.call(sigdrift_texc, cf)
    ## predict response at xmid and texc using the full fitted model
    fitted <- do.call(sigmoidal_drift, c(list(c(cf$xmid, texc_val)), cf))

    return(build_fit_results(
        list2DF(c(
            as.list(coefs),
            list(
                texc = texc_val,
                xmid_fitted = fitted[[1L]],
                texc_fitted = fitted[[2L]]
            )
        ))[na_cols],
        model,
        x_fit,
        t_fit,
        valid,
        env = ctx$env
    ))
}
