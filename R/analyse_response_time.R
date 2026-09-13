#' Fractional response time
#'
#' @description
#' Estimate the time at which a numeric vector reaches a specified fraction
#' of its total response amplitude relative to a baseline, e.g.
#' *half-response time* at `response_fraction = 0.5`. Vector-level companion
#' to [analyse_kinetics()] with `method = "response_time"`.
#'
#' @param start_time A numeric value in units of `t` specifying the response
#'   onset. Samples where `t <= start_time` define the baseline window.
#'   *Default* is `0`.
#' @param response_fraction A numeric vector in the range `[0, 1]` specifying
#'   the fractional response amplitude(s) to detect. Defaults to `0.5` (50%
#'   response, i.e. half-response time). Multiple values (e.g.
#'   `c(0.5, 0.632)`) return one element per fraction.
#' @param ... Additional arguments.
#' @inheritParams replace_invalid
#' @inheritParams find_kinetics_idx
#' @inheritParams validate_mnirs
#'
#' @details
#' A non-parametric approach (estimated directly from the observed data without
#' assuming a specific mathematical shape). `response_fraction = 0.5`
#' approximates the inflection point (`xmid`) of a symmetric sigmoid function.
#' `response_fraction = 0.632` approximates the time constant (`tau`;
#' \eqn{\tau}) of a monoexponential function, or `xmid` of a left-Gompertz
#' function. `response_fraction = 0.368` approximates `xmid` of a
#' right-Gompertz function. This is a good fallback estimation method if
#' parametric methods are not successfully fit.
#'
#' ## Method
#'
#' The target response value is: `fitted = A + (B - A) * response_fraction`
#'
#' Where `A` is the mean baseline value (`t <= start_time`) and `B` is the
#' extreme (peak or trough) value after `start_time`. `response_value` is the
#' first observed sample equal to or greater/lesser than the target `fitted`
#' value (above for *"positive"*, below for *"negative"* `direction`).
#' `response_time` is the elapsed time from `start_time` to `response_value`.
#'
#' [analyse_kinetics()] first trims `x` to `end_window` past the first extreme,
#' so `B` there is the first local extreme with no greater/lesser values
#' within `end_window`. Called directly, `B` is the global extreme of `x`
#' after `start_time`.
#'
#' ## Direction
#'
#' `direction` is detected automatically by default as either *"positive"*
#' (upward) or *"negative"* (downward) response, from the dominant excursion
#' of `x` above or below its initial baseline (the median of the earliest
#' samples). When tied, the greater absolute extreme decides. `B` is the
#' maximum for *"positive"* or the minimum for *"negative"*, and can be
#' overwritten manually.
#'
#' ## Baseline
#'
#' When no samples exist where `t <= start_time`, the first sample `x[1]` is
#' used as the baseline `A` with a warning. `start_time` must be within the
#' range of `t`.
#'
#' @returns A named list containing:
#'   \item{`A`}{The mean baseline value of `x` where `t <= start_time`.}
#'   \item{`B`}{The extreme (maximum or minimum) value of `x` after
#'   `start_time`.}
#'   \item{`response_time`}{The elapsed time from `start_time` to the
#'   fractional response, in units of `t`; one element per
#'   `response_fraction`.}
#'   \item{`response_value`}{The observed value of `x` at the response index;
#'   one element per `response_fraction`.}
#'   \item{`fitted`}{The target fractional response value
#'   `A + (B - A) * response_fraction`; one element per `response_fraction`.}
#'   \item{`baseline_idx`}{Integer indices where `t <= start_time`.}
#'   \item{`response_idx`}{Integer index at each `response_value`.}
#'   \item{`extreme_idx`}{Integer index at the extreme value `B`.}
#'
#' @seealso [analyse_kinetics()], [peak_slope()], [monoexponential()]
#'
#' @examples
#' ## create an exponential curve with random noise
#' set.seed(13)
#' t <- 0:60
#' x <- monoexponential(t, A = 20, B = 60, tau = 8, TD = 10) +
#'     rnorm(length(t), 0, 1)
#'
#' ## half-response time (0.5) and time constant approximation (0.632 ~= tau)
#' RT <- response_time(x, t, start_time = 10, response_fraction = c(0.5, 0.632))
#' RT$response_time
#'
#' plot(t, x, type = "l", col = "grey60", xlab = "t", ylab = "x")
#' ## mean baseline `A` across the baseline window
#' segments(
#'     t[min(RT$baseline_idx)], RT$A,
#'     t[max(RT$baseline_idx)], RT$A,
#'     col = "red", lwd = 2
#' )
#' ## response values at 0.5 (red) and 0.632 (blue), and the extreme `B`
#' points(
#'     t[RT$response_idx],
#'     RT$response_value,
#'     col = c("red", "blue"),
#'     pch = 19
#' )
#' points(t[RT$extreme_idx], RT$B, col = "red", pch = 19)
#'
#' @export
response_time <- function(
    x,
    t = seq_along(x),
    start_time = 0,
    response_fraction = 0.5,
    direction = c("auto", "positive", "negative"),
    verbose = TRUE,
    ...
) {
    ## internal callers pass `env` through `...` to report conditions
    ## as coming from the user-facing function
    env <- list(...)$env %||% environment()
    validate_numeric(
        response_fraction, Inf, c(0, 1), msg2 = "between {col_blue('[0, 1]')}.",
        env = env
    )
    direction <- match.arg(direction)
    args <- list(...)

    if (!(args$bypass_checks %||% FALSE)) {
        validate_x_t(x, t, allow_na = TRUE, env = env)
        if (missing(verbose)) {
            verbose <- getOption("mnirs.verbose", default = TRUE)
        }
        ## detect direction from net trend, fallback to abs magnitude
        direction <- detect_direction(x, t, direction = direction)
    }

    baseline_idx <- which(t <= start_time)

    if (!(args$bypass_checks %||% FALSE)) {
        validate_numeric(start_time, 1L, env = env)
        if (length(baseline_idx) == 0L) {
            if (verbose) {
                cli_warn(c(
                    "!" = "No observations where {.arg t} <= \\
                    {.arg start_time} = {.val {start_time}}.",
                    "i" = "{.code x[1]} used as response baseline."
                ), call = warn_call(env))
            }
            baseline_idx <- 1L
            start_time <- t[baseline_idx]
        }
        if (start_time > t[length(t)]) {
            cli_abort(c(
                "x" = "No observations in {.arg t} before {.arg start_time}.",
                "i" = "{.arg start_time} must be specified within the \\
                range of {.arg t}."
            ), call = env)
        }
    }

    ## process =====================================================
    ## look for extreme after start_time
    x_valid <- c(rep(NA_real_, length(baseline_idx)), x[t > start_time])
    extreme_idx <- if (direction == "positive") {
        which.max(x_valid)
    } else {
        which.min(x_valid)
    }
    ## no valid samples after start_time: which.max/min return integer(0),
    ## which would collapse B and downstream 1-row results to length 0
    if (length(extreme_idx) == 0L) {
        extreme_idx <- NA_integer_
    }

    A <- mean(x[baseline_idx], na.rm = TRUE)
    B <- x[extreme_idx]
    response_fitted <- A + (B - A) * response_fraction
    compare_fn <- if (direction == "positive") `>=` else `<=`
    ## first sample reaching each fractional response value
    response_idx <- vapply(response_fitted, \(.f) {
        which(compare_fn(x_valid, .f))[1L]
    }, integer(1))

    if (anyNA(response_idx)) {
        if (verbose) {
            cli_warn(c(
                "!" = "No valid {.val {direction}} extremes after \\
                {.arg start_time}. Returning {.val {NA}}."
            ), call = warn_call(env))
        }
        response_fitted[is.na(response_idx)] <- NA_real_
    }

    return(list(
        A = A,
        B = B,
        response_time = t[response_idx] - start_time, ## real
        response_value = x[response_idx], ## real
        fitted = response_fitted, ## predicted
        baseline_idx = baseline_idx, ## all baseline samples
        response_idx = response_idx, ## mid sample
        extreme_idx = extreme_idx ## end sample
    ))
}


#' Analyse fractional kinetics response time across NIRS channels
#'
#' Internal channel-level dispatch for
#' `analyse_kinetics(method = "response_time")`. Computes the fractional
#' response time for each `nirs_channel` within a single *"mnirs"* data
#' frame. See [analyse_kinetics()] for user-facing documentation.
#'
#' @inheritParams validate_mnirs
#' @inheritParams analyse_kinetics
#' @inheritParams response_time
#'
#' @returns A `data.frame` with one row per `nirs_channel` per
#'   `response_fraction` and columns `nirs_channels`, `response_fraction`,
#'   `A`, `B`, `response_time`, `response_value`, `fitted`, `idx`. Per-channel
#'   metadata are attached as attributes:
#'   - `"model"`: `NULL` (no parametric model is fitted).
#'   - `"fitted_data"`: a named list of per-channel data frames with
#'     columns `window_idx` and `fitted`, containing the baseline,
#'     response, and extreme key points.
#'   - `"diagnostics"`: a `data.frame` with one row per `nirs_channel`
#'     containing model fit diagnostics.
#'   - `"channel_args"`: a `data.frame` with one row per `nirs_channel`
#'     recording the resolved arguments used.
#'
#' @seealso [analyse_kinetics()], [response_time()]
#'
#' @keywords internal
analyse_response_time <- function(
    data,
    nirs_channels = NULL,
    time_channel = NULL,
    start_time = NULL,
    response_fraction = 0.5,
    direction = c("auto", "positive", "negative"),
    end_window = Inf,
    verbose = TRUE,
    ...,
    env = rlang::caller_env()
) {
    ## validation ==============================================
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    args <- list(...)
    ## interval label; falls back to the `data` argument name when unsupplied
    interval_name <- args$interval_name %||% deparse(substitute(data))

    ## shared prologue: validate data, resolve channels/time, broadcast and
    ## validate per-channel args
    setup <- setup_kinetics_worker(
        data,
        enquo(nirs_channels),
        enquo(time_channel),
        arg_list = mget(
            c("start_time", "response_fraction", "direction", "end_window")
        ),
        choices = list(direction = c("auto", "positive", "negative")),
        verbose = verbose,
        env = env
    )
    ## method-specific fit: fractional response time (no model fit)
    response_time_fit <- function(x, t, valid, .a, ctx) {
        x_fit <- x[valid$idx]
        t_fit <- t[valid$idx]
        ## quote = TRUE so `env` (a defused call object for condition
        ## attribution) is passed as-is, not evaluated by do.call
        ## `t_fit` is elapsed from start_time, so the baseline splits at 0.
        ## `.a` itself is left intact for `channel_args` reporting.
        ## verbose = TRUE so fit warnings always signal; the capture handler
        ## in analyse_kinetics_channels() governs console emission
        response <- do.call(response_time, c(
            list(x = x_fit, t = t_fit),
            replace(.a, "start_time", 0),
            list(verbose = TRUE, bypass_checks = TRUE, env = env),
            args
        ), quote = TRUE)

        ## `response_time()` indexes the fit window; map the response sample
        ## back to its original data frame row, which differs whenever
        ## non-finite samples were dropped. NA idx maps to NA
        ## one row per response_fraction; scalar A/B recycle across rows
        coefs <- data.frame(
            response_fraction = .a$response_fraction,
            A = response$A,
            B = response$B,
            response_time = response$response_time,
            response_value = response$response_value,
            fitted = response$fitted,
            idx = valid$idx[response$response_idx]
        )

        ## bind baseline vec with `A`, and response and extreme scalars,
        ## mapping fit-window positions back to original data frame rows
        fitted_data <- data.frame(
            window_idx = valid$idx[c(
                response$baseline_idx,
                response$response_idx,
                response$extreme_idx
            )],
            fitted = c(
                rep(response$A, length(response$baseline_idx)),
                response$fitted,
                response$B
            )
        )
        ## omit NA idx
        fitted_data <- fitted_data[!is.na(fitted_data$window_idx), ]

        list(
            coefs = coefs,
            model = NULL,
            fitted_data = fitted_data,
            diag = compute_diagnostics(
                x = x_fit[1L:3L], ## placeholder
                t = t_fit[1L:3L], ## placeholder
                fitted = c(coefs$A[1L], coefs$fitted[1L], coefs$B[1L]),
                n_params = 0L, ## invalid for response time method
                env = env
            )
        )
    }

    return(analyse_kinetics_channels(
        data,
        setup$nirs_channels,
        setup$time_channel,
        setup$per_channel,
        response_time_fit,
        verbose,
        interval_name,
        extra_args = args,
        env = env
    ))
}
