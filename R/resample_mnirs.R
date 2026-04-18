#' Re-sample an *{mnirs}* data frame
#'
#' Up- or down-sample an *"mnirs"* data frame to a new sample rate, filling
#' new samples via nearest-neighbour matching or interpolation.
#'
#' @param resample_rate An *optional* sample rate (Hz) for the output data
#'   frame. If `NULL` (*default*) resamples to the existing `sample_rate`,
#'  which regularises any irregular samples without changing the rate.
#'
#' @param method A character string specifying how new samples are filled.
#'   Default is *"none"*. Filling must be opted into explicitly (see *Details*):
#'   \describe{
#'     \item{`"none"`}{Matches each new sample to the nearest original
#'     `time_channel` value without any interpolation, to within tolerance of
#'     half a sample-interval. New samples are returned as `NA`.}
#'     \item{`"locf"`}{(*"Last observation carried forward"*). Fills new and
#'     missing samples with the most recent valid non-`NA` value to the left,
#'     or the nearest valid value to the right for leading `NA`s. Safe for
#'     numeric, integer, and character columns.}
#'     \item{`"linear"`}{Fills new and missing samples via linear interpolation
#'     using [stats::approx()]. Suitable for numeric columns only; non-numeric
#'     columns will fall back to `"locf"` behaviour.}
#'   }
#'
#' @inheritParams validate_mnirs
#'
#' @details
#' This function uses [replace_missing()] (based on [stats::approx()]) to
#' interpolate across new samples in the resampled data range.
#'
#' ## Sample rate and time channel
#'
#' `time_channel` and `sample_rate` are retrieved automatically from `data`
#' of class *"mnirs"*, if not defined explicitly.
#'
#' Otherwise, `sample_rate` will be estimated from the values in `time_channel`.
#' However, this may return unexpected values, and it is safer to define
#' `sample_rate` explicitly or retrieve it from *"mnirs"* metadata.
#'
#' ## Default behaviour
#'
#' When `resample_rate` is omitted, the output has the same `sample_rate` as
#' the input but with a regular, evenly-spaced `time_channel`. This is useful
#' for regularising data that contains missing or repeated samples without
#' changing the nominal rate.
#'
#' ## Column handling
#'
#' Numeric columns are interpolated according to `method` (see
#' `?replace_missing`). Non-numeric columns (e.g. character event labels,
#' integer lap numbers) are always filled by last-observation-carried-forward,
#' regardless of `method`:
#'
#' - For `method = "none"`, existing rows are matched to the nearest original
#'   values of `time_channel` without interpolation or filling, meaning newly
#'   created samples and any `NA`s in the original data are returned as `NA`.
#' - When down-sampling, numeric columns use time-weighted averaging. 
#'   Non-numeric columns use the first valid value in each output bin.
#'
#' @returns A [tibble][tibble::tibble-package] of class `"mnirs"`. Metadata are
#'   stored as attributes and can be accessed with `attributes(data)`.
#'
#' @examples
#' ## read example data
#' data <- read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"),
#'     nirs_channels = c(smo2 = "SmO2 Live"),
#'     time_channel = c(time = "hh:mm:ss"),
#'     verbose = TRUE
#' )
#'
#' ## note warning about irregular sampling
#' data
#'
#' data_resampled <- resample_mnirs(
#'     data,               ## blank channels will be retrieved from metadata
#'     resample_rate = 2,  ## blank by default will resample to `sample_rate`
#'     method = "linear",  ## linear interpolation across resampled indices
#'     verbose = TRUE      
#' )
#'
#' ## note the altered `time` values resolving the above warning
#' data_resampled
#'
#' @export
resample_mnirs <- function(
    data,
    time_channel = NULL,
    sample_rate = NULL,
    resample_rate = sample_rate, ## placeholder indicating default condition
    method = c("none", "linear", "locf"),
    verbose = TRUE
) {
    ## validation ====================================
    validate_mnirs_data(data)
    metadata <- attributes(data)
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    time_channel <- validate_time_channel(enquo(time_channel), data)
    sample_rate <- validate_sample_rate(
        data, time_channel, sample_rate, verbose
    )
    validate_numeric(
        resample_rate, 1, c(0, Inf), FALSE, msg1 = "one-element positive"
    )
    ## assign default resample_rate as sample_rate
    if (is.null(resample_rate)) {
        resample_rate <- sample_rate
    }
    resample_time <- 1 / resample_rate
    method <- match.arg(method)

    ## calculate resampling parameters  ===================================
    colnames <- names(data)
    time_vec <- data[[time_channel]]
    time_range <- floor(range(time_vec, na.rm = TRUE) * sample_rate) /
        sample_rate
    resampled_times <- seq(time_range[1L], time_range[2L], by = resample_time)
    result <- setNames(data.frame(resampled_times), time_channel)

    ## identify numeric cols; integers and `time_channel` excluded
    numeric_cols <- vapply(data, \(.x) {
        is.numeric(.x) && !is.integer(.x)
    }, logical(1))
    numeric_cols[time_channel] <- FALSE
    non_numeric_cols <- colnames[!numeric_cols & colnames != time_channel]

    ## nearest-match index for "none" method =============================
    if (method == "none") {
        idx <- max.col(
            -abs(outer(resampled_times, time_vec, `-`)),
            ties.method = "first"
        )
        tol <- resample_time * 0.5
        match_ok <- abs(resampled_times - time_vec[idx]) < tol
        idx[!match_ok] <- NA
    }

    ## interpolate numeric columns (locf/linear only) ==================
    if (method != "none" && any(numeric_cols)) {
        result[colnames[numeric_cols]] <- lapply(
            data[numeric_cols], \(.x) {
                replace_missing(
                    .x, time_vec, method = method, xout = resampled_times
                )
            }
        )
    }

    ## index-based fill for non-numeric columns ===========================
    idx_cols <- if (method == "none") {
        colnames[colnames != time_channel]
    } else {
        non_numeric_cols
    }

    if (length(idx_cols) > 0) {
        if (method == "none") {
            ## tolerance-matched index computed above
            result[idx_cols] <- lapply(data[idx_cols], \(.x) .x[idx])
        } else if (resample_rate < sample_rate) {
            ## down-sample: first non-NA per bin
            bin <- findInterval(
                time_vec, resampled_times, rightmost.closed = FALSE
            )
            bin[bin == 0] <- 1
            result[idx_cols] <- lapply(data[idx_cols], \(.x) {
                unname(unsplit(
                    lapply(split(.x, bin), \(.v) .v[!is.na(.v)][1L]),
                    unique(bin)
                ))
            })
        } else {
            ## up-sample / regularise: forward fill
            idx <- findInterval(
                resampled_times, time_vec, rightmost.closed = FALSE
            )
            idx[idx == 0] <- 1
            result[idx_cols] <- lapply(data[idx_cols], \(.x) .x[idx])
        }
    }

    ## Metadata =================================
    metadata$time_channel <- time_channel
    metadata$sample_rate <- resample_rate

    if (verbose) {
        cli_inform(c(
            "i" = "Output is resampled at {.val {resample_rate}} Hz."
        ))
    }

    ## column order same as input
    result <- result[colnames]
    return(create_mnirs_data(result, metadata))
}
