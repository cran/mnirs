#' Shift data range
#'
#' Move the range of data channels in a data frame up or down, while preserving
#' the absolute amplitude/dynamic range of each channel, and the relative
#' scaling across channels. e.g. shift the minimum data value to zero for all
#' positive values, or shift the mean of the first time span in a recording
#' to zero.
#'
#' @param to A numeric value in units of `nirs_channels` to which the data
#'   channels will be shifted, e.g. shift the minimum value to zero.
#' @param by A numeric value in units of `nirs_channels` by which the data
#'   channels will be shifted, e.g. shift all values up by 10 units.
#' @param position Indicates where the reference values will be shifted from.
#'   \describe{
#'      \item{`"min"`}{(The *default*) will shift the minimum value(s) `to`
#'      or `by` the specified value.}
#'      \item{`"max"`}{Will shift the maximum value(s) `to` or `by` the
#'      specified values.}
#'      \item{`"first"`}{Will shift first value(s) `to` or `by` the specified
#'      values.}
#'   }
#' @inheritParams validate_mnirs
#' @inheritParams replace_mnirs
#' @inheritParams rescale_mnirs
#'
#' @details
#' `nirs_channels = list()` can be used to group data channels (column names)
#'   to preserve absolute or relative scaling.
#'
#' - Channels grouped together in a vector (e.g. `list(c("A", "B"))`) will be
#'   shifted to a common value, and the relative scaling within that group
#'   will be preserved.
#'
#' - Channels in separate list vectors (e.g. `list("A", "B")`) will be
#'   shifted independently, and relative scaling between groups will be lost.
#'
#' - A single vector of channel names (e.g. `c("A", "B")`) will group
#'   channels together.
#'
#' - Channels (columns) in `data` not explicitly defined in `nirs_channels`
#'   will be passed through untouched to the output data frame.
#'
#' `nirs_channels` and `time_channel` can be retrieved automatically from
#'   `data` of class *"mnirs"* which has been processed with `{mnirs}`,
#'   if not defined explicitly. This will default to returning all
#'   `nirs_channels` grouped together, and should be defined explicitly
#'   for other grouping arrangements.
#'
#' Only one of either `to` or `by` and one of either `width` or `span` should
#'   be defined. If both of either pairing are defined, `to` will be preferred
#'   over `by`, and `width` will be preferred over `span`.
#'
#' @returns
#' A [tibble][tibble::tibble-package] of class *"mnirs"* with metadata
#'   available with `attributes()`.
#'
#' @examples
#' ## read example data
#' data <- read_mnirs(
#'     file_path = example_mnirs("moxy_ramp"),
#'     nirs_channels = c(smo2_left = "SmO2 Live",
#'                       smo2_right = "SmO2 Live(2)"),
#'     time_channel = c(time = "hh:mm:ss"),
#'     verbose = FALSE
#' ) |>
#'     shift_mnirs(        ## un-grouped nirs channels to shift separately 
#'         nirs_channels = list(smo2_left, smo2_right), 
#'         to = 0,         ## NIRS values will be shifted to zero
#'         span = 120,     ## shift the *first* 120 sec of data to zero
#'         position = "first"
#'     )
#'
#' data
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         plot(data, time_labels = TRUE) +
#'             ggplot2::geom_hline(yintercept = 0, linetype = "dotted")
#'     }
#' }  
#'
#' @export
shift_mnirs <- function(
    data,
    nirs_channels = list(),
    time_channel = NULL,
    to = NULL,
    by = NULL,
    width = NULL,
    span = NULL,
    position = c("min", "max", "first"),
    verbose = TRUE
) {
    ## validation =============================================
    ## do nothing condition
    if (is.null(c(to, by))) {
        cli_abort(c(
            "x" = "Shift value undefined",
            "i" = "One of {.arg to} or {.arg by} must be defined."
        ))
    }

    validate_mnirs_data(data)
    metadata <- attributes(data)
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    nirs_channels <- validate_nirs_channels(
        enquo(nirs_channels), data, verbose, as_list = TRUE
    )
    time_channel <- validate_time_channel(enquo(time_channel), data)
    validate_numeric(to, 1, msg1 = "one-element")
    validate_numeric(by, 1, msg1 = "one-element")
    if (!is.null(to) && !is.null(by)) {
        by <- NULL
        if (verbose) {
            cli_inform(c("i" = "{.arg to} = {.val {to}} overrides {.arg by}."))
        }
    }

    if (
        verbose &&
            is.null(attr(data, "nirs_channels")) &&
            !is.list(nirs_channels)
    ) {
        cli_inform(c(
            "!" = "{.fn shift_mnirs} accepts {.arg nirs_channels} = \\
            {col_blue('list()')} for channel grouping. See `?shift_mnirs`."
        ))
    }
    nirs_listed <- make_list(nirs_channels)
    nirs_unlisted <- unlist(nirs_listed, use.names = FALSE)

    ## Metadata =================================
    metadata$nirs_channels <- unique(nirs_unlisted)
    metadata$time_channel <- time_channel

    ## shift_by ====================================================
    ## shift_by does not require grouping or calculating positions
    if (!is.null(by) && is.null(to)) {
        data[nirs_unlisted] <- data[nirs_unlisted] + by
        return(create_mnirs_data(data, metadata))
    }

    ## calculate shift_to values ====================================
    ## validate
    position <- match.arg(position)
    validate_width_span(width, span, verbose, "for `shift_mnirs()`.")
    time_vec <- data[[time_channel]]

    if (position == "first") {
        ## for span, take data <= first time_channel value + span, assuming sorted
        width <- width %||% rev(which(time_vec <= (time_vec[1L] + span)))[1L]
        ## drop = FALSE to avoid reducing to vector with one `nirs_unlisted`
        shift_values <- colMeans(
            data[seq_len(width), nirs_unlisted, drop = FALSE],
            na.rm = TRUE
        )
    } else if (position %in% c("min", "max")) {
        ## find local windows within width/span centred around idx
        window_idx <- compute_local_windows(
            t = time_vec, width = width, span = span,
        )
        shift_fun <- match.fun(position)
        ## compute min or max along local means
        ## return named vec of min/max for each nirs_channel
        shift_values <- vapply(data[nirs_unlisted], \(.x) {
            shift_fun(
                compute_local_fun(.x, window_idx, mean, na.rm = TRUE), 
                na.rm = TRUE
            )
        }, numeric(1))
    }

    ## apply shifts ==============================================
    ## find the value to shift to per nirs_channel group
    group_shift_values <- vapply(nirs_listed, \(.cols) {
        switch(
            position,
            min = min(shift_values[.cols]),
            max = max(shift_values[.cols]),
            first = mean(shift_values[.cols])
        )
    }, numeric(1))

    ## expand grouped shift values for each nirs_channel
    channel_shifts <- rep(group_shift_values, lengths(nirs_listed))
    ## move the shift_value to the `to` value for each nirs_channel
    data[nirs_unlisted] <- lapply(seq_along(nirs_unlisted), \(.i) {
        data[[nirs_unlisted[.i]]] - channel_shifts[.i] + to
    })

    return(create_mnirs_data(data, metadata))
}
