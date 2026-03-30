#' Re-scale data range
#'
#' Expand or reduce the range (min and max values) of data channels to a new
#' amplitude/dynamic range, e.g. re-scale the range of NIRS data to `c(0, 100)`.
#'
#' @usage
#' rescale_mnirs(
#'   data,
#'   nirs_channels = list(NULL),
#'   range,
#'   verbose = TRUE
#' )
#' 
#' @param nirs_channels A `list()` of character vectors indicating grouping
#'   structure of mNIRS channel names to operate on (see *Details*). Must
#'   match column names in `data` exactly. Retrieved from metadata if not
#'   defined explicitly.
#'   \describe{
#'      \item{`list("A", "B", "C")`}{Will operate on each channel independently,
#'      losing the relative scaling between channels.}
#'      \item{`list(c("A", "B", "C"))`}{Will operate on all channels together,
#'      preserving the relative scaling between channels.}
#'      \item{`list(c("A", "B"), c("C", "D"))`}{Will operate on channels `A`
#'      & `B` in one group, and `C` & `D` in another group, preserving
#'      relative scaling within, but not between groups.}
#'   }
#' @param range A numeric vector in the form `c(min, max)`, indicating the
#'   range of output values to which data channels will be re-scaled.
#' @inheritParams validate_mnirs
#'
#' @details
#' `nirs_channels = list()` can be used to group data channels (column names)
#'   to preserve absolute or relative scaling.
#'
#' - Channels grouped together in a vector (e.g. `list(c("A", "B"))`) will be
#'   re-scaled to a common range, and the relative scaling within that group
#'   will be preserved.
#'
#' - Channels in separate list vectors (e.g. `list("A", "B")`) will be
#'   re-scaled independently, and relative scaling between groups will be lost.
#'
#' - A single vector of channel names (e.g. `c("A", "B")`) will group
#'   channels together.
#'
#' - Channels (columns) in `data` not explicitly defined in `nirs_channels`
#'   will be passed through untouched to the output data frame.
#'
#' `nirs_channels` can be retrieved automatically from `data` of class
#'   *"mnirs"* which has been processed with `{mnirs}`, if not defined
#'   explicitly. This will default to returning all `nirs_channels` grouped
#'   together, and should be defined explicitly for other grouping arrangements.
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
#'     rescale_mnirs(        ## un-grouped nirs channels to rescale separately 
#'         nirs_channels = list(smo2_left, smo2_right), 
#'         range = c(0, 100) ## rescale to a 0-100% functional exercise range
#'     )
#' 
#' data
#'
#' \donttest{
#'     if (requireNamespace("ggplot2", quietly = TRUE)) {
#'         plot(data, time_labels = TRUE) +
#'             ggplot2::geom_hline(yintercept = c(0, 100), linetype = "dotted")
#'     }
#' }
#'
#' @export
rescale_mnirs <- function(
    data,
    nirs_channels = list(NULL),
    range,
    verbose = TRUE
) {
    ## validate =================================
    validate_mnirs_data(data, ncol = 1)
    metadata <- attributes(data)
    if (missing(verbose)) {
        verbose <- getOption("mnirs.verbose", default = TRUE)
    }
    nirs_channels <- validate_nirs_channels(
        enquo(nirs_channels), data, verbose
    )
    validate_numeric(
        range, 2, 
        msg1 = "two-element", 
        msg2 = "between {col_blue('range[1], range[2]]')}."
    )

    ## rescale range ================================
    ## this is actually a fast modify-in-place for loop
    for (group in nirs_channels) {
        group_data <- as.matrix(data[, group, drop = FALSE])
        group_data_range <- range(group_data, na.rm = TRUE)
        range_diff <- diff(group_data_range)

        if (range_diff != 0) {
            # fmt: skip
            data[, group] <- (group_data - group_data_range[1]) /
                range_diff * diff(range) + range[1]
        }
    }

    ## Metadata =================================
    metadata$nirs_channels <- unique(
        c(metadata$nirs_channels, unlist(nirs_channels))
    )

    return(create_mnirs_data(data, metadata))
}
