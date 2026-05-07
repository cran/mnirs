## -----------------------------------------------------------------------------
#| label: setup
#| warning: false

# pak::pak("jemarnold/mnirs") ## install development version
library(ggplot2)   ## for plotting
library(mnirs)


## -----------------------------------------------------------------------------
nirs_channels = c(
    renamed1 = "original_name1", 
    renamed2 = "original_name2"
)


## -----------------------------------------------------------------------------
## {mnirs} includes sample files from a few NIRS devices
example_mnirs()

## partial matching will error if matches multiple
try(example_mnirs("moxy"))

data_raw <- read_mnirs(
    file_path = example_mnirs("moxy_ramp"), ## call an example data file
    nirs_channels = c(
        smo2_left = "SmO2 Live",    ## identify and rename channels
        smo2_right = "SmO2 Live(2)"
    ),
    time_channel = c(time = "hh:mm:ss"), ## date-time format will be converted to numeric
    event_channel = NULL,           ## leave blank if unused
    sample_rate = NULL,             ## if blank, will be estimated from time_channel
    add_timestamp = FALSE,          ## omit a date-time timestamp column
    zero_time = TRUE,               ## recalculate time values from zero
    keep_all = FALSE,               ## return only the specified data channels
    verbose = TRUE                  ## show warnings & messages
)

## Note the above info message that sample_rate was estimated correctly at 2 Hz ☝
## ignore the warnings about irregular sampling for now, we will resample later

data_raw


## -----------------------------------------------------------------------------
## note the `time_labels` plot argument to display time values as `h:mm:ss`
plot(
    data_raw,
    points = FALSE,
    time_labels = TRUE,
    n.breaks = 5,
    na.omit = FALSE
)


## -----------------------------------------------------------------------------
## view metadata
attr(data_raw, "nirs_channels")

attr(data_raw, "time_channel")

attr(data_raw, "sample_rate")


## -----------------------------------------------------------------------------
data_resampled <- resample_mnirs(
    data_raw,            ## blank channels will be retrieved from metadata
    time_channel = time, ## channels can be left blank or specified explicitly
    sample_rate = NULL,  ## blank by default will be retrieved from metadata
    resample_rate = 2,   ## blank by default will resample to `sample_rate`
    method = "linear"    ## linear interpolation across resampled indices
)

## note the altered "time" values from the original data frame 👇
data_resampled


## -----------------------------------------------------------------------------
data_cleaned <- replace_mnirs(
    data_resampled,     ## blank channels will be retrieved from metadata
    invalid_values = 0, ## known invalid values in the data
    invalid_above = 90, ## remove data spikes above 90
    outlier_cutoff = 3, ## recommended default value
    width = 7,          ## window to detect and replace outliers/missing values
    method = "linear"   ## linear interpolation over `NA`s
)

plot(data_cleaned, time_labels = TRUE)


## -----------------------------------------------------------------------------
data_filtered <- filter_mnirs(
    data_cleaned,           ## blank channels will be retrieved from metadata
    method = "butterworth", ## Butterworth digital filter is a common choice
    order = 2,              ## filter order number
    W = 0.02,               ## filter fractional critical frequency `[0, 1]`
    type = "low",           ## specify a "low-pass" filter
    na.rm = TRUE            ## explicitly ignore NAs
)

## we will add the non-filtered data back to the plot to compare
plot(data_filtered, time_labels = TRUE) +
    geom_line(
        data = data_cleaned, 
        aes(y = smo2_left, colour = "smo2_left"), alpha = 0.4
    ) +
    geom_line(
        data = data_cleaned, 
        aes(y = smo2_right, colour = "smo2_right"), alpha = 0.4
    )


## -----------------------------------------------------------------------------
data_shifted <- shift_mnirs(
    data_filtered,     ## un-grouped nirs channels to shift separately 
    nirs_channels = list(smo2_left, smo2_right), ## 👈 channels grouped separately
    to = 0,            ## NIRS values will be shifted to zero
    span = 120,        ## shift the *first* 120 sec of data to zero
    position = "first"
)

plot(data_shifted, time_labels = TRUE) +
    geom_hline(yintercept = 0, linetype = "dotted")


## -----------------------------------------------------------------------------
data_rescaled <- rescale_mnirs(
    data_filtered,    ## un-grouped nirs channels to rescale separately 
    nirs_channels = list(smo2_left, smo2_right), ## 👈 channels grouped separately
    range = c(0, 100) ## rescale to a 0-100% functional exercise range
)

plot(data_rescaled, time_labels = TRUE) +
    geom_hline(yintercept = c(0, 100), linetype = "dotted")


## -----------------------------------------------------------------------------
options(mnirs.verbose = FALSE)

nirs_data <- read_mnirs(
    example_mnirs("train.red"),
    nirs_channels = c(
        smo2_left = "SmO2 unfiltered",
        smo2_right = "SmO2 unfiltered"
    ),
    time_channel = c(time = "Timestamp (seconds passed)"),
    zero_time = TRUE
) |>
    resample_mnirs(method = "linear") |> ## default settings will resample to the same `sample_rate`
    replace_mnirs(
        invalid_above = 73,
        outlier_cutoff = 3,
        span = 7
    ) |>
    filter_mnirs(
        method = "butterworth",
        order = 2,
        W = 0.02,
        na.rm = TRUE
    ) |>
    shift_mnirs(
        nirs_channels = list(smo2_left, smo2_right), ## 👈 channels grouped separately
        to = 0,
        span = 60,
        position = "first"
    ) |>
    rescale_mnirs(
        nirs_channels = list(c(smo2_left, smo2_right)), ## 👈 channels grouped together
        range = c(0, 100)
    )

plot(nirs_data, time_labels = TRUE)


## -----------------------------------------------------------------------------
#| fig.width: 8
#| fig.asp: 0.55

## return each interval independently with `event_groups = "distinct"`
distinct <- extract_intervals(
    nirs_data,                  ## channels blank for "distinct" grouping
    start = by_time(348, 1064), ## manually identified interval start times
    end = by_time(458, 1174),   ## interval end time (start + 150 sec)
    event_groups = "distinct",  ## return a list of data frames for each (2) event
    span = c(0, 0),             ## no boundary modification
    zero_time = FALSE           ## return original time values
)

plot(distinct, time_labels = TRUE)


## -----------------------------------------------------------------------------
## ensemble average both intervals with `event_groups = "ensemble"`
ensemble <- extract_intervals(
    nirs_data,                  ## channels recycled to all intervals by default
    nirs_channels = c(smo2_left, smo2_right),
    start = by_time(368, 1084), ## alternatively specify start times + span
    event_groups = "ensemble",  ## ensemble-average across two intervals
    span = c(-20, 90),          ## span recycled to all intervals by default
    zero_time = TRUE            ## re-calculate common time to start from `0`
)

plot(ensemble, time_labels = TRUE) + 
    geom_vline(xintercept = 0, linetype = "dotted")

