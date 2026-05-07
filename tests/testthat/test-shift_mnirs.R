test_that("shift_mnirs requires either to or by", {
    data <- tibble(
        time = 1:10,
        ch1 = rnorm(10)
    )

    expect_error(
        shift_mnirs(data, nirs_channels = list("ch1"), time_channel = "time"),
        "to.*by"
    )
})

test_that("shift_mnirs shifts by constant correctly", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch1", "ch2"),
        time_channel = "time",
        by = 5,
        verbose = FALSE
    )

    expect_equal(result$ch1, data$ch1 + 5)
    expect_equal(result$ch2, data$ch2 + 5)
})

test_that("shift_mnirs preserves relative scaling within groups", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list(c("ch1", "ch2")),
        time_channel = "time",
        to = 0,
        width = 1,
        position = "min",
        verbose = FALSE
    )

    # ch1 min is 1, so shift both by -1
    expect_equal(result$ch1, 0:9)
    expect_equal(result$ch2, 10:19)
    # Relative difference preserved
    expect_equal(result$ch2[1] - result$ch1[1], 10)
})

test_that("shift_mnirs loses relative scaling across groups", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch1", "ch2"),
        time_channel = "time",
        to = 0,
        width = 1,
        position = "min",
        verbose = FALSE
    )

    # Each shifted independently to min = 0
    expect_equal(result$ch1, 0:9)
    expect_equal(result$ch2, 0:9)
    # Relative scaling lost
    expect_equal(result$ch2[1] - result$ch1[1], 0)
})

test_that("shift_mnirs handles position = 'max' correctly", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch1"),
        time_channel = "time",
        to = 100,
        width = 1,
        position = "max",
        verbose = FALSE
    )

    # Max is 10, shift to 100 means add 90
    expect_equal(result$ch1, 91:100)
})

test_that("shift_mnirs handles position = 'first' with span and width", {
    data <- tibble(
        time = seq(0, 9, by = 1),
        ch1 = c(rep(10, 3), 1:7)
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch1"),
        time_channel = "time",
        to = 0,
        position = "first",
        span = 2,
        verbose = FALSE
    )

    # Mean of first 3 values (time 0-2) is 10
    expect_equal(result$ch1, c(rep(0, 3), -9:-3))

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch1"),
        time_channel = "time",
        to = 0,
        position = "first",
        width = 3,
        verbose = FALSE
    )
    expect_equal(result$ch1, c(rep(0, 3), -9:-3))
})

test_that("shift_mnirs preserves unshifted channels", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20,
        ch3 = 21:30
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch1"),
        time_channel = "time",
        by = 5,
        width = 1,
        verbose = FALSE
    )

    expect_equal(result$ch1, 6:15)
    expect_equal(result$ch2, 11:20)
    expect_equal(result$ch3, 21:30)
})

test_that("shift_mnirs updates metadata correctly", {
    data <- tibble(time = 1:10, ch1 = 1:10, ch2 = 1:10)

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch1"),
        time_channel = "time",
        by = 5,
        width = 1,
        verbose = FALSE
    )

    expect_true("ch1" %in% attr(result, "nirs_channels"))
    expect_equal(attr(result, "time_channel"), "time")

    attr(data, "time_channel") <- NULL

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch2"),
        time_channel = "time",
        by = 5,
        width = 1,
        verbose = FALSE
    )

    expect_true("ch2" %in% attr(result, "nirs_channels"))
    expect_equal(attr(result, "time_channel"), "time")
})

test_that("shift_mnirs handles mixed channel groups", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10,
        ch2 = 11:20,
        ch3 = 31:40,
        ch4 = 41:50
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list(c("ch1", "ch2"), c("ch3", "ch4")),
        time_channel = "time",
        to = 0,
        width = 1,
        position = "min",
        verbose = FALSE
    )

    # Group 1: min is 1, both shifted by -1
    expect_equal(result$ch1, 0:9)
    expect_equal(result$ch2, 10:19)
    # Group 2: min is 31, both shifted by -31
    expect_equal(result$ch3, 0:9)
    expect_equal(result$ch4, 10:19)
    # Scaling preserved within, not between groups
    expect_equal(result$ch2[1] - result$ch1[1], 10)
    expect_equal(result$ch4[1] - result$ch3[1], 10)
})

test_that("shift_mnirs handles NA values correctly", {
    data <- tibble(
        time = 1:10,
        ch1 = c(NA, 2:10)
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list("ch1"),
        time_channel = "time",
        to = 0,
        width = 1,
        position = "min",
        verbose = FALSE
    )

    # Min of non-NA values is 2
    expect_true(is.na(result$ch1[1]))
    expect_equal(result$ch1[2:10], 0:8)
})

test_that("shift_mnirs informs when nirs_channels is not a list()", {
    data <- tibble(time = 1:5, ch1 = 1:5)

    ## fires: verbose=TRUE, no metadata, non-list channels
    expect_message(
        shift_mnirs(
            data,
            nirs_channels = "ch1",
            time_channel = "time",
            by = 1,
            verbose = TRUE
        ),
        "list\\(\\).*channel grouping"
    )

    ## silent: verbose=FALSE
    expect_no_message(
        shift_mnirs(
            data,
            nirs_channels = "ch1",
            time_channel = "time",
            by = 1,
            verbose = FALSE
        )
    )

    ## silent: nirs_channels already a list()
    expect_no_message(
        shift_mnirs(
            data,
            nirs_channels = list("ch1"),
            time_channel = "time",
            by = 1,
            verbose = TRUE
        )
    )

    ## silent: data has nirs_channels metadata
    attr(data, "nirs_channels") <- "ch1"
    expect_no_message(
        shift_mnirs(
            data,
            nirs_channels = "ch1",
            time_channel = "time",
            by = 1,
            verbose = TRUE
        )
    )
})

test_that("shift_mnirs prioritises to over by", {
    data <- tibble(
        time = 1:10,
        ch1 = 1:10
    )

    expect_message(
        result <- shift_mnirs(
            data,
            nirs_channels = list("ch1"),
            time_channel = "time",
            to = 0,
            by = 100,
            width = 1,
            position = "min",
            verbose = TRUE
        ),
        "to.*overrides.*by"
    )

    # Should use 'to', not 'by'
    expect_equal(result$ch1, 0:9)
})

test_that("shift_mnirs handles unevenly sampled data with span", {
    set.seed(13)
    # 0.5 Hz target (2 second intervals) with ±0.1s error
    base_times <- seq(0, 10, by = 2)
    time_error <- runif(length(base_times), -0.1, 0.1)

    data <- data.frame(
        time = round(base_times + time_error, 2),
        ch1 = c(10, 5, 20, 15, 25, 30)
    )

    # Test with 4-second span (should include ~2 samples per window)
    result <- shift_mnirs(
        data,
        list("ch1"),
        time_channel = "time",
        to = 0,
        span = 4,
        position = "first"
    )

    # First window should average first 3 values within 4s of start
    first_window_idx <- which(data$time <= data$time[1] + 4)
    expect_equal(which(data$time <= 4), first_window_idx)
    expected_mean <- mean(data$ch1[first_window_idx])
    expect_equal(result$ch1[1], data$ch1[1] - expected_mean)

    ## width = 3
    result <- shift_mnirs(
        data,
        list("ch1"),
        time_channel = "time",
        to = 0,
        width = 3,
        position = "first"
    )

    # First window should average first 3 values
    expect_equal(result$ch1[1], data$ch1[1] - mean(data$ch1[1:3]))

    ## unequal sampling and position = "minimum"
    result_min <- shift_mnirs(
        data,
        list("ch1"),
        time_channel = "time",
        to = 0,
        span = 4,
        position = "min"
    )

    x <- data$time
    min_mean <- Inf
    idx <- NA_real_
    for (i in seq_along(x)) {
        window_idx <- which(x >= x[i] & x <= x[i] + 4)
        if (length(window_idx) > 1) {
            # print(window_idx)
            window_mean <- mean(data$ch1[window_idx])
            # print(window_mean)
            min_mean <- min(min_mean, window_mean)
            if (min_mean == window_mean) idx <- window_idx
        }
    }
    expect_equal(idx, first_window_idx)
})

test_that("shift_mnirs handles multiple channel groups", {
    data <- data.frame(
        time = 1:2,
        ch1 = c(10, 20),
        ch2 = c(15, 25),
        ch3 = c(5, 35)
    )
    nirs_channels <- list(c("ch1", "ch2"), "ch3")
    result <- shift_mnirs(data, nirs_channels, time, to = 0, width = 1)

    ## check grouping together: min shuold come from ch1 and ch3
    expect_true(any(result$ch1 == 0, na.rm = TRUE))
    expect_false(any(result$ch2 == 0, na.rm = TRUE))
    expect_true(any(result$ch3 == 0, na.rm = TRUE))
    ## check both shifted together maintaining relative scaling
    expect_equal(result$ch1 - result$ch2, data$ch1 - data$ch2)
    ## check both shifted independently
    expect_false(isTRUE(all.equal(
        result$ch1 - result$ch3,
        data$ch1 - data$ch3
    )))
})

test_that("shift_mnirs preserves non-channel columns", {
    data <- data.frame(
        time = 1:3,
        ch1 = c(10, 20, 30),
        other = c("A", "B", "C")
    )
    result <- shift_mnirs(data, list("ch1"), "time", by = 5)

    expect_equal(result$time, c(1, 2, 3))
    expect_equal(result$other, c("A", "B", "C"))
    expect_equal(result$ch1, c(15, 25, 35))
})

test_that("shift_mnirs validates position argument", {
    data <- data.frame(ch1 = c(10, 20))
    expect_error(
        shift_mnirs(data, list("ch1"), position = "invalid"),
        "to.*by"
    )
})

test_that("shift_mnirs handles empty channel list", {
    data <- data.frame(time = 1:3, value = c(10, 20, 30))
    expect_error(shift_mnirs(data, list(), by = 5), "`nirs_channels` not detected")
})

## shift_mnirs() tidy eval integration ========================================
## helper to create test data with metadata
create_test_data <- function(
    time_max = 10,
    sample_rate = 10,
    add_metadata = TRUE
) {
    time <- seq(0, time_max, by = 1 / sample_rate)
    nrow <- length(time)
    data <- tibble(
        time = time,
        nirs1 = rnorm(nrow, 50, 5),
        nirs2 = rnorm(nrow, 60, 5),
        nirs3 = rnorm(nrow, 80, 5),
        event = c(1, rep(NA, nrow - 2), 2),
    )
    class(data) <- c("mnirs", class(data))

    if (add_metadata) {
        attr(data, "time_channel") <- "time"
        attr(data, "nirs_channels") <- c("nirs1", "nirs2")
        attr(data, "event_channel") <- "event"
        attr(data, "sample_rate") <- sample_rate
    }

    return(data)
}

test_that("shift_mnirs() works with quoted character strings", {
    data <- create_test_data()
    result <- shift_mnirs(
        data,
        nirs_channels = c("nirs1", "nirs2"),
        time_channel = "time",
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )

    expect_s3_class(result, "mnirs")
    expect_true(all(c("nirs1", "nirs2") %in% names(result)))
    expect_equal(
        result$nirs1[1],
        data$nirs1[1] - mean(c(data$nirs1[1:5], data$nirs2[1:5]))
    )
})

test_that("shift_mnirs() works with bare symbol column names", {
    data <- create_test_data()
    result <- shift_mnirs(
        data,
        nirs_channels = c(nirs1, nirs2),
        time_channel = time,
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_true(all(c("nirs1", "nirs2") %in% names(result)))
    expect_equal(
        result$nirs1[1],
        data$nirs1[1] - mean(c(data$nirs1[1:5], data$nirs2[1:5]))
    )
})

test_that("shift_mnirs() works with external character vector", {
    data <- create_test_data()
    channels <- c("nirs1", "nirs2")
    time_col <- "time"
    result <- shift_mnirs(
        data,
        nirs_channels = channels,
        time_channel = time_col,
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_true(all(c("nirs1", "nirs2") %in% names(result)))
    first_mean <- mean(c(data$nirs1[1:5], data$nirs2[1:5]))
    expect_equal(
        c(result$nirs1[1], result$nirs2[1]),
        c(
            data$nirs1[1] - first_mean,
            data$nirs2[1] - first_mean
        )
    )
})

test_that("shift_mnirs() works with external list object", {
    data <- create_test_data()
    channels <- list(c("nirs1", "nirs2"))
    result <- shift_mnirs(
        data,
        nirs_channels = channels,
        time_channel = "time",
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    first_mean <- mean(c(data$nirs1[1:5], data$nirs2[1:5]))
    expect_equal(
        c(result$nirs1[1], result$nirs2[1]),
        c(
            data$nirs1[1] - first_mean,
            data$nirs2[1] - first_mean
        )
    )
})

test_that("shift_mnirs() works with list() grouping and bare symbols", {
    data <- create_test_data()
    result <- shift_mnirs(
        data,
        nirs_channels = list(c(nirs1, nirs2)),
        time_channel = time,
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    first_mean <- mean(c(data$nirs1[1:5], data$nirs2[1:5]))
    expect_equal(
        c(result$nirs1[1], result$nirs2[1]),
        c(
            data$nirs1[1] - first_mean,
            data$nirs2[1] - first_mean
        )
    )
})

test_that("shift_mnirs() works with list() separate groups", {
    data <- create_test_data()
    result <- shift_mnirs(
        data,
        nirs_channels = list(nirs1, nirs2),
        time_channel = time,
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_equal(
        result$nirs1[1],
        data$nirs1[1] - mean(c(data$nirs1[1:5]))
    )
    expect_equal(
        result$nirs2[1],
        data$nirs2[1] - mean(c(data$nirs2[1:5]))
    )
})

test_that("shift_mnirs() works with tidyselect starts_with()", {
    data <- create_test_data()
    result <- shift_mnirs(
        data,
        nirs_channels = tidyselect::starts_with("nirs"),
        time_channel = time,
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_true(all(c("nirs1", "nirs2") %in% names(result)))
    first_mean <- mean(c(data$nirs1[1:5], data$nirs2[1:5], data$nirs3[1:5]))
    expect_equal(
        c(result$nirs1[1], result$nirs2[1], result$nirs3[1]),
        c(
            data$nirs1[1] - first_mean,
            data$nirs2[1] - first_mean,
            data$nirs3[1] - first_mean
        )
    )
})

test_that("shift_mnirs() works with tidyselect in list()", {
    data <- data.frame(
        time = 1:10,
        smo2_left = runif(10, 50, 70),
        smo2_right = runif(10, 50, 70),
        thb = runif(10, 12, 14)
    )
    data <- create_mnirs_data(
        data,
        list(
            nirs_channels = c("smo2_left", "smo2_right", "thb"),
            time_channel = "time"
        )
    )

    result <- shift_mnirs(
        data,
        nirs_channels = list(tidyselect::starts_with("smo2"), thb),
        time_channel = time,
        to = 0,
        width = 5,
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_true(all(c("smo2_left", "smo2_right", "thb") %in% names(result)))
})

test_that("shift_mnirs() uses metadata when channels NULL", {
    data <- create_test_data()
    result <- shift_mnirs(
        data,
        nirs_channels = NULL,
        time_channel = NULL,
        to = 0,
        width = 5,
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_true(all(c("nirs1", "nirs2") %in% names(result)))
})

test_that("shift_mnirs() errors with non-existent column name", {
    data <- create_test_data()
    expect_error(
        shift_mnirs(
            data,
            nirs_channels = nonexistent,
            time_channel = time,
            to = 0,
            width = 5,
            verbose = FALSE
        ),
        "not detected"
    )
})

test_that("shift_mnirs() mixed quoted and unquoted in list()", {
    data <- create_test_data()
    result <- shift_mnirs(
        data,
        nirs_channels = list("nirs1", nirs2),
        time_channel = "time",
        to = 0,
        width = 5,
        verbose = FALSE
    )
    expect_s3_class(result, "mnirs")
    expect_true(all(c("nirs1", "nirs2") %in% names(result)))
})

test_that("shift_mnirs() preserves grouping with external list", {
    data <- data.frame(
        time = 1:10,
        nirs1 = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
        nirs2 = c(15, 25, 35, 45, 55, 65, 75, 85, 95, 105),
        nirs3 = c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95)
    )
    data <- create_mnirs_data(
        data,
        list(
            nirs_channels = c("nirs1", "nirs2", "nirs3"),
            time_channel = "time"
        )
    )

    channels <- list(c("nirs1", "nirs2"), "nirs3")
    result <- shift_mnirs(
        data,
        channels,
        "time",
        to = 0,
        width = 5,
        position = "first",
        verbose = FALSE
    )

    ## nirs1 and nirs2 grouped:  shifted by same amount
    ## ch3 independent: shifted separately
    expect_equal(
        result$nirs1[1] - result$nirs2[1],
        data$nirs1[1] - data$nirs2[1]
    )
    first_mean <- mean(c(data$nirs1[1:5], data$nirs2[1:5]))
    expect_equal(
        c(result$nirs1[1], result$nirs2[1]),
        c(
            data$nirs1[1] - first_mean,
            data$nirs2[1] - first_mean
        )
    )
    expect_equal(
        result$nirs3[1],
        data$nirs3[1] - mean(c(data$nirs3[1:5]))
    )
})


## integration tests ================================================
test_that("shift_mnirs works on Moxy", {
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        verbose = FALSE
    ) |>
        dplyr::mutate(
            dplyr::across(
                dplyr::matches("smo2"),
                \(.x) {
                    replace_invalid(
                        .x,
                        invalid_values = c(0, 100),
                        method = "none"
                    )
                }
            )
        )

    data_shifted <- shift_mnirs(
        data,
        nirs_channels = c("smo2_left", "smo2_right"),
        time_channel = NULL,
        to = 0,
        by = NULL,
        span = 0,
        position = c("min", "max", "first"),
        verbose = FALSE
    )

    # plot(data) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))
    # plot(data_shifted) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))

    ## check grouping together: min value should come from smo2_right
    expect_false(any(data_shifted$smo2_left == 0, na.rm = TRUE))
    expect_true(any(data_shifted$smo2_right == 0, na.rm = TRUE))
    ## check both shifted together maintaining relative scaling
    expect_equal(
        data_shifted$smo2_left - data_shifted$smo2_right,
        data$smo2_left - data$smo2_right
    )
})

test_that("shift_mnirs(position = 'first') works on Moxy", {
    data <- read_mnirs(
        file_path = example_mnirs("moxy_ramp.xlsx"),
        nirs_channels = c(smo2_left = "SmO2 Live", smo2_right = "SmO2 Live(2)"),
        time_channel = c(time = "hh:mm:ss"),
        zero_time = TRUE,
        verbose = FALSE
    ) |>
        dplyr::mutate(
            dplyr::across(
                dplyr::matches("smo2"),
                \(.x) {
                    replace_invalid(
                        .x,
                        invalid_values = c(0, 100),
                        method = "none"
                    )
                }
            ),
            dplyr::across(
                dplyr::matches("smo2"),
                \(.x) replace_missing(.x, )
            )
        )

    data_shifted <- shift_mnirs(
        data,
        nirs_channels = c("smo2_left", "smo2_right"),
        time_channel = NULL,
        to = 0,
        by = NULL,
        span = 120,
        position = "first",
        verbose = FALSE
    )

    # plot(data) + ggplot2::ylim(0, 100)
    # plot(data_shifted) + ggplot2::geom_ribbon(
    #     data = ~ dplyr::filter(.x, time <= 120),
    #     aes(ymax = 10, ymin = -10), alpha = 0.1, linewidth = NA) +
    #     ggplot2::geom_segment(aes(x = -Inf, xend = 120, y = 0, yend = 0),
    #                           linetype = "dotted", colour = "black")

    first_mean <- data_shifted |>
        dplyr::filter(time <= 120) |>
        dplyr::summarise(
            mean = mean(c(smo2_left, smo2_right), na.rm = TRUE)
        ) |>
        dplyr::pull(mean)

    expect_equal(first_mean, 0)
})

test_that("shift_mnirs works on Train.Red", {
    data <- read_mnirs(
        file_path = example_mnirs("train.red_intervals.csv"),
        nirs_channels = c(
            smo2_left = "SmO2 unfiltered",
            smo2_right = "SmO2 unfiltered",
            o2hb_left = "O2HB unfiltered",
            o2hb_right = "O2HB unfiltered"
        ),
        time_channel = c(time = "Timestamp (seconds passed)"),
        verbose = FALSE
    )

    data_shifted <- shift_mnirs(
        data,
        nirs_channels = list(
            "smo2_left",
            "smo2_right",
            c("o2hb_left", "o2hb_right")
        ),
        time_channel = NULL,
        to = 0,
        by = NULL,
        span = 0, ## should default to 1 sample?
        position = "min",
        verbose = FALSE
    )

    data$time[which(diff(data$time) < 0)]

    # plot(data) + ggplot2::ylim(0, 100)
    # plot(data_shifted) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0))

    ## check grouping together: min value should come from each group
    expect_true(any(data_shifted$smo2_left == 0, na.rm = TRUE))
    expect_true(any(data_shifted$smo2_right == 0, na.rm = TRUE))
    expect_true(any(data_shifted$o2hb_left == 0, na.rm = TRUE))
    expect_false(any(data_shifted$o2hb_right == 0, na.rm = TRUE))
    ## check both shifted together maintaining relative scaling
    expect_equal(
        data_shifted$o2hb_left - data_shifted$o2hb_right,
        data$o2hb_left - data$o2hb_right
    )
    ## check both shifted independently
    expect_false(isTRUE(all.equal(
        data_shifted$smo2_left - data_shifted$smo2_right,
        data$smo2_left - data$smo2_right
    )))
})
