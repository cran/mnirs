test_that("rescale_mnirs rescales single channel correctly", {
    data <- tibble(A = c(0, 50, 100), B = c(10, 20, 30))

    result <- rescale_mnirs(
        data,
        nirs_channels = list("A"),
        range = c(0, 1)
    )

    expect_equal(result$A, c(0, 0.5, 1))
    expect_equal(result$B, data$B) # unchanged
})

test_that("rescale_mnirs preserves relative scaling across grouped channels", {
    data <- tibble(A = c(0, 50, 100), B = c(25, 25, 50))

    result <- rescale_mnirs(
        data,
        nirs_channels = list(c("A", "B")),
        range = c(0, 1)
    )

    expect_equal(result$A, c(0, 0.5, 1))
    expect_equal(result$B, c(0.25, 0.25, 0.5))
})

test_that("rescale_mnirs handles multiple separate groups", {
    data <- tibble(A = c(10, 100), B = c(0, 50), C = c(10, 200))

    result <- rescale_mnirs(
        data,
        nirs_channels = list(c("A", "B"), "C"),
        range = c(0, 10)
    )

    expect_equal(result$A, c(1, 10))
    expect_equal(result$B, c(0, 5))
    # C scaled independently
    expect_equal(result$C, c(0, 10))
})

test_that("rescale_mnirs handles negative ranges", {
    data <- tibble(A = c(0, 50, 100))

    result <- rescale_mnirs(
        data,
        nirs_channels = list("A"),
        range = c(-1, 1)
    )

    expect_equal(result$A, c(-1, 0, 1))
})

test_that("rescale_mnirs handles NA values", {
    data <- tibble(A = c(0, NA, 100))

    result <- rescale_mnirs(
        data,
        nirs_channels = list("A"),
        range = c(0, 1)
    )

    expect_equal(result$A, c(0, NA, 1))
})

test_that("rescale_mnirs errors for nirs_channels", {
    data <- tibble(A = c(0, 50, 100))

    expect_error(
        rescale_mnirs(data, nirs_channels = NULL, range = c(0, 1)),
        "nirs_channels.*not detected"
    )
    expect_error(
        rescale_mnirs(data, range = c(0, 1)),
        "nirs_channels.*not detected"
    )
    expect_error(
        rescale_mnirs(data, nirs_channels = "doesn't exist", range = c(0, 1)),
        "nirs_channels.*match exactly"
    )
})

test_that("rescale_mnirs errors with invalid range", {
    data <- tibble(A = c(0, 50, 100))

    expect_error(
        rescale_mnirs(data, nirs_channels = list("A"), range = c(0, 1, 2)),
        "range.*numeric"
    )
    expect_error(
        rescale_mnirs(data, nirs_channels = list("A"), range = "invalid"),
        "range.*numeric"
    )
})

test_that("rescale_mnirs returns unmodified column when values are constant", {
    data <- tibble(A = c(50, 50, 50), B = c(0, 100, 200))

    ## ungrouped should not change
    result <- rescale_mnirs(
        data,
        nirs_channels = list("A", "B"),
        range = c(0, 100)
    )

    expect_equal(result$A, c(50, 50, 50))
    expect_equal(result$B, c(0, 50, 100))

    ## grouped should change
    result <- rescale_mnirs(
        data,
        nirs_channels = list(c("A", "B")),
        range = c(0, 100)
    )

    # When grouped, if one channel is constant and another varies,
    # the constant channel should remain unchanged
    expect_equal(result$A, data$A / 2)
    expect_equal(result$B, data$B / 2)
})

test_that("rescale_mnirs informs when nirs_channels is not a list()", {
    data <- tibble(A = c(0, 50, 100))

    ## fires: verbose=TRUE, no metadata, non-list channels
    expect_message(
        rescale_mnirs(
            data,
            nirs_channels = "A",
            range = c(0, 1),
            verbose = TRUE
        ),
        "list\\(\\).*channel grouping"
    )

    ## silent: verbose=FALSE
    expect_no_message(
        rescale_mnirs(
            data,
            nirs_channels = "A",
            range = c(0, 1),
            verbose = FALSE
        )
    )

    ## silent: nirs_channels already a list()
    expect_no_message(
        rescale_mnirs(
            data,
            nirs_channels = list("A"),
            range = c(0, 1),
            verbose = TRUE
        )
    )

    ## silent: data has nirs_channels metadata
    attr(data, "nirs_channels") <- "A"
    expect_no_message(
        rescale_mnirs(
            data,
            nirs_channels = "A",
            range = c(0, 1),
            verbose = TRUE
        )
    )
})

test_that("rescale_mnirs updates metadata correctly", {
    data <- tibble(A = c(50, 50, 50), B = c(0, 100, 200))

    result <- rescale_mnirs(
        data,
        nirs_channels = list("A", "B"),
        range = c(0, 100),
        verbose = FALSE
    )

    expect_true(all(c("A", "B") %in% attr(result, "nirs_channels")))

    attr(data, "nirs_channels") <- "A"
    result <- rescale_mnirs(
        data,
        nirs_channels = list("A", "B"),
        range = c(0, 100),
        verbose = FALSE
    )

    expect_true(all(c("A", "B") %in% attr(result, "nirs_channels")))
})

test_that("rescale_mnirs works on Moxy", {
    file_path <- example_mnirs("moxy_ramp.xlsx")

    df <- read_mnirs(
        file_path = file_path,
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

    result <- rescale_mnirs(
        df,
        nirs_channels = list(c("smo2_left", "smo2_right")),
        range = c(0, 100)
    )

    # plot(df) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))
    # plot(result) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))

    ## check grouping together: min & max value should come from smo2_right
    expect_false(any(result$smo2_left %in% c(0, 100), na.rm = TRUE))
    expect_true(any(result$smo2_right %in% c(0, 100), na.rm = TRUE))

    ## check grouped apart
    result <- rescale_mnirs(
        df,
        nirs_channels = list("smo2_left", "smo2_right"),
        range = c(0, 100)
    )

    # plot(df) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))
    # plot(result) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))

    ## check grouping together: min value should come from smo2_right
    expect_true(any(result$smo2_left %in% c(0, 100), na.rm = TRUE))
    expect_true(any(result$smo2_right %in% c(0, 100), na.rm = TRUE))
})

test_that("rescale_mnirs works on Train.Red", {
    file_path <- example_mnirs("train.red_intervals.csv")

    df <- read_mnirs(
        file_path = file_path,
        nirs_channels = c(
            smo2_left = "SmO2 unfiltered",
            smo2_right = "SmO2 unfiltered",
            o2hb_left = "O2HB unfiltered",
            o2hb_right = "O2HB unfiltered"
        ),
        time_channel = c(time = "Timestamp (seconds passed)"),
        verbose = FALSE
    )

    result <- rescale_mnirs(
        df,
        nirs_channels = list(
            "smo2_left",
            "smo2_right",
            c("o2hb_left", "o2hb_right")
        ),
        range = c(0, 100)
    )

    # plot(df) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))
    # plot(result) + ggplot2::ylim(0, 100) + geom_hline(yintercept = c(0, 100))

    ## check grouping together: min value should come from each group
    expect_true(any(result$smo2_left %in% c(0, 100), na.rm = TRUE))
    expect_true(any(result$smo2_right %in% c(0, 100), na.rm = TRUE))
    expect_true(any(result$o2hb_left == 0, na.rm = TRUE))
    expect_false(any(result$o2hb_right == 0, na.rm = TRUE))
    expect_false(any(result$o2hb_left == 100, na.rm = TRUE))
    expect_true(any(result$o2hb_right == 100, na.rm = TRUE))
})
