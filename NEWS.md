# mnirs 0.6.0

## Prepare for CRAN submission

* *NEWS.md* is truncated to the most recent relevant release updates. The full *NEWS.md* remains in the `dev` branch.

* Some development functions have been omitted from the package build in expectation of CRAN review. They are still present in `dev` branch and can be installed with `pak::pak("jemarnold/mnirs")`.

### Updated core functions

* `filter_mnirs()`: Method-specific arguments (e.g. `order` for `method = "butterworth"`, or `width` for `method = "moving_average"`) removed from the generic function call. Continue to be passed to the appropriate method via `...`. Documentation and info/warning/abort messages updated.

* `filter_ma()` better separates effects of `partial` and `na.rm` args:
    * `partial = FALSE` by default returns NA at edges where insufficient number of samples are available compared to the specified `width` or `span`.
    * `partial = TRUE` calculates mean values at edges, as long as one valid non-`NA` sample is available.
    * `na.rm = FALSE` by default behaves as expected with `mean(na.rm = FALSE)`, propagating any `NA`s in the local window to the calculated mean with a warning.
        * **NOTE** This differs from the behaviour of `na.rm = FALSE` in `filter_butter()`, which errors if there are any internal `NA`s present. This has not been changed.
    * `na.rm = TRUE` ignores `NA`s and calculates local means as long as one valid sample is present.

* `plot.mnirs()`
    * Add new arg: `points = TRUE` will plot points in addition to lines, as a useful quick shortcut.
    * Update `na.omit` now omits non-valid `c(NA, NaN, Inf, -Inf)` values from plotting, not just `NA`.

* `palette_mnirs()` can now accept either a single numeric value specifying the number of colours to return, or any number of (valid) character colour names.


# mnirs 0.5.2

### Shiny app

* Update online shiny app hosted at <https://jemarnold-mnirs-app.share.connect.posit.cloud/> with basic reading and pre-processing functionality.

# mnirs 0.5.1

### Updated core functions

* `read_mnirs()` now has better automated channel detection logic for *"Artinis" / "Oxysoft"* file formats. `nirs_channels` and `time_channel` can be left blank to automatically *"sample"*, *"time"* (from `sample_rate`), and nirs channel *"2"*.

* `replace_outliers()` now accepts numeric `outlier_cutoff` values for more precise outlier detection thresholds (previously was integers only). Documentation also updated.


# mnirs 0.5.0

### Updated core functions

* `read_mnirs()` expands `event_channel` to work with integer *"lap"* numbers, or character event label as previous.
    * Should now work on more .csv file formats; previously read errors may have occured where the file contained header rows above the data table, resulting in improper detection of columns.
    * `event_channel` can now be specified as an integer `lap` column, in addition to a character column as previous.
    * Other {mnirs} functions may expect `event_channel` to be either character or integer-ish.

* `extract_intervals()`:
    * Function arguments `start` and `end` are used to specify one or both of a start and end point to the target interval.
    * Specify `start`/`end` values with helper functions `by_time()`, `by_label()`, `by_lap()`, and `by_sample()`.
    * Numeric values are automatically coerced to "time" values; Explicit integer values (e.g. `2L`) are coerced to "lap"; Character strings are coerced to event "label".

### Package resources

* A rough draft *{mnirs}* hex icon and package cheatsheet have been added.


# mnirs 0.4.2 and prior

* Full previous changelog is available on the github `dev` branch: https://github.com/jemarnold/mnirs/blob/dev/NEWS.md