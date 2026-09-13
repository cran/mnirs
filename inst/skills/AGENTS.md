# `{mnirs}` Agent Reference

v0.8.0 | R (>= 4.1) | MIT | https://jemarnold.github.io/mnirs/

Muscle NIRS (mNIRS) processing: read `.xls(x)`/`.csv`/`.txt`/`.tsv`/`.ftn(2)` → resample → clean → filter → (transform) → extract intervals → analyse kinetics → plot.

---

## 1. `"mnirs"` class — tibble subclass + attributes

`attributes(data)`; `attr(data, "nirs_channels")`.

| Attribute | Type | Set by |
|---|---|---|
| `nirs_device` | chr(1) | `read_mnirs()` auto-detect |
| `nirs_channels` | chr | NIRS signal cols |
| `time_channel` | chr(1) | time col |
| `event_channel` | chr(1) | event/lap col |
| `sample_rate` | num(1) Hz | estimated if NULL |
| `start_timestamp` | POSIXct | absolute start |
| `interval_times` | `list(start, end)` | `extract_intervals()` |
| `interval_span` | num(2) | `extract_intervals()` |

`verbose = TRUE` on all functions. Omitted → `getOption("mnirs.verbose", TRUE)`. Explicit value overrides option.

---

## 2. Pipeline

```
read_mnirs()
├─ plot()                           # visualise at any step
└─ resample_mnirs()                 # regularise time grid, up/down-sample
   └─ replace_mnirs()               # invalid/outliers/NA
      └─ filter_mnirs()             # smooth
         ├─ shift_mnirs()           # shift level, keep amplitude
         ├─ rescale_mnirs()         # normalise range
         ├─ correct_blood_volume()  # optional: normalise THb changes
         └─ extract_intervals()     # list of interval dfs
            └─ analyse_kinetics()   # "mnirs_kinetics"
               └─ plot()            # fit + markers
```

Order matters: `resample → replace → filter`. Constraints §5.

---

## 3. Function Reference

**`data` input** (all df-level fns): single df → df; list of dfs → list; grouped df (`dplyr::group_by()`) → list per group. `extract_intervals()` flattens list input (§3.7). `plot.mnirs()` facets list input.

**Channel args**: `nirs_channels`, `time_channel`, `event_channel` default `NULL` → metadata. Accept `{tidyselect}`: `"smo2"`, `smo2`, `starts_with("time")`.

**Per-channel args**: named `list()` keyed by channel. `span = list(smo2 = 10)` → others default. `span = list(hhb = 10, 5)` → one unnamed = fallback. Unknown names warn, ignored. Also keyable by `group_channels` group (§3.5) or `group_intervals` group (§3.7, §3.8).

### 3.1 Read

```r
read_mnirs(
    file_path,
    nirs_channels = NULL,   # chr; rename c(new = "old")
    time_channel  = NULL,   # chr(1); rename c(time = "Timestamp")
    event_channel = NULL,   # chr(1); optional
    sample_rate   = NULL,   # Hz; estimated if NULL
    add_timestamp = FALSE,  # add POSIXct "timestamp" col
    zero_time     = FALSE,  # time[1] → 0
    keep_all      = FALSE,  # keep all cols
    verbose       = TRUE
)
example_mnirs(file = NULL)   # NULL = list all; partial match
create_mnirs_data(data, ...) # constructor; ... = named metadata
```

- `NULL` channels → auto-detect device, header row, channels, time col. No channels + known device → all cols returned.
- Duplicate col names → `_1` suffix. Unnamed cols → `col_n` (col number).
- Artinis Oxysoft: channels from Legend block, cleaned lowercase (`rx1_tx1_o2hb`); rename by number, cleaned name, or legend name. `(Sample number)` → `sample` + derived `time`; `(Event)` → `event`; trailing label col → `labels` (`event_channel = c(event = "labels")`).
- PIONIRS `.ftn(2)`: `Time`, `TagLabel`, `StO2*`; `Iteration`/`Tag` companion cols.
- Companion cols (`labels`, `Iteration`, `Tag`) only with `keep_all = TRUE`.
- POSIXct `time_channel` → numeric, rebased 0 regardless of `zero_time`.
- Sample-index `time_channel` → `sample_rate` mis-estimated 1 Hz; specify explicitly.
- Warns on irregular sampling (non-monotonic/repeated/unequal).
- Example files: `artinis_intervals.xlsx`, `moxy_intervals.csv`, `moxy_ramp.xlsx`, `pionirs_occlusion.ftn`, `portamon_oxcap.xlsx`, `train.red_intervals.csv`.

### 3.2 Resample

```r
resample_mnirs(data, time_channel = NULL, sample_rate = NULL,
    resample_rate = sample_rate,        # default = regularise only
    method = c("none", "linear", "locf"))
```

`"none"`: nearest match, new samples `NA`. `"linear"`/`"locf"`: `stats::approx()` numeric cols. Non-numeric: locf (up) / first-in-bin (down).

### 3.3 Clean

```r
replace_mnirs(data, nirs_channels = NULL, time_channel = NULL,
    invalid_values = NULL,  # exact values
    invalid_above  = NULL,  # x >= threshold
    invalid_below  = NULL,  # x <= threshold
    outlier_cutoff = NULL,  # Hampel MAD multiplier; NULL = skip; 3 ≈ 3σ, 2 ≈ Tukey
    width = NULL, span = NULL,  # window: samples XOR time; width wins
    method = c("linear", "median", "locf", "none"))
```

Order: invalid → outliers → NA.

Vector-level:
```r
replace_invalid(x, t, invalid_values, invalid_above, invalid_below, width, span, method = c("median", "none"))
replace_outliers(x, t, outlier_cutoff = 3, width, span, method = c("median", "none"))
replace_missing(x, t, width, span, method = c("linear", "median", "locf"))
```

### 3.4 Filter

```r
filter_mnirs(data, nirs_channels = NULL, time_channel = NULL,
    method = c("smooth_spline", "butterworth", "moving_average"),
    na.rm = FALSE, ...,
    spar = NULL,                          # smooth_spline; NULL = GCV
    order = 2L, W = NULL, fc = NULL, sample_rate = NULL,  # butterworth
    type = c("low", "high", "stop", "pass"), edges = c("rev", "rep1", "none"),
    width = NULL, span = NULL, partial = FALSE)           # moving_average
```

Method matching case/separator-insensitive. `smooth_spline` errors on NA/duplicated time. `butterworth` needs `{signal}`, errors on NA; `W` normalised or `fc` Hz.

Vector-level:
```r
filter_moving_average(x, t, width, span, partial = FALSE, na.rm = FALSE)  # alias filter_ma()
filter_butterworth(x, order = 2L, W, type = "low", edges = "rev", na.rm = FALSE)  # alias filter_butter()
```

### 3.5 Transform

```r
shift_mnirs(data, nirs_channels = NULL, time_channel = NULL,
    group_channels = c("ensemble", "distinct"),  # or list()
    to = NULL,    # target level; overrides `by`
    by = NULL,    # +/- shift
    width = NULL, span = NULL,
    position = c("min", "max", "first"))  # reference value
rescale_mnirs(data, nirs_channels = NULL,
    group_channels = c("ensemble", "distinct"),
    range)        # c(min, max)
```

`group_channels`: `"ensemble"` = one reference for all (relative scaling kept); `"distinct"` = per channel; `list(g1 = c("A", "B"), "C")` = custom groups, per-channel args keyable by group name.

### 3.6 Correct Blood Volume

```r
correct_blood_volume(data, oxy_channel = NULL, deoxy_channel = NULL,
    total_channel = NULL)  # derived from pair if NULL
```

Vectors paired by position for multiple pairs. Ryan 2012; Beever & Tripp 2020. Negative values handled via `shift_mnirs()`. Corrected `total` → 0.

### 3.7 Extract Intervals

```r
extract_intervals(data, nirs_channels = NULL, time_channel = NULL,
    event_channel = NULL, sample_rate = NULL,
    group_intervals = c("distinct", "ensemble"),  # or list()
    group_channels = NULL,      # channels ensemble-averaged per group
    start = NULL, end = NULL,   # by_*(); end NULL = span around start
    span = list(c(-60, 60)),    # c(before_start, after_end)
    zero_time = FALSE)          # rebase per interval
## event_groups = deprecated()
## → flat named list of "mnirs" dfs
by_time(...); by_label(..., ignore_case = FALSE, fixed = FALSE)  # regex default
by_lap(...); by_sample(...)   # by_label/by_lap need event_channel
```

- `start`/`end`: `list()` of mixed `by_*()` allowed; times concatenated in order.
- `span`: negative = earlier, positive = later. Scalar recycled by sign: `60` → `c(0, 60)`; `-60` → `c(-60, 0)`. `list()` per interval group, recycled.
- Names: `interval_<df>.<interval>`; ensemble/custom → `<name>_<df>`.
- `group_intervals`: `"distinct"` = one df per interval; `"ensemble"` = one averaged df (needs regular time grid); `list(g1 = c(1, 2), c(3, 4))` = one df per group.

### 3.8 Analyse Kinetics

```r
analyse_kinetics(data, nirs_channels = NULL, time_channel = NULL,
    method = c("response_time", "peak_slope", "monoexponential",
               "exponential_drift", "biexponential", "sigmoidal", "sigmoidal_drift"),
    start_time = NULL,   # fit t=0; NULL → interval_times metadata → t[1] → 0
    direction = c("auto", "positive", "negative"),
    end_window = Inf,    # fit end = first extreme with none greater within window; Inf = global
    group_intervals = "ensemble",  # or list() of row indices
    zero_time = FALSE,   # rebase per interval/group; shifts interval_times
    ...,                 # control = list()/nls.control() for nls methods; global only
    response_fraction = 0.5,                          # response_time; vectorised
    width = NULL, span = NULL, align = c("centre", "left", "right"),
    partial = FALSE, na.rm = FALSE,                   # peak_slope
    use_TD = TRUE,                                    # exp methods
    shape = c("symmetric", "gompertz", "gompertz_left"),  # sigmoidal methods
    drift_fraction = NULL,                            # drift methods; NULL → 0.95, range (0.5, 1)
    fix = NULL)
## analyze_kinetics() alias
```

- Method aliases case/separator-insensitive + shorthand: `hrt`, `slope`, `mrt`, `tau`, `gompertz`, `xmid`.
- `method` global only. All other args per-channel (`list(smo2 = 10)`) and per-interval (`list(interval_1 = list(smo2 = 10))`) except `control`.
- `direction` constrains fitted-amplitude sign for nls methods; `NA` if unsatisfiable.
- `use_TD = TRUE`: 4-param → 3-param fallback (monoexp), 6 → 5 (biexp).
- `"exponential_drift"`: monoexp + hinge-linear drift from onset `TD - tau·log(1 - drift_fraction)` (= `TD + 3·tau` default; `expdrift_onset()`; held constant). `texc` = `max(onset, TD + tau·log(|B-A| / (|slope_B|·tau)))`. Falls back to monoexp on failure or `|slope_B|·(t_end - onset) < 2·rmse`.
- `"biexponential"`: stage 1 = monoexp on `end_window` (`Inf` → 30 time units past first extreme); stage 2 = full biexp, `A`/`tau`/`TD` box-bounded near stage 1 (undocumented `tau_flex = 1/3`, `TD_flex = 2`, `A_flex = NULL` → 2·sd resid), `B`/`B2`/`tau2` free. Falls back exp_drift → monoexp on failure, monotonic `texc`, `tau2 >= 2·span`, or `|B2 - B| < 2·rmse`. Coef cols = union of chain.
- `"sigmoidal"`: `shape`: `"symmetric"` = `SSlogistic()`; `"gompertz"` early inflection; `"gompertz_left"` late.
- `"sigmoidal_drift"`: sigmoid + hinge-linear drift from onset (analytic inverse per shape; `sigdrift_onset()`). `texc` = `max(onset, t where |S'(t)| = |slope_B|)` (`sigdrift_texc()`; `uniroot()` for Gompertz). Falls back to sigmoidal (same shape) on failure or `|slope_B|·(t_end - texc) < 2·rmse`.
- Fallbacks: warning + `model` column. Undocumented `model_fallback = FALSE` keeps raw fit.
- `control`: merged over `maxiter = 500, warnOnly = TRUE` on every `nls()` call incl. refits/fallbacks. Unknown names abort.

**`group_intervals`** (rows; no `"distinct"`): `"ensemble"` = all rows per df; `list(trial1 = 1:12, trial2 = 13:24)` = each group separate interval; unnamed → `interval_<n>`; multi-df → `<group>_<df>`. Ungrouped rows dropped (message); overlaps warn. Main use = recursive fit on `"mnirs_kinetics"` coefs:
```r
analyse_kinetics(result, nirs_channels = slope, time_channel = peak_slope_time,
    method = "monoexp", group_intervals = list(trial1 = 1:12, trial2 = 13:24))
```

**`fix`**: `list(A = 0)` global; `list(smo2 = list(A = 0))` per channel; `list(interval_1 = list(smo2 = list(A = 0)))` per interval×channel. Finite scalars; names = model params; cannot fix all. Fixed params excluded from `n_params`. Fixed `TD` needs `use_TD = TRUE` and disables reduced-param fallback; prefer `use_TD = FALSE` over `fix = list(TD = 0)`.

**Return `"mnirs_kinetics"`** (list; prints coef table; one row per channel per interval):

| Element | Content |
|---|---|
| `method` | canonical method |
| `model` | nested list of `lm`/`nls` per interval per channel |
| `coefficients` | df |
| `data` | input dfs + `<channel>_fitted` cols |
| `interval_times` | df: `interval`, `start_times`, `end_times` (if present) |
| `diagnostics` | `n_obs`, `n_params` (free, excl. `fix`), `r2`, `adj_r2`, `rmse`, `snr`, `cv_rmse`, `aic`, `aicc`, `bic` |
| `channel_args` | resolved args |
| `warnings` | df `type` = `"warning"`/`"error"`; captured regardless of `verbose` |
| `call` | matched call |

**Coefficients** (all prefixed `interval`, `nirs_channels`, `start_time`; times elapsed from `start_time`; `*_fitted` = predicted value at that time):

| Method | Columns |
|---|---|
| `response_time` | `response_fraction` (row per value), `A` baseline mean, `B` extreme, `response_time`, `response_value`, `fitted` (`A + (B-A)·f`), `idx` |
| `peak_slope` | `slope`, `intercept`, `fitted`, `peak_slope_time`, `idx` (at `align`) |
| `monoexponential` | `A`, `B`, `tau`, `k` (`1/tau`), `TD`, `MRT` (`TD+tau`), `HRT` (`TD+tau·ln2`), `MRT_fitted`, `HRT_fitted` |
| `exponential_drift` | monoexp + `texc`, `slope_B`, `drift_fraction`, `texc_fitted`, `model` |
| `biexponential` | `A`, `B`, `tau` (fast), `MRT`, `texc` (`NA` if monotonic), `B2`, `tau2` (slow), `TD`, `MRT_fitted`, `texc_fitted`, `model` + fallback cols (`NA` unless used) |
| `sigmoidal` | `A`, `B`, `xmid` (inflection; midpoint only if symmetric), `slope` (at `xmid`), `xmid_fitted` |
| `sigmoidal_drift` | sigmoidal + `texc`, `slope_B`, `drift_fraction`, `texc_fitted`, `model` |

`aic`/`bic` comparable only within matching `n_obs` + `n_params`.

**Vector-level / model fns:**
```r
response_time(x, t = seq_along(x), start_time = 0, response_fraction = 0.5, direction)
## → A, B, response_time, response_value, fitted, baseline_idx, response_idx, extreme_idx
peak_slope(x, t = seq_along(x), width, span, align, direction, partial = FALSE, na.rm = FALSE)
## → slope, intercept, y, t, idx, fitted, window_idx, model

monoexponential(t, A, B, tau, TD = NULL)   # A + (B-A)(1 - exp(-pmax(t-TD,0)/tau))
exponential_drift(t, A, B, tau, slope_B, drift_fraction, TD = NULL)
biexponential(t, A, B, tau, B2, tau2, TD = NULL)
## A + (B-A)(1-exp(-ts/tau)) + (B2-B)(1-exp(-ts/tau2)), ts = pmax(t-TD, 0)
## texc = TD + log(r)/(1/tau - 1/tau2), r = -(B-A)tau2/((B2-B)tau); NA if r <= 1
logistic(t, A, B, xmid, slope, asym = NULL)  # 4-param; asym → 5-param Richards
gompertz(t, A, B, xmid, slope); gompertz_left(t, A, B, xmid, slope)
sigmoidal_drift(t, A, B, xmid, slope, slope_B, drift_fraction, shape)
## S(t) + slope_B·pmax(t - onset, 0); onset = xmid + u/k
## symmetric: u = log(f/(1-f)), k = 4·slope/(B-A)
## gompertz: u = -log(-log f); gompertz_left: u = log(-log(1-f)); k = slope·e/(B-A)

## selfStart for nls(); port algorithm where bounded
nls(x ~ SSmonoexponential(t, A, B, tau, TD), data = df)         # drop TD → 3-param
nls(x ~ SSexponential_drift(t, A, B, tau, slope_B, drift_fraction, TD), data = df)
nls(x ~ SSbiexponential(t, A, B, tau, B2, tau2, TD), data = df,
    algorithm = "port", lower = c(-Inf, -Inf, 0, -Inf, 0, 0))
nls(x ~ SSlogistic(t, A, B, xmid, slope, asym), data = df)       # 5-param fragile
nls(x ~ SSgompertz(t, A, B, xmid, slope), data = df)
nls(x ~ SSgompertz_left(t, A, B, xmid, slope), data = df)
nls(x ~ SSsigmoidal_drift(t, A, B, xmid, slope, slope_B, drift_fraction = 0.95,
    shape = "gompertz"), data = df, algorithm = "port", control = nls.control(warnOnly = TRUE))
```

### 3.9 Plot

Needs `{ggplot2}`; `{scales}` for axis formatting.

```r
plot.mnirs(x, points = FALSE, time_labels = FALSE, na.omit = FALSE, ...)
## x = df or list (faceted); ... = facet_wrap args, n.breaks, breaks
plot.mnirs_kinetics(x, fitted = TRUE, markers = TRUE, labels = TRUE, ...)
## fitted = dashed curve (not response_time); markers = onset line + coef points;
## labels = panel annotation; ... = label_size, plot.mnirs() args
theme_mnirs(base_size = 14, base_family = "sans", border = c("partial", "full"),
    ink = "black", paper = "white", accent = "#0080ff", ...)
palette_mnirs(...)   # () = all 12; (4); ("red", "blue")
scale_colour_mnirs(..., aesthetics = "colour")  # alias scale_color_mnirs()
scale_fill_mnirs(..., aesthetics = "fill")
breaks_timespan(unit = c("secs", "mins", "hours", "days", "weeks"), n = 5)
format_hmmss(x)      # secs → "mm:ss" / "h:mm:ss"
```

---

## 4. Dependencies

Imports: `cli`, `data.table`, `lifecycle`, `readxl`, `rlang`, `stats`, `tibble`, `tidyselect`, `utils`.

Suggests: `signal` (butterworth), `ggplot2` + `scales` (plot/theme/scales), `dplyr` (grouped df input), `knitr` + `quarto` (vignettes), `zoo` + `testthat` (tests).

---

## 5. Constraints

| Constraint | Detail |
|---|---|
| `read_mnirs()` irregular-sample warning | fires before downstream `resample_mnirs()`; verify output |
| Pipeline order | `resample → replace → filter`; other orders change results |
| `extract_intervals(group_intervals = "ensemble")` | needs regular time grid; warns if irregular |
| nls convergence | weak fits return with warnings; check `warnings` |
| `direction` bound | `NA` coefs if unsatisfiable; check `biexponential` |
| `biexponential` identifiability | fast phase from stage 1 on `end_window`; too long → slow stage 1 → fallback; check `warnings` + `model` |

---

## 6. Source Map

| File | Contents |
|---|---|
| `R/read_mnirs.R` | `read_mnirs()`, `example_mnirs()`, `create_mnirs_data()` |
| `R/read_mnirs_helpers.R` | device/header/channel detection |
| `R/resample_mnirs.R` | `resample_mnirs()` |
| `R/replace_mnirs.R` | `replace_mnirs()`, `replace_invalid/outliers/missing()` |
| `R/filter_mnirs.R` | `filter_mnirs()`, `filter_moving_average()`, `filter_butterworth()` + aliases |
| `R/rolling_helpers.R` | `compute_window_bounds()`, `compute_outliers()` |
| `R/shift_mnirs.R`, `R/rescale_mnirs.R` | `shift_mnirs()`, `rescale_mnirs()` |
| `R/correct_blood_volume.R` | `correct_blood_volume()` |
| `R/extract_intervals.R` | `extract_intervals()` |
| `R/extract_interval_helpers.R` | `by_*()`, boundary resolution, `validate_interval_channels/groups()` |
| `R/analyse_kinetics.R` | `analyse_kinetics()`/`analyze_kinetics()`, S3 methods per `method` |
| `R/aanalyse_kinetics_helpers.R` | (`aa` = load order) `method_aliases`, `kinetics_fallbacks`, `analyse_kinetics_intervals()` → `analyse_<method>()` → `analyse_kinetics_channels()` → `fit_<method>()`; `detect_direction()`, `enforce_direction()`, `compute_diagnostics()`, `kinetics_warnings_df()`, `validate_kinetics_args()` |
| `R/analyse_response_time.R` | `response_time()` |
| `R/analyse_peak_slope.R` | `peak_slope()`, `rolling_slope()` |
| `R/analyse_monoexponential.R` | `monoexponential()`, `SSmonoexponential()` |
| `R/analyse_exponential_drift.R` | `exponential_drift()`, `expdrift_onset()`, `expdrift_start()`, `SSexponential_drift()` |
| `R/analyse_biexponential.R` | `biexponential()`, `biexp_texc()`, `biexp_init()`, `SSbiexponential()` |
| `R/analyse_sigmoidal.R` | `logistic()`, `gompertz()`, `gompertz_left()`, `sigmoid_core()`, `SS*()` |
| `R/analyse_sigmoidal_drift.R` | `sigmoidal_drift()`, `sigdrift_onset()`, `sigdrift_texc()`, `sigdrift_start()`, `SSsigmoidal_drift()` |
| `R/plot.mnirs.R` | `plot.mnirs()`, `plot.mnirs_kinetics()`, `theme_mnirs()`, `palette_mnirs()`, scales, `breaks_timespan()`, `format_hmmss()` |
| `R/mnirs_methods.R` | `print.mnirs()`, `print.mnirs_kinetics()` |
| `R/channel_args.R` | `resolve_channel_args()`, `validate_group_channels()` |
| `R/as_data_list.R` | `as_data_list()`, `map_mnirs_intervals()` — list/grouped dispatch |
| `R/validate_mnirs.R` | `validate_numeric/mnirs_data/nirs_channels/time_channel/event_channel/sample_rate/width_span/x_t/start_time/fix/findInt()` |
| `R/signif_trailing.R` | `signif_trailing()`, `seq_range()` |
| `R/data.R`, `R/mnirs-package.R` | example file docs, package roxygen |

Messages: `cli_abort()`/`cli_warn()`/`cli_inform()`. Roxygen2 markdown. pkgdown: `_pkgdown.yml`.
