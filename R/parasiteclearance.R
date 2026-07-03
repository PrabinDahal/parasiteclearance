#' Estimate parasite clearance half-life for one patient profile
#'
#' Estimates the parasite clearance rate constant and slope half-life from
#' serial parasitaemia-time data using a Flegg et al.-style algorithm.
#'
#' @param data A data frame containing time and parasitaemia columns.
#' @param time_col Name of the time column, in hours.
#' @param parasite_col Name of the parasitaemia column, parasites per microlitre.
#' @param detection_limit Detection limit in parasites per microlitre.
#' @param initial_minimum Minimum acceptable initial parasitaemia. Default is 1000.
#' @param final_maximum Maximum acceptable final parasitaemia for estimation. Default is 1000.
#' @param min_first24 Minimum number of measurements required in first 24 hours to estimate lag phase.
#' @param max_gap_first24 Maximum allowed gap between measurements in first 24 hours to estimate lag phase.
#' @param clean Logical. If TRUE, perform automated data cleaning.
#' @param return_data Logical. If TRUE, return cleaned data and fitted values.
#'
#' @return A list with clearance estimates, selected model, lag phase, half-life, and status.
#'
#' @export
estimate_clearance <- function(data,
                               time_col = "time",
                               parasite_col = "parasitaemia",
                               detection_limit = 16,
                               initial_minimum = 1000,
                               final_maximum = 1000,
                               min_first24 = 3,
                               max_gap_first24 = 14,
                               clean = TRUE,
                               return_data = TRUE) {
  profile <- prepare_profile(
    data = data,
    time_col = time_col,
    parasite_col = parasite_col
  )

  if (clean) {
    cleaned <- clean_profile(
      profile,
      detection_limit = detection_limit
    )
  } else {
    cleaned <- profile
    cleaned$flag <- "kept"
    cleaned$used <- TRUE
  }

  analysis_data <- cleaned[cleaned$used, , drop = FALSE]
  analysis_data <- analysis_data[order(analysis_data$time), , drop = FALSE]

  checks <- check_estimable(
    analysis_data,
    initial_minimum = initial_minimum,
    final_maximum = final_maximum
  )

  if (!checks$estimable) {
    return(make_result(
      status = "not_estimable",
      reason = checks$reason,
      cleaned_data = if (return_data) cleaned else NULL
    ))
  }

  lag_estimable <- check_lag_estimable(
    analysis_data,
    min_first24 = min_first24,
    max_gap_first24 = max_gap_first24
  )

  if (!lag_estimable) {
    fit <- fit_linear_or_tobit(analysis_data)
    k <- -coef(fit$model)["time"]

    return(make_result(
      status = "estimated",
      reason = "lag_phase_not_estimable_linear_model_used",
      model_type = fit$model_type,
      k = as.numeric(k),
      slope_half_life = log(2) / as.numeric(k),
      tlag = 0,
      cleaned_data = if (return_data) cleaned else NULL,
      fit = fit
    ))
  }

  fit_candidates <- fit_candidate_models(analysis_data)
  best_fit <- select_best_model(fit_candidates)

  lag_result <- estimate_lag_and_k(
    data = analysis_data,
    best_fit = best_fit
  )

  make_result(
    status = "estimated",
    reason = "ok",
    model_type = best_fit$model_type,
    k = lag_result$k,
    slope_half_life = log(2) / lag_result$k,
    tlag = lag_result$tlag,
    cleaned_data = if (return_data) cleaned else NULL,
    fit = best_fit,
    linear_data = lag_result$linear_data
  )
}


#' Clean one parasite-time profile
#'
#' @param profile Data frame with columns `time` and `parasitaemia`.
#' @param detection_limit Detection limit in parasites per microlitre.
#'
#' @return Data frame with cleaning flags and `used` column.
#' @export
clean_profile <- function(profile, detection_limit = 16) {
  x <- profile[order(profile$time), , drop = FALSE]
  x$flag <- "kept"
  x$used <- TRUE

  x <- remove_recurrence_data(x)
  x <- remove_trailing_zeros(x)
  x <- remove_tails(x)
  x <- replace_first_sustained_zero(x, detection_limit)
  x <- remove_extreme_values(x)
  x <- remove_outliers(x)

  x
}


remove_recurrence_data <- function(x) {
  if (nrow(x) == 0) return(x)

  remove_after <- rep(FALSE, nrow(x))

  # Remove measurements after 7 days, i.e. after 168 hours.
  remove_after <- remove_after | x$time > 168

  # Remove data after two consecutive measurements are more than 24 hours apart.
  gaps <- diff(x$time)
  large_gap_index <- which(gaps > 24)

  if (length(large_gap_index) > 0) {
    first_cut <- large_gap_index[1] + 1
    remove_after[first_cut:nrow(x)] <- TRUE
  }

  x$flag[remove_after] <- "recurrence_or_late_data"
  x$used[remove_after] <- FALSE
  x
}


remove_trailing_zeros <- function(x) {
  used <- which(x$used)
  if (length(used) == 0) return(x)

  y <- x[used, , drop = FALSE]
  zero_idx <- which(y$parasitaemia == 0)

  if (length(zero_idx) == 0) return(x)

  last_positive <- max(which(y$parasitaemia > 0), na.rm = TRUE)
  trailing_zero_idx <- zero_idx[zero_idx > last_positive]

  if (length(trailing_zero_idx) > 1) {
    remove_local <- trailing_zero_idx[-1]
    remove_global <- used[remove_local]
    x$flag[remove_global] <- "trailing_zero"
    x$used[remove_global] <- FALSE
  }

  x
}


remove_tails <- function(x) {
  used <- which(x$used)
  if (length(used) < 3) return(x)

  y <- x[used, , drop = FALSE]

  low <- y$parasitaemia > 0 & y$parasitaemia < 100

  if (sum(low, na.rm = TRUE) < 2) return(x)

  low_idx <- which(low)

  # Tail: repeated low parasitaemias near detection limit at terminal end.
  # Conservative implementation: identify first low value after which all
  # remaining positive values are < 100 or zero.
  for (idx in low_idx) {
    remaining <- y$parasitaemia[idx:nrow(y)]
    positive_remaining <- remaining[remaining > 0]

    if (length(positive_remaining) >= 2 && all(positive_remaining < 100)) {
      remove_local <- idx:nrow(y)
      remove_global <- used[remove_local]

      # Keep the first point of tail as boundary? Flegg removes repeated
      # parasitaemias below 100 and data between them. For estimation,
      # exclude the terminal tail.
      x$flag[remove_global] <- "tail"
      x$used[remove_global] <- FALSE
      break
    }
  }

  x
}


replace_first_sustained_zero <- function(x, detection_limit) {
  used <- which(x$used)
  if (length(used) == 0) return(x)

  y <- x[used, , drop = FALSE]

  zero_idx <- which(y$parasitaemia == 0)
  if (length(zero_idx) == 0) return(x)

  for (zi in zero_idx) {
    if (zi > 1 && y$parasitaemia[zi - 1] > 0) {
      global_i <- used[zi]
      x$parasitaemia[global_i] <- detection_limit
      x$flag[global_i] <- "zero_replaced_detection_limit"
      break
    }
  }

  x
}


remove_extreme_values <- function(x) {
  bad <- x$used & (
    is.na(x$time) |
      is.na(x$parasitaemia) |
      x$time < 0 |
      x$parasitaemia < 0 |
      x$parasitaemia > 3e6
  )

  x$flag[bad] <- "extreme_or_invalid"
  x$used[bad] <- FALSE
  x
}


remove_outliers <- function(x) {
  used <- which(x$used)
  if (length(used) < 4) return(x)

  y <- x[used, , drop = FALSE]
  logp <- log(y$parasitaemia)
  times <- y$time

  slopes <- diff(logp) / diff(times)
  avg_slope <- (tail(logp, 1) - logp[1]) / (tail(times, 1) - times[1])

  if (!is.finite(avg_slope) || avg_slope == 0) return(x)

  norm_slope <- slopes / avg_slope

  remove_local_points <- integer(0)

  for (i in seq_len(length(norm_slope) - 1)) {
    t_mid <- y$time[i + 1]

    ns_i <- norm_slope[i]
    ns_next <- norm_slope[i + 1]

    if (!is.finite(ns_i) || !is.finite(ns_next)) next

    if (t_mid <= 12) {
      if (ns_i < -20 && ns_next > 10) {
        remove_local_points <- c(remove_local_points, i + 1)
      }
    } else {
      if ((ns_i < -7.5 && ns_next > 10) ||
          (ns_i < -40 && ns_next > 3.75)) {
        remove_local_points <- c(remove_local_points, i + 1)
      }
    }

    if ((ns_i > 2 && ns_next < -10) ||
        (ns_i > 10 && ns_next < -2) ||
        (ns_i > 1 && ns_next < -20) ||
        (ns_i > 50 && ns_next < 0.4)) {
      remove_local_points <- c(remove_local_points, i + 1)
    }
  }

  # Last-point outlier rule.
  n <- nrow(y)
  if (n >= 2) {
    second_last <- y$parasitaemia[n - 1]
    last <- y$parasitaemia[n]

    if (second_last < 200 && last > 3 * second_last && last > 100) {
      remove_local_points <- c(remove_local_points, n)
    }
  }

  remove_local_points <- unique(remove_local_points)

  if (length(remove_local_points) > 0) {
    remove_global <- used[remove_local_points]
    x$flag[remove_global] <- "outlier"
    x$used[remove_global] <- FALSE
  }

  x
}


prepare_profile <- function(data, time_col, parasite_col) {
  if (!time_col %in% names(data)) {
    stop("`time_col` not found in data.", call. = FALSE)
  }

  if (!parasite_col %in% names(data)) {
    stop("`parasite_col` not found in data.", call. = FALSE)
  }

  tibble::tibble(
    time = as.numeric(data[[time_col]]),
    parasitaemia = as.numeric(data[[parasite_col]])
  )
}


check_estimable <- function(data,
                            initial_minimum = 1000,
                            final_maximum = 1000) {
  positive_or_replaced <- data$parasitaemia > 0

  if (sum(positive_or_replaced) < 3) {
    return(list(
      estimable = FALSE,
      reason = "fewer_than_three_nonzero_measurements"
    ))
  }

  d <- data[positive_or_replaced, , drop = FALSE]
  d <- d[order(d$time), , drop = FALSE]

  if (d$parasitaemia[1] < initial_minimum) {
    return(list(
      estimable = FALSE,
      reason = "initial_parasitaemia_too_low"
    ))
  }

  if (tail(d$parasitaemia, 1) >= final_maximum) {
    return(list(
      estimable = FALSE,
      reason = "final_parasitaemia_too_high"
    ))
  }

  list(estimable = TRUE, reason = "ok")
}


check_lag_estimable <- function(data,
                                min_first24 = 3,
                                max_gap_first24 = 14) {
  first24 <- data[data$time <= 24, , drop = FALSE]

  if (nrow(first24) < min_first24) return(FALSE)

  gaps <- diff(sort(first24$time))

  if (length(gaps) > 0 && any(gaps > max_gap_first24)) {
    return(FALSE)
  }

  TRUE
}


make_result <- function(status,
                        reason,
                        model_type = NA_character_,
                        k = NA_real_,
                        slope_half_life = NA_real_,
                        tlag = NA_real_,
                        cleaned_data = NULL,
                        fit = NULL,
                        linear_data = NULL) {
  structure(
    list(
      status = status,
      reason = reason,
      model_type = model_type,
      clearance_rate_constant = k,
      slope_half_life = slope_half_life,
      tlag = tlag,
      cleaned_data = cleaned_data,
      fit = fit,
      linear_data = linear_data
    ),
    class = "parasite_clearance_result"
  )
}


#' @export
print.parasite_clearance_result <- function(x, ...) {
  cat("Parasite clearance estimate\n")
  cat("Status:", x$status, "\n")
  cat("Reason:", x$reason, "\n")

  if (identical(x$status, "estimated")) {
    cat("Model:", x$model_type, "\n")
    cat("Lag phase:", round(x$tlag, 3), "hours\n")
    cat("Clearance rate constant:", round(x$clearance_rate_constant, 5), "per hour\n")
    cat("Slope half-life:", round(x$slope_half_life, 3), "hours\n")
  }

  invisible(x)
}

fit_linear_or_tobit <- function(data) {
  d <- add_log_parasitaemia(data)
  has_censored <- any(d$flag == "zero_replaced_detection_limit", na.rm = TRUE)

  if (has_censored) {
    model <- fit_tobit_model(d, degree = 1)
    return(list(
      model = model,
      model_type = "tobit_linear",
      degree = 1,
      aic = stats::AIC(model),
      data = d
    ))
  }

  model <- stats::lm(log_parasitaemia ~ time, data = d)

  list(
    model = model,
    model_type = "linear",
    degree = 1,
    aic = stats::AIC(model),
    data = d
  )
}


fit_candidate_models <- function(data) {
  d <- add_log_parasitaemia(data)
  n <- nrow(d)
  has_censored <- any(d$flag == "zero_replaced_detection_limit", na.rm = TRUE)

  candidates <- list()

  if (n == 3) {
    candidates$linear <- fit_model_by_degree(d, degree = 1, censored = has_censored)
    return(candidates)
  }

  if (n == 4) {
    candidates$linear <- fit_model_by_degree(d, degree = 1, censored = has_censored)
    candidates$quadratic <- fit_model_by_degree(d, degree = 2, censored = has_censored)

    if (d$parasitaemia[2] > 1.25 * d$parasitaemia[1]) {
      d2 <- d[-1, , drop = FALSE]
      candidates$maximum_regression <- fit_model_by_degree(
        d2,
        degree = 1,
        censored = has_censored,
        model_type = if (has_censored) "tobit_maximum_regression" else "maximum_regression"
      )
    }

    return(candidates)
  }

  if (n == 5 && has_censored) {
    candidates$linear <- fit_model_by_degree(d, degree = 1, censored = TRUE)
    candidates$quadratic <- fit_model_by_degree(d, degree = 2, censored = TRUE)

    if (d$parasitaemia[2] > 1.25 * d$parasitaemia[1]) {
      d2 <- d[-1, , drop = FALSE]
      candidates$maximum_regression <- fit_model_by_degree(
        d2,
        degree = 1,
        censored = TRUE,
        model_type = "tobit_maximum_regression"
      )
    }

    return(candidates)
  }

  candidates$linear <- fit_model_by_degree(d, degree = 1, censored = has_censored)
  candidates$quadratic <- fit_model_by_degree(d, degree = 2, censored = has_censored)
  candidates$cubic <- fit_model_by_degree(d, degree = 3, censored = has_censored)

  candidates
}


fit_model_by_degree <- function(d,
                                degree,
                                censored = FALSE,
                                model_type = NULL) {
  if (is.null(model_type)) {
    prefix <- if (censored) "tobit_" else ""
    model_type <- paste0(
      prefix,
      c("linear", "quadratic", "cubic")[degree]
    )
  }

  if (censored) {
    model <- fit_tobit_model(d, degree)
  } else {
    form <- polynomial_formula(degree)
    model <- stats::lm(form, data = d)
  }

  list(
    model = model,
    model_type = model_type,
    degree = degree,
    aic = stats::AIC(model),
    data = d
  )
}


select_best_model <- function(candidates) {
  aics <- vapply(candidates, function(z) z$aic, numeric(1))
  candidates[[which.min(aics)]]
}


add_log_parasitaemia <- function(data) {
  d <- data[data$parasitaemia > 0, , drop = FALSE]
  d$log_parasitaemia <- log(d$parasitaemia)
  d
}


polynomial_formula <- function(degree) {
  switch(
    as.character(degree),
    "1" = log_parasitaemia ~ time,
    "2" = log_parasitaemia ~ time + I(time^2),
    "3" = log_parasitaemia ~ time + I(time^2) + I(time^3),
    stop("Unsupported degree.", call. = FALSE)
  )
}


fit_tobit_model <- function(d, degree) {
  # Left-censored at log(detection limit) for points flagged as replaced zero.
  # survreg with type = "left" handles left-censored observations.
  #
  # For exact observations:
  #   Surv(y, y, type = "interval2")
  # For left-censored y <= limit:
  #   Surv(NA, y, type = "interval2")

  y <- d$log_parasitaemia
  censored <- d$flag == "zero_replaced_detection_limit"

  left <- y
  right <- y

  left[censored] <- NA_real_

  surv_y <- survival::Surv(left, right, type = "interval2")

  if (degree == 1) {
    form <- surv_y ~ time
  } else if (degree == 2) {
    form <- surv_y ~ time + I(time^2)
  } else if (degree == 3) {
    form <- surv_y ~ time + I(time^2) + I(time^3)
  } else {
    stop("Unsupported degree.", call. = FALSE)
  }

  survival::survreg(form, data = d, dist = "gaussian")
}


predict_fit <- function(fit, newdata = NULL) {
  if (is.null(newdata)) {
    newdata <- fit$data
  }

  as.numeric(stats::predict(fit$model, newdata = newdata, type = "response"))
}


estimate_lag_and_k <- function(data, best_fit) {
  d <- best_fit$data
  d <- d[order(d$time), , drop = FALSE]

  if (best_fit$degree == 1 || is_concave_or_linear(best_fit)) {
    linear_fit <- fit_linear_or_tobit(d)
    k <- -coef(linear_fit$model)["time"]

    return(list(
      k = as.numeric(k),
      tlag = 0,
      linear_data = d,
      linear_fit = linear_fit
    ))
  }

  lag <- identify_lag_phase(d, best_fit)

  if (lag$tlag <= 0) {
    linear_fit <- fit_linear_or_tobit(d)
    k <- -coef(linear_fit$model)["time"]

    return(list(
      k = as.numeric(k),
      tlag = 0,
      linear_data = d,
      linear_fit = linear_fit
    ))
  }

  linear_data <- d[d$time >= lag$tlag, , drop = FALSE]

  # Need at least 3 points in final linear segment.
  if (nrow(linear_data) < 3) {
    linear_fit <- fit_linear_or_tobit(d)
    k <- -coef(linear_fit$model)["time"]

    return(list(
      k = as.numeric(k),
      tlag = 0,
      linear_data = d,
      linear_fit = linear_fit
    ))
  }

  linear_fit <- fit_linear_or_tobit(linear_data)
  k <- -coef(linear_fit$model)["time"]

  list(
    k = as.numeric(k),
    tlag = lag$tlag,
    linear_data = linear_data,
    linear_fit = linear_fit
  )
}


is_concave_or_linear <- function(fit) {
  # Evaluate curvature over observed time domain.
  # For quadratic: second derivative = 2 * beta2.
  # For cubic: second derivative = 2 * beta2 + 6 * beta3 * t.
  co <- coef(fit$model)
  degree <- fit$degree

  if (degree == 1) return(TRUE)

  times <- fit$data$time

  if (degree == 2) {
    beta2 <- unname(co[grep("I\\(time\\^2\\)", names(co))])
    if (length(beta2) == 0 || !is.finite(beta2)) return(TRUE)
    return(beta2 <= 0)
  }

  if (degree == 3) {
    beta2 <- unname(co[grep("I\\(time\\^2\\)", names(co))])
    beta3 <- unname(co[grep("I\\(time\\^3\\)", names(co))])

    if (length(beta2) == 0 || length(beta3) == 0) return(TRUE)

    curvature <- 2 * beta2 + 6 * beta3 * times

    # If mostly concave/non-convex, do not estimate lag.
    return(mean(curvature > 0, na.rm = TRUE) < 0.5)
  }

  TRUE
}


identify_lag_phase <- function(d, fit) {
  pred <- predict_fit(fit, newdata = d)

  slopes <- diff(pred) / diff(d$time)

  if (length(slopes) < 2 || all(!is.finite(slopes))) {
    return(list(tlag = 0, lag_indices = integer(0)))
  }

  smax <- min(slopes, na.rm = TRUE)

  if (!is.finite(smax) || smax >= 0) {
    return(list(tlag = 0, lag_indices = integer(0)))
  }

  normalized <- slopes / smax

  # Flegg diagram logic:
  # If normalized slopes are negative or < 1/5, inspect whether this occurs
  # at the beginning of the profile only. Since smax is negative, smaller
  # normalized values represent flatter slopes relative to maximum clearance.
  lag_slope_idx <- which(normalized < 0.2 | normalized < 0)

  if (length(lag_slope_idx) == 0) {
    return(list(tlag = 0, lag_indices = integer(0)))
  }

  # Lag must be contiguous from the beginning.
  expected <- seq_len(max(lag_slope_idx))

  if (!all(expected %in% lag_slope_idx)) {
    return(list(tlag = 0, lag_indices = integer(0)))
  }

  # Lag phase ends at the time of the last measurement whose following slope
  # is flat. The boundary point is included in the final linear fit.
  tlag <- d$time[max(lag_slope_idx) + 1]

  list(
    tlag = as.numeric(tlag),
    lag_indices = seq_len(max(lag_slope_idx))
  )
}


#' Estimate parasite clearance for multiple patients
#'
#' @param data Data frame containing patient ID, time, and parasitaemia.
#' @param id_col Patient identifier column.
#' @param time_col Time column, in hours.
#' @param parasite_col Parasitaemia column.
#' @param detection_limit Detection limit. Either a scalar or a column name.
#' @param ... Additional arguments passed to [estimate_clearance()].
#'
#' @return A tibble with one row per patient.
#' @export
estimate_clearance_batch <- function(data,
                                     id_col = "id",
                                     time_col = "time",
                                     parasite_col = "parasitaemia",
                                     detection_limit = 16,
                                     ...) {
  if (!id_col %in% names(data)) {
    stop("`id_col` not found in data.", call. = FALSE)
  }

  ids <- unique(data[[id_col]])

  results <- lapply(ids, function(id) {
    d <- data[data[[id_col]] == id, , drop = FALSE]

    dl <- if (length(detection_limit) == 1 && detection_limit %in% names(d)) {
      unique(d[[detection_limit]])[1]
    } else {
      detection_limit
    }

    res <- estimate_clearance(
      data = d,
      time_col = time_col,
      parasite_col = parasite_col,
      detection_limit = dl,
      ...
    )

    tibble::tibble(
      id = id,
      status = res$status,
      reason = res$reason,
      model_type = res$model_type,
      tlag = res$tlag,
      clearance_rate_constant = res$clearance_rate_constant,
      slope_half_life = res$slope_half_life
    )
  })

  dplyr::bind_rows(results)
}




