# JT's functions:

jthomggtheme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14,
                              hjust = 0.5),
    axis.ticks = element_line(),
    # Y-axis:
    axis.title.y = element_text(angle = 0,
                                size = 12,
                                vjust = 0.5),
    axis.text.y  = element_text(size = 11),
    # X-axis:
    axis.title.x = element_text(angle = 0,
                                size = 12,
                                vjust = 0.5),
    axis.text.x  = element_text(angle = 0,
                                size = 11))


calc_rmse <- function(actual, predicted) {
  
  sqrt(mean((actual - predicted)^2))
  
}


count_na <- function(df) {
  n_rows <- nrow(df)
  df_na <- sapply(df, function(x) sum(is.na(x))) # or colSums(is.na(df))
  if (any(df_na > 0)) {
    if ("tidyverse" %in% .packages()) {
      df_na |>
        tibble::as_tibble(rownames = c("Variable")) |>
        dplyr::rename(NA_count = value) |>
        dplyr::filter(NA_count > 0) |>
        dplyr::mutate(Percent = round(100 * NA_count / n_rows, 2)) |>
        dplyr::arrange(dplyr::desc(Percent))
    } else {
      df_na <- data.frame(Variable = names(df_na), NA_count = df_na)
      df_na <- subset(df_na, df_na$NA_count > 0)
      df_na$Percent <- round(100 * df_na$NA_count / n_rows, 2)
      df_na <- df_na[order(df_na$NA_count, decreasing = TRUE)]
      rownames(df_na) <- NULL
      df_na
    }
  } else {
    print("No missing values (NA) found.")
  }
}

replace_na_with_median <- function(x) {
  ifelse(is.na(x), median(x, na.rm = TRUE), x)
}

inv_boxcox <- function(x, lambda) {
  if (lambda == 0) {
    exp(x)
  } else {
    (lambda * x + 1)^(1 / lambda)
  }
}

find_mode <- function(x) {
  unique_values <- unique(x)
  tab <- tabulate(match(x, unique_values))
  unique_values[tab == max(tab)]
}


#' Calculation of the Mutual Information (MI) criterion
#' between Response and Predictors.
#' Requires the "infotheo" package!
#'
#' @author Jorge Thomas
#' @param data a data frame with Response (Target) and Predictors (features) columns.
#' @param target a string with the column name of the Response.
#' @return a data.frame ordered with the MI calculation.
calc_mi_score <- function(data, target) {
  
  mi_list <- vector(mode = "list")
  disc_method <- "equalwidth" # "equalfreq" or "globalequalwidth"
  numbins <- sqrt(nrow(data))
  m_data <- infotheo::discretize(data.matrix(data),
                                 disc = disc_method,
                                 nbins = numbins)
  
  for (col_name in colnames(m_data)){
    
    mi_list[col_name] <- infotheo::mutinformation(m_data[, target], m_data[, col_name])
    
  }
  
  res_vec <- sort(unlist(mi_list, use.names = TRUE), decreasing = TRUE)
  
  df_mi <- data.frame(Variable = names(res_vec),
                      nats = res_vec,
                      Percent = round(100 * res_vec / res_vec[1], 4))
  
  df_mi  
}


#' Calculation of Yeo-Johnson transform for one vector:
#' 
#' 
#' 

calc_yeo_johnson <- function(y, data)  {
  
  # Specify the recipe
  rec <- recipe(~ y, data = data)  |> 
  step_YeoJohnson(y)  # Apply the Yeo-Johnson transformation to column y
  
  # Prepare the recipe
  prepped_rec <- prep(rec, training = data)
  y_yj <- bake(prepped_rec, new_data = data)
  
  y_yj
  
}


#' Find and Compare Best-Fit Distributions for Repair Times
#'
#' Fits Weibull, Lognormal, and Exponential distributions to a vector of repair times (duration),
#' compares their goodness-of-fit, and plots the empirical and fitted distributions.
#'
#' @param turbine A vector or factor indicating the turbine (for labeling purposes).
#' @param duration A numeric vector of repair times (e.g., Time To Repair or Time To Failure).
#'
#' @details
#' The function fits Weibull, Lognormal, and Exponential distributions to the provided duration data
#' using the \code{fitdistrplus} package. It then overlays the fitted probability density functions
#' on a histogram and kernel density estimate of the empirical data. Goodness-of-fit statistics are printed
#' to help select the best distribution.
#'
#' @return Invisibly returns the goodness-of-fit statistics (as a list) and displays a plot.
#' @import ggplot2
#' @importFrom fitdistrplus fitdist gofstat
#' @examples
#' \dontrun{
#' find_best_fit(turbine = my_data$Turbine, duration = my_data$TTR)
#' }
#' @export
find_best_fit <- function(turbine, duration) {

  # paste("Fit for", turbine, "with", length(duration), "samples.")

  df <- data.frame(turbine, duration)
  # Fit distributions
  fit_weibull <- fitdistrplus::fitdist(df$duration, "weibull")
  fit_lnorm <- fitdistrplus::fitdist(df$duration, "lnorm")
  fit_exp <- fitdistrplus::fitdist(df$duration, "exp")

  # Plot empirical vs theoretical distributions
  # plt <- ggplot(df, aes(x = duration)) +
  #   geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "#1E88E5", alpha = 0.7) +
  #   geom_density(color = "#FF5252", linewidth = 1) +
  #   stat_function(fun = dweibull, 
  #                 args = list(shape = fit_weibull$estimate[1], 
  #                             scale = fit_weibull$estimate[2]),
  #                 color = "darkgreen", linetype = "dashed") +
  #   stat_function(fun = dlnorm,
  #                 args = list(meanlog = fit_lnorm$estimate[1],
  #                             sdlog = fit_lnorm$estimate[2]),
  #                 color = "blue", linetype = "dashed") +
  #   stat_function(fun = dexp,
  #                 args = list(rate = fit_exp$estimate[1]),
  #                 color = "orange", linetype = "dashed") +

  #   labs(title = paste("Distribution of Duration Samples with Fitted Curves", df$turbine[1]),
  #        x = "Duration (hours)",
  #        y = "Density") +
  #   theme_minimal() 
  
  gof_stats <- fitdistrplus::gofstat(list(fit_weibull, fit_lnorm, fit_exp))

  return(gof_stats)
}


#' Monte Carlo Chronological Simulation of Failure and Repair Times
#'
#' Simulates a chronological sequence of alternating failure (up) and repair (down) events for a wind turbine,
#' using exponential and lognormal distributions for time-to-failure (TTF) and time-to-repair (TTR), respectively.
#'
#' @param sim_parameters A list containing the simulation parameters: \code{frate} (failure rate, numeric, [h^-1]),
#'   \code{meanlog} (meanlog parameter for lognormal repair time), and \code{sdlog} (sdlog parameter for lognormal repair time).
#' @param n_events Integer. Number of failure-repair cycles to simulate.
#' @param dt_start character. Start date and time for the simulation in "YYYY-MM-DD HH" format.
#' @param  truncate_dist Logical. If TRUE, truncates the distributions to avoid unrealisticly short durations. Ideal for combining  components.
#' @param ttf_min_duration Numeric. Minimum duration for time-to-failure (TTF) in hours.
#' @param ttr_min_duration Numeric. Minimum duration for time-to-repair (TTR) in hours.
#' @details
#' The function generates \code{n_events} time-to-failure values from an exponential distribution with rate \code{frate},
#' and \code{n_events} time-to-repair values from a lognormal distribution with parameters \code{meanlog} and \code{sdlog}.
#' It then intercalates these into a chronological sequence, builds a timestamped data frame starting from "2023-01-01 00",
#' and assigns a state column (1 = up, 0 = down) for each interval.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{timestamp}: POSIXct, start time of each interval
#'     \item \code{state}: integer, 1 (operating) or 0 (repair)
#'     \item \code{duration}: numeric, duration of the interval in hours
#'   }
#'
#' @examples
#' \dontrun{
#' sim_parameters <- list(frate = 0.01, meanlog = 1, sdlog = 0.5)
#' df_sim <- mcs_events(sim_parameters, n_events = 100)
#' head(df_sim)
#' }
#' @importFrom dplyr tibble select
#' @importFrom lubridate ymd_h
#' @imoprtFrom truncdist rtrunc
#' @export
mcs_events <- function(sim_parameters, n_events, dt_start, truncate_dist = FALSE, ttf_min_duration = 1/6, ttr_min_duration = 1/6) {
  # Check if the required parameters are provided
  if (missing(sim_parameters) || missing(n_events)) {
    stop("Please provide both sim_parameters and n_events.")
  }

  # Check if sim_parameters is a list and contains the required elements
  if (!is.list(sim_parameters) || !all(c("frate", "meanlog", "sdlog") %in% names(sim_parameters))) {
    stop("sim_parameters must be a list containing frate, meanlog, and sdlog.")
  }

  # Extract parameters from the list
  lambda <- sim_parameters$frate      # Failure rate [h^-1]
  meanlog <- sim_parameters$meanlog   # Lognormal meanlog for TTR
  sdlog <- sim_parameters$sdlog       # Lognormal sdlog for TTR

  set.seed(1982)  # For reproducibility

  if (truncate_dist) {
    # Simulate time-to-failure (exponential distribution)
    ttf_sim <- truncdist::rtrunc(n = n_events,
                                 spec = "exp",
                                 a = ttf_min_duration,
                                 b = Inf,
                                 rate = lambda)
    # Simulate time-to-repair (lognormal distribution)
    ttr_sim <- truncdist::rtrunc(n = n_events,
                                 spec = "lnorm",
                                 a = ttr_min_duration,
                                 b = Inf,
                                 meanlog = meanlog,
                                 sdlog = sdlog)

  } else {
    # Simulate time-to-failure (exponential distribution)
    ttf_sim <- rexp(n_events, rate = lambda)
    # Simulate time-to-repair (lognormal distribution)
    ttr_sim <- rlnorm(n_events, meanlog = meanlog, sdlog = sdlog)
  }

  df_sim <- dplyr::tibble(
    event = seq(1:n_events),
    ttf = ttf_sim,
    ttr = ttr_sim,
    failure_time = cumsum(ttf_sim),
    total_time = cumsum(ttf_sim + ttr_sim)
  )

  df_sim$A <- cumsum(df_sim$ttf) / df_sim$total_time
  df_sim$U <- 1-df_sim$A 

  # Intercalate ttf and ttr into a single vector
  ttf_ttr_intercalated <- as.vector(rbind(df_sim$ttf, df_sim$ttr))
  A_U_intercalated <- as.vector(rbind(df_sim$A, df_sim$U))

  # Build a data frame
  results <- dplyr::tibble(
    duration = ttf_ttr_intercalated,
    state = rep(c(1, 0), length.out = length(ttf_ttr_intercalated)),
    timestamp = lubridate::ymd_h(dt_start) + seconds(c(0, cumsum(duration[-length(duration)])*3600)),
    A_U = A_U_intercalated
  )

  results <- results |> dplyr::select(timestamp, state, duration, A_U)
  results
}


# get mode function
# df_all |>
#   count(MasVnrType, name = "Count", sort = TRUE) |>
#   filter(Count == max(Count, na.rm = TRUE)) |>
#   pull(MasVnrType)