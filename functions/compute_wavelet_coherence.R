#' wrapper function of `biwavelet::wt()` to compute the continuous wavelet transform of a data frame with containing a one dimensional time series
#'
#' @author Lotte Pohl (lotte.pohl@@gmail.com)
#'
#' @param values numeric vector of values to be taken as input for the CWT calculation.
#' @param dt The desired time step between two values that are taken as input. Can either be the original dt between two values or a larger time distance, then the values are downsampled. Has to be in the unit of the desired periods of the CWT, i.e. hours, days, minutes, years, ...
#' @param factor_smallest_scale The smallest scale of changes that are to be detected by the CWT. Alternative description: Resolution.


library(biwavelet)
library(dplyr)
library(gridExtra)

# test ###

values1 <- signal_12 %>% 
  dplyr::filter(date_time %>% 
                  between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(10),
                          (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(10))) %>%
  select(depth_m)

values2 <- signal_12_24 %>% 
  dplyr::filter(date_time %>% 
                  between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(10), 
                          (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(10))) %>%
  select(depth_m)

dt <- dt_hours * 30 # every 600 min
factor_smallest_scale <- 4

# 1. compute wavelet coherence ####

compute_wavelet_coherence <- function(values1, values2, dt, factor_smallest_scale = 8){
  # To Do: test that class(values) %in% c("tbl_df", "tbl", data.frame")
  
  # downsample raw values
  values1_downsampled <- values1 %>% 
    dplyr::ungroup() %>%
    dplyr::filter(dplyr::row_number() %% ((dt * 60) / 2) == 0) # change 60 to allow multiple time scales, explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt_smallperiods = 1: ...%% 30 (i.e. get every 30th val)
  
  values2_downsampled <- values2 %>% 
    dplyr::ungroup() %>%
    dplyr::filter(dplyr::row_number() %% ((dt * 60) / 2) == 0) # change 60 to allow multiple time scales, explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt_smallperiods = 1: ...%% 30 (i.e. get every 30th val)
  
  # make time vector according to nrow() of parameter and given dt
  timevector <- base::seq(from = 0, to = (base::nrow(values1_downsampled) * dt) - dt, by = dt)
  wtc_input1 <- base::cbind(timevector, values1_downsampled) %>% base::as.matrix()
  wtc_input2 <- base::cbind(timevector, values2_downsampled) %>% base::as.matrix()
  # make wt result
  wtc <- biwavelet::wtc(d1 = wtc_input1,
                        d2 = wtc_input2,
                        s0 = factor_smallest_scale * dt) # this specifies the scale at which periods are looked at
  
  return(cwt)
}

plot(wtc, plot.phase = TRUE)

gridExtra::grid.arrange(signal_12_CWT_plot, signal_12_24_CWT_plot, plot(wtc, plot.phase = TRUE))

