#' wrapper function of `biwavelet::wt()` to compute the continuous wavelet transform of a data frame with containing a one dimensional time series
#'
#' @author Lotte Pohl (lotte.pohl@@gmail.com)
#'
#' @param values numeric vector of values to be taken as input for the CWT calculation.
#' @param dt The desired time step between two values that are taken as input. Can either be the original dt between two values or a larger time distance, then the values are downsampled. Has to be in the unit of the desired periods of the CWT, i.e. hours, days, minutes, years, ...
#' @param factor_smallest_scale The smallest scale of changes that are to be detected by the CWT. Alternative description: Resolution.


library(biwavelet)


# 1. compute CWT ####

compute_CWT <- function(values, dt, factor_smallest_scale = 8){
  # To Do: test that class(values) %in% c("tbl_df", "tbl", data.frame")
  # downsample raw values
  values_downsampled <- values %>% 
    dplyr::ungroup() %>%
    dplyr::filter(dplyr::row_number() %% ((dt * 60) / 2) == 0) # change 60 to allow multiple time scales, explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt_smallperiods = 1: ...%% 30 (i.e. get every 30th val)
  
  # make time vector according to nrow() of parameter and given dt
  timevector <- base::seq(from = 0, to = (base::nrow(values_downsampled) * dt) - dt, by = dt)
  wt_input <- base::cbind(timevector, values_downsampled) %>% base::as.matrix()
  # make wt result
  cwt <- biwavelet::wt(d = wt_input,
                       dt = dt,
                       do.sig = T,
                       s0 = factor_smallest_scale * dt) # this specifies the scale at which periods are looked at
  
  return(cwt)
}
