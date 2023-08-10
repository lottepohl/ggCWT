#' function to transform an output of `biwavelet::wt()` into a data frame that can be taken as input for a `ggplot2` object
#' 
#' @author Lotte Pohl (lotte.pohl@@gmail.com)
#'
#' `values` and `date_times` must have the same length, i.e. one value should correspond to one element of `date_times`.
#'
#' @param values numeric vector of values to be taken as input for the wavelet_result calculation.
#' @param dt The desired time step between two values that are taken as input. Can either be the original dt between two values or a larger time distance, then the values are downsampled. Has to be in the unit of the desired periods of the wavelet_result, i.e. hours, days, minutes, years, ...
#' @param date_times vector containing date and time that correspond to the values.
#' @param object of class "biwavelet" from `biwavelet::wt()`.


library(biwavelet)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

##### test ###
# values <- signal_12_24 %>% select(depth_m)
# dt <- dt_hours * 15
# date_times <- signal_12_24 %>% select(date_time)
# wavelet_result <- signal_12_24_cwt

make_wavelet_df <- function(wavelet_result, date_times, dt){

  # downsample date_times
  date_times_downsampled <- date_times %>% 
    dplyr::ungroup() %>%
    dplyr::filter(dplyr::row_number() %% ((dt * 60) / 2) == 0) %>% # change 60 to allow multiple time scales, explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt_smallperiods = 1: ...%% 30 (i.e. get every 30th val)
    `colnames<-`("date_times_downsampled") # round brackets == obligatory
  
  # make time vector according to nrow() of parameter and given dt
  timevector <- base::seq(from = 0, to = (nrow(date_times_downsampled) * dt) - dt, by = dt)
  
  dates <- date_times_downsampled %>% 
    dplyr::mutate(t = timevector,
                  date_times_character= date_times_downsampled %>% as.character())
  
  # extract important results
  period <- wavelet_result$period %>% base::as.data.frame() %>% `colnames<-`("period")
  xaxis <- wavelet_result$xaxis %>% base::as.data.frame() %>% `colnames<-`("time")
  signif <- wavelet_result$signif %>% base::as.data.frame() %>%
    purrr::set_names(dates$date_times_character) %>%
    base::cbind(period) %>%
    tidyr::pivot_longer(cols = -last_col(offset = 0), names_to = "date_times_character") %>% #don't pivot the last column
    dplyr::rename(significance = value) 
  
  wavelet_df <- wavelet_result$power %>% base::as.data.frame() %>%
    purrr::set_names(dates$date_times_character) %>%
    base::cbind(period) %>%
    dplyr::arrange(desc(period))
  
  wavelet_df <- wavelet_df %>%
    tidyr::pivot_longer(cols = -last_col(offset = 0), names_to = "date_times_character") %>% 
    dplyr::rename(power = value) %>% 
    dplyr::left_join(dates, by = "date_times_character", multiple = "all") %>%
    dplyr::left_join(signif, by = c("period", "date_times_character"), multiple = "all") %>% #View()
    dplyr::mutate(date = date_times_downsampled %>% lubridate::date() %>% as.POSIXct.Date()) %>%
    dplyr::group_by(date, period) %>%
    dplyr::summarise(power = power %>% max(),
                     significance = significance %>% max(),
                     t = t %>% max()) %>%
    dplyr::mutate(power_log = log2(power),
                  power_log_scale = power_log %>% scale(),
                  period_log = log2(period),
                  t = sprintf("%04f", t %>% as.numeric()),
                  sig = ifelse(significance >= 1, 1, 0)) %>%
    dplyr::ungroup()
  
  return(wavelet_df)
}

# make_wavelet_df <- function(values, date_times, dt, wavelet_result){
#   # To Do: test that class(values) %in% c("tbl_df", "tbl", data.frame")
#   
#   # downsample raw values
#   values_downsampled <- values %>% 
#     dplyr::ungroup() %>%
#     dplyr::filter(dplyr::row_number() %% ((dt * 60) / 2) == 0) # change 60 to allow multiple time scales, explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt_smallperiods = 1: ...%% 30 (i.e. get every 30th val)
#  
#   # downsample date_times
#   date_times_downsampled <- date_times %>% 
#     dplyr::ungroup() %>%
#     dplyr::filter(dplyr::row_number() %% ((dt * 60) / 2) == 0) # change 60 to allow multiple time scales, explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt_smallperiods = 1: ...%% 30 (i.e. get every 30th val)
#   
#   # combine downsampled date_times and values
#   data_downsampled <- dplyr::tibble(values_downsampled,
#                                     date_times_downsampled) %>%
#     `colnames<-`(c("values_downsampled", "date_times_downsampled"))
#   
#   
#   # make time vector according to nrow() of parameter and given dt
#   timevector <- base::seq(from = 0, to = (nrow(data_downsampled) * dt) - dt, by = dt)
#   
#   dates <- data_downsampled %>% dplyr::select(date_times_downsampled) %>% 
#     dplyr::mutate(t = timevector,
#            date = date_times_downsampled %>% lubridate::date() %>% as.POSIXct.Date(), #sum by day
#            date_times_character= date_times_downsampled %>% as.character())
#   
#   # extract important results
#   period <- wavelet_result$period %>% base::as.data.frame() %>% `colnames<-`("period")
#   xaxis <- wavelet_result$xaxis %>% base::as.data.frame() %>% `colnames<-`("time")
#   signif <- wavelet_result$signif %>% base::as.data.frame() %>%
#     purrr::set_names(base::as.character(data_downsampled$date_times_downsampled)) %>%
#     # purrr::set_names(as.character(date_times$date)) %>% #sum by day, but not yet
#     base::cbind(period) %>%
#     tidyr::pivot_longer(cols = -last_col(offset = 0), names_to = "date_times_character") %>% #don't pivot the last column
#     dplyr::rename(significance = value) #%>%
#     # mutate(date_times_downsampled = date_times_downsampled %>% base::as.POSIXct())
#   
#   # wavelet_df <- wavelet_result$power.corr %>% base::as.data.frame() %>%
#   wavelet_df <- wavelet_result$power %>% base::as.data.frame() %>%
#     purrr::set_names(base::as.character(data_downsampled$date_times_downsampled)) %>%
#     # purrr::set_names(as.character(date_times$date)) %>% #sum by day
#     base::cbind(period) %>%
#     # mutate(period_log = log2(period))# %>%
#     dplyr::arrange(desc(period))
#   
#   wavelet_df <- wavelet_df %>%
#     # tidyr::pivot_longer(cols = -c(last_col(offset = 1), last_col(offset = 0)), names_to = "date_time") %>% #don't pivot the two last columns
#     tidyr::pivot_longer(cols = -last_col(offset = 0), names_to = "date_times_character") %>% 
#     dplyr::rename(power = value) %>% 
#     # relocate(date, period, power) %>%
#     dplyr::left_join(dates, by = "date_times_character", multiple = "all") %>%
#     dplyr::left_join(signif, by = c("period", "date_times_character"), multiple = "all") %>% #View()
#     dplyr::group_by(date, period) %>%
#     dplyr::summarise(power = power %>% max(),
#               significance = significance %>% max(),
#               t = t %>% max()) %>%
#     dplyr::mutate(power_log = log2(power),
#            power_log_scale = power_log %>% scale(),
#            period_log = log2(period),
#            t = sprintf("%04f", t %>% as.numeric()),
#            sig = ifelse(significance >= 1, 1, 0)) %>%
#     dplyr::ungroup()
#   
#   return(wavelet_df)
# }
