# Script to conduct and plot wavelet analyses

library(biwavelet)
library(dplyr)
library(plotly)
library(tibble)
library(tidyverse)
library(purrr)
library(RColorBrewer)

# rm(list = ls())

# paste0(getwd(), "/01_code/02_load_data/load_depth_temp_logs.R") %>% base::source()
# source(paste0(getwd(), "/01_code/02_load_data/load_dst_summarystatistics.R"))

# daily summaries ####

compute_wavelet <- function(parameter, dt, factor_smallest_scale){
  # make time vector according to nrow() of parameter and given dt
  timevector <- seq(from = 0, to = (nrow(parameter) * dt) - dt, by = dt)
  wt_input <- cbind(timevector, parameter) %>% as.matrix()
  # make wt result
  wt_output <- biwavelet::wt(d = wt_input,
                             dt = dt,
                             do.sig = T,
                             s0 = factor_smallest_scale * dt) # this specifies the scale at which periods are looked at
  
  return(wt_output)
}
# test with 321 date vector
# dates <- long_dst_date %>% dplyr::filter(tag_serial_number == "1293321") %>% dplyr::select(date) 
# t1 <- cbind(1:100, rnorm(100))
# wt_output <- t1 %>% wt()
# dates <- t1[,1] %>% as.data.frame()

wavelet_output_compare <- function(dates, wt_output){
  dates <- dates %>% `colnames<-`("date") %>% mutate(t = seq(from = 1, to = nrow(dates)),
                                                     date_posicxt = date,
                                                     date = date %>% as.character())
  # extract important results
  period <- wt_output$period %>% as.data.frame() %>% `colnames<-`("period")
  xaxis <- wt_output$xaxis %>% as.data.frame() %>% `colnames<-`("time")
  signif <- wt_output$signif %>% as.data.frame() %>%
    purrr::set_names(as.character(dates$date)) %>%
    cbind(period) %>%
    pivot_longer(cols = -last_col(offset = 0), names_to = "date") %>% #don't pivot the two last columns
    rename(significance = value)
  # wt_df <- wt_output$power.corr %>% as.data.frame() %>%
  wt_df <- wt_output$power %>% as.data.frame() %>%
    purrr::set_names(as.character(dates$date)) %>%
    cbind(period) %>%
    # mutate(period_log = log2(period))# %>%
    arrange(desc(period)) %>%
    mutate(height = period - dplyr::lead(period), #height = log2(period) - dplyr::lead(log2(period)
           height = height + 0.15) 
  wt_df$height[wt_df$height %>% is.na()] <- wt_df$period[nrow(wt_df) - 1] - wt_df$period[nrow(wt_df)] # fill the NA created by `dplyr::lead()`
  wt_df <- wt_df %>%
    pivot_longer(cols = -c(last_col(offset = 1), last_col(offset = 0)), names_to = "date") %>% #don't pivot the two last columns
    rename(power = value) %>%
    # relocate(date, period, power) %>%
    left_join(signif, by = join_by(period, date), multiple = "all") %>%
    mutate(date = date, # %>% as.POSIXct()
           # power_scaled = power %>% scale(),
           power_log = log2(power),
           power_log_scale = power_log %>% scale(),
           period_log = log2(period)) %>%
    left_join(dates, by = "date", multiple = "all") %>%
    mutate(t = sprintf("%03d", t %>% as.numeric()),
           sig = ifelse(significance >= 1, 1, 0),
           height2 = period - dplyr::lead(period), #height = log2(period) - dplyr::lead(log2(period)
           height2 = height2 + 0.15)
  wt_df$height2[wt_df$height2 %>% is.na()] <- wt_df$period[nrow(wt_df) - 1] - wt_df$period[nrow(wt_df)] # fill the NA created by `dplyr::lead()`

  return(wt_df)
}

# raw depth log ####

compute_wavelet_hrperiod <- function(depthlog = masterias_depth_temp, tag_serial_num = "1293308", parameter_name = "depth_m", dt = dt_smallperiods, factor_smallest_scale = 8){
  # downsample raw depthlog
  depthlog_downsampled <- depthlog %>% 
    ungroup() %>%
    dplyr::filter(tag_serial_number == tag_serial_num,
           row_number() %% ((dt_smallperiods * 60) / 2) == 0) %>%
    dplyr::select(parameter_name)# explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt_smallperiods = 1: ...%% 30 (i.e. get every 30th val)
  
  # make time vector according to nrow() of parameter and given dt
  timevector <- seq(from = 0, to = (nrow(depthlog_downsampled) * dt) - dt, by = dt)
  wt_input <- cbind(timevector, depthlog_downsampled) %>% as.matrix()
  # make wt result
  wt_output <- biwavelet::wt(d = wt_input,
                             dt = dt,
                             do.sig = T,
                             s0 = factor_smallest_scale * dt) # this specifies the scale at which periods are looked at
  
  return(wt_output)
}

wavelet_output_compare_hrperiod <- function(depthlog = masterias_depth_temp, tag_serial_num = "1293308", parameter_name = "depth_m", dt = dt_smallperiods, wt_output = wt_output){
  
  # downsample raw depthlog
  depthlog_downsampled <- depthlog %>% 
    ungroup() %>%
    dplyr::filter(tag_serial_number == tag_serial_num,
           row_number() %% ((dt_smallperiods * 60) / 2) == 0) %>%
    dplyr::select(c(parameter_name, date_time))# explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt_smallperiods = 1: ...%% 30 (i.e. get every 30th val)
  
  # make time vector according to nrow() of parameter and given dt
  timevector <- seq(from = 0, to = (nrow(depthlog_downsampled) * dt) - dt, by = dt)
  
  dates <- depthlog_downsampled %>% dplyr::select(date_time) %>% 
    mutate(t = timevector,
           date = date_time %>% lubridate::date() %>% as.POSIXct.Date(), #sum by day
           date_posicxt = date_time,
           date_time = date_time %>% as.character())
  # extract important results
  period <- wt_output$period %>% as.data.frame() %>% `colnames<-`("period")
  xaxis <- wt_output$xaxis %>% as.data.frame() %>% `colnames<-`("time")
  signif <- wt_output$signif %>% as.data.frame() %>%
    purrr::set_names(as.character(dates$date_time)) %>%
    # purrr::set_names(as.character(dates$date)) %>% #sum by day, hier noch nicht
    cbind(period) %>%
    pivot_longer(cols = -last_col(offset = 0), names_to = "date_time") %>% #don't pivot the two last columns
    rename(significance = value)
  # wt_df <- wt_output$power.corr %>% as.data.frame() %>%
  wt_df <- wt_output$power %>% as.data.frame() %>%
    purrr::set_names(as.character(dates$date_time)) %>%
    # purrr::set_names(as.character(dates$date)) %>% #sum by day
    cbind(period) %>%
    # mutate(period_log = log2(period))# %>%
    arrange(desc(period))
  wt_df <- wt_df %>%
    # pivot_longer(cols = -c(last_col(offset = 1), last_col(offset = 0)), names_to = "date_time") %>% #don't pivot the two last columns
    pivot_longer(cols = -last_col(offset = 0), names_to = "date_time") %>%
    rename(power = value) %>%
    # relocate(date, period, power) %>%
    left_join(signif, by = join_by(period, date_time), multiple = "all") %>%
    left_join(dates, by = "date_time", multiple = "all") %>%
    group_by(date, period) %>%
    summarise(power = power %>% max(),
              significance = significance %>% max(),
              t = t %>% max()) %>%
    mutate(
      # date_time = date_time, # %>% as.POSIXct()
      # power_scaled = power %>% scale(),
      power_log = log2(power),
      power_log_scale = power_log %>% scale(),
      period_log = log2(period)) %>%
    mutate(t = sprintf("%04f", t %>% as.numeric()),
           sig = ifelse(significance >= 1, 1, 0))
  
  return(wt_df)
}
