# Script with function to transform an output of `biwavelet::wt()` into a data frame that can be taken as input for a `ggplot` object
# Author: Lotte Pohl, Date: 2023-08-08

library(biwavelet)
library(dplyr)
library(purrr)
library(ggplot2)

wavelet_output_compare_hrperiod <- function(depthlog = masterias_depth_temp, tag_serial_num = "1293308", parameter_name = "depth_m", dt = dt, cwt = cwt){
  
  # downsample raw depthlog
  depthlog_downsampled <- depthlog %>% 
    dplyr::ungroup() %>%
    dplyr::filter(tag_serial_number == tag_serial_num,
                  dplyr::row_number() %% ((dt * 60) / 2) == 0) %>%
    dplyr::select(c(parameter_name, date_time))# explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt = 1: ...%% 30 (i.e. get every 30th val)
  
  # make time vector according to nrow() of parameter and given dt
  timevector <- base::seq(from = 0, to = (nrow(depthlog_downsampled) * dt) - dt, by = dt)
  
  dates <- depthlog_downsampled %>% dplyr::select(date_time) %>% 
    dplyr::mutate(t = timevector,
           date = date_time %>% lubridate::date() %>% as.POSIXct.Date(), #sum by day
           date_posicxt = date_time,
           date_time = date_time %>% as.character())
  # extract important results
  period <- cwt$period %>% as.data.frame() %>% `colnames<-`("period")
  xaxis <- cwt$xaxis %>% as.data.frame() %>% `colnames<-`("time")
  signif <- cwt$signif %>% as.data.frame() %>%
    purrr::set_names(as.character(dates$date_time)) %>%
    # purrr::set_names(as.character(dates$date)) %>% #sum by day, but not yet
    cbind(period) %>%
    pivot_longer(cols = -last_col(offset = 0), names_to = "date_time") %>% #don't pivot the two last columns
    rename(significance = value)
  # wt_df <- cwt$power.corr %>% as.data.frame() %>%
  wt_df <- cwt$power %>% as.data.frame() %>%
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
    left_join(signif, by = c("period", "date_time"), multiple = "all") %>%
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
