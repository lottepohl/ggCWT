# Script to carry out cross wavelet analysis of data storage tag depth logs from two *Mustelus asterias* individuals from 2018/2019.
# Author: Lotte Pohl
# Date: 2023-08-10

library(dplyr)
library(scales)
library(tidyr)
library(ggplot2)
library(biwavelet)

# load functions ####

base::paste0(base::getwd(), "/functions/compute_CWT.R") %>% base::source()
base::paste0(base::getwd(), "/functions/make_wavelet_df.R") %>% base::source()
base::paste0(base::getwd(), "/functions/ggplot_wavelet.R") %>% base::source()
base::paste0(base::getwd(), "/functions/compute_bivariate_wavelet_analysis.R") %>% base::source()

# generic signals ####
# sample rate same as the data from *Mustelus asterias*, i.e. one value every 2 min

# VERSTEHEN: warum f in days? Warum nicht hours?
f = 1/720 #because if 30 samples per hour are to be taken, these amount to 30 * 24 = 720 samples per day

signal_start_date <- "2018-07-01 01:00:00" %>% base::as.POSIXct(tz = "UTC")
signal_length_days <- 500

sampling_rate_hours <- 30 # 30 samples per hour
dt_hours <- 1/sampling_rate_hours # 2 mins expressed in the unit of hours, i.e., 2 mins / 60 mins == 1/30 hours ~ 0.033h

## 12 h period: RESTING ####

signal_12 <- dplyr::tibble(t = seq(0,signal_length_days - f, f),
                              depth_m = 10 * (sin((4 * pi)*t) - 1),
                              date_time = seq(signal_start_date,
                                              signal_start_date + lubridate::days(signal_length_days),
                                              length.out = base::length(t)))

## 24h period: DIEL VERTICAL MIGRATION ####

signal_24 <- dplyr::tibble(t = seq(0,signal_length_days - f, f),
                              depth_m = 10 * (sin((2 * pi)*t) - 1),
                              date_time = seq(signal_start_date,
                                              signal_start_date + lubridate::days(signal_length_days),
                                              length.out = base::length(t)))

## combine generic signals ####

generic_signals <- signal_12 %>%
  dplyr::rename(depth_m_12 = depth_m) %>%
  dplyr::left_join(signal_24, by = c("t", "date_time")) %>%
  dplyr::rename(depth_m_24 = depth_m)

rm(signal_12, signal_24)

## plot signals ####
generic_signals_plot_subset <- ggplot2::ggplot(data = generic_signals %>% dplyr::filter(date_time %>% between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(10), 
                                                                                                              (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(10)))) +
  geom_line(aes(x = date_time, y = depth_m_12), size = .75, colour = 'darkorange') + #, colour = "#357984"
  geom_line(aes(x = date_time, y = depth_m_24), size = .75, colour = 'darkgreen') +
  labs(x = "", y = "Depth in m") +
  scale_x_datetime(
    date_minor_breaks = "1 day",
    date_breaks = "1 week",
    date_labels = "%b %d", #  
    expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

generic_signals_plot_subset + theme_bw()

# test: XWT of 12 and 24h signals

xwt_12_24_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                          values1 = generic_signals %>% dplyr::select(depth_m_12),
                                                          values2 = generic_signals %>% dplyr::select(depth_m_24),
                                                          dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = generic_signals %>% dplyr::select(date_time),
                  dt = 15 * dt_hours) %>%
  ggplot_wavelet(max_period = max(.$period)) # TODO: understand why argument for max_period is needed
  
xwt_12_24_plot #%>% View()

cwt_12_plot <- compute_CWT(values = generic_signals %>% dplyr::select(depth_m_12),
                           dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = generic_signals %>% dplyr::select(date_time),
                  dt = 15 * dt_hours) %>%
  ggplot_wavelet(max_period = max(.$period))

cwt_12_plot

cwt_24_plot <- compute_CWT(values = generic_signals %>% dplyr::select(depth_m_24),
                           dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = generic_signals %>% dplyr::select(date_time),
                  dt = 15 * dt_hours) %>%
  ggplot_wavelet(max_period = max(.$period))

cwt_24_plot

gridExtra::grid.arrange(cwt_12_plot, cwt_24_plot, xwt_12_24_plot, ncol = 1)
