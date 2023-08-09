# Script with example to calculate and plot a continuous wavelet transform of a simulated time series.
# Author: Lotte Pohl, Date: 2023-08-08

library(dplyr)
library(ggplot2)
# library(tibble)

# load functions

base::paste0(base::getwd(), "/functions/compute_CWT.R") %>% base::source()

# plot language == English
Sys.setlocale("LC_TIME", "English")

# 1. create example signal ####

## 12 and 24h frequencies ####

f = 1/720 #because if 30 samples per hour are to be taken, these amount to 30 * 24 = 720 samples per day

signal_start_date <- "2018-07-01 01:00:00" %>% base::as.POSIXct(tz = "UTC")
signal_length_days <- 500

sampling_rate_hours <- 30 # 30 samples per hour
dt_hours <- 1/sampling_rate_hours #2 mins expressed in the unit of hours, i.e., 2 mins / 60 mins == 1/30 hours ~ 0.033h


signal <- dplyr::tibble(t = seq(0,signal_length_days - f, f),
                 depth_m = 10 * (sin((4 * pi)*t) - 1),
                 date_time = seq(signal_start_date,
                                 signal_start_date + lubridate::days(signal_length_days),
                                 length.out = base::length(t))) %>%
  dplyr::mutate(depth_m = ifelse(date_time > signal_start_date + lubridate::days(signal_length_days / 2),
                                 10 * (sin((2 * pi)*t) - 1),
                                 depth_m))


### plot signal ####
signal_plot_subset <- ggplot2::ggplot() +
  geom_line(data = signal %>% dplyr::filter(date_time %>% between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(10), 
                                                                  (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(10)))
            , aes(x = date_time, y = depth_m), size = .75) + #, colour = "#357984"
  labs(x = "", y = "Depth in m") +
  scale_x_datetime(
    date_minor_breaks = "1 day",
    date_breaks = "1 week",
    date_labels = "%b %d", #  
    expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

signal_plot_subset + theme_bw()

# 2. `compute_CWT` example

signal_CWT <- compute_CWT(values = signal %>% select(depth_m),
                          dt = dt_hours * 15, # every 30 min
                          factor_smallest_scale = 8) # test different vals

#### test ###
values <- signal %>% select(depth_m)
dt <- dt_hours * 15 # every 30 min