# Script with example to calculate and plot a continuous wavelet transform of a simulated time series.
# Author: Lotte Pohl, Date: 2023-08-08

library(dplyr)
library(ggplot2)
# library(tibble)

# 0. load functions

base::paste0(base::getwd(), "/functions/compute_CWT.R") %>% base::source()
base::paste0(base::getwd(), "/functions/make_CWT_df.R") %>% base::source()
base::paste0(base::getwd(), "/functions/plot_CWT_ggplot.R") %>% base::source()
base::paste0(base::getwd(), "/functions/compute_bivariate_wavelet_analysis.R") %>% base::source()

# plot language == English
Sys.setlocale("LC_TIME", "English")

## set plot theme ####

## theme ####

plot_theme <- ggplot2::theme(
  plot.title = element_text(size = 11, face = "bold"),
  plot.subtitle = element_text(size = 11),
  axis.title = element_text(size = 11),
  axis.text = element_text(size = 9), #5.5
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 9),
  legend.key = element_rect(fill = "transparent", colour = "transparent"),
  # legend.key.width = unit(2, "cm"),
  legend.margin = margin(t = -10, b = -10, r = -10, l = -10),
  plot.tag = element_text(face = "bold", size = 12),
  # plot.tag.position =  c(0.065, 0.96), #"topleft", 
  plot.tag.position =  c(0.01, 0.98), #"topleft", 
  plot.background = element_blank(),
  panel.background = element_blank(),
  legend.box.background = element_blank(),
  legend.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
  # panel.background = element_rect(fill = "transparent"),
  panel.grid.major = element_line(color = "gray40", linetype = "solid", size = 0.5),
  panel.grid.minor = element_line(color = "gray60", linetype = "dashed", size = 0.35),
  strip.background = element_blank(),
  strip.background.x = element_blank(),
  strip.background.y = element_blank(),
  strip.text = element_text() #face = "bold", 
)

# Set the theme as the default for all plots
ggplot2::theme_set(plot_theme)

# 1. create example signals ####

f = 1/720 #because if 30 samples per hour are to be taken, these amount to 30 * 24 = 720 samples per day

signal_start_date <- "2018-07-01 01:00:00" %>% base::as.POSIXct(tz = "UTC")
signal_length_days <- 500

sampling_rate_hours <- 30 # 30 samples per hour
dt_hours <- 1/sampling_rate_hours # 2 mins expressed in the unit of hours, i.e., 2 mins / 60 mins == 1/30 hours ~ 0.033h

## 12h period ####

signal_12 <- dplyr::tibble(t = seq(0,signal_length_days - f, f),
                              depth_m = 10 * (sin((4 * pi)*t) - 1),
                              date_time = seq(signal_start_date,
                                              signal_start_date + lubridate::days(signal_length_days),
                                              length.out = base::length(t))) 

## 12 and 24h periods ####

signal_12_24 <- dplyr::tibble(t = seq(0,signal_length_days - f, f),
                 depth_m = 10 * (sin((4 * pi)*t) - 1),
                 date_time = seq(signal_start_date,
                                 signal_start_date + lubridate::days(signal_length_days),
                                 length.out = base::length(t))) %>%
  dplyr::mutate(depth_m = ifelse(date_time > signal_start_date + lubridate::days(signal_length_days / 2),
                                 10 * (sin((2 * pi)*t) - 1),
                                 depth_m))

### plot signals ####

#### 12h period ####

signal_12_plot_subset <- ggplot2::ggplot() +
  geom_line(data = signal_12 %>% dplyr::filter(date_time %>% between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(10), 
                                                                        (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(10)))
            , aes(x = date_time, y = depth_m), size = .75) + #, colour = "#357984"
  labs(x = "", y = "Depth in m") +
  scale_x_datetime(
    date_minor_breaks = "1 day",
    date_breaks = "1 week",
    date_labels = "%b %d", #  
    expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

signal_12_plot_subset + theme_bw()

#### 12 and 24h periods ####

signal_12_24_plot_subset <- ggplot2::ggplot() +
  geom_line(data = signal_12_24 %>% dplyr::filter(date_time %>% between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(10), 
                                                                  (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(10)))
            , aes(x = date_time, y = depth_m), size = .75) + #, colour = "#357984"
  labs(x = "", y = "Depth in m") +
  scale_x_datetime(
    date_minor_breaks = "1 day",
    date_breaks = "1 week",
    date_labels = "%b %d", #  
    expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

signal_12_24_plot_subset + theme_bw()

# 2. `compute_CWT()` example ####

signal_12_CWT <- compute_CWT(values = signal_12 %>% select(depth_m),
                                dt = dt_hours * 15, # every 30 min
                                factor_smallest_scale = 2) # test different vals

signal_12_24_CWT <- compute_CWT(values = signal_12_24 %>% select(depth_m),
                          dt = dt_hours * 15, # every 30 min
                          factor_smallest_scale = 2) # test different vals

# 3. `make_CWT_df()` example ####

signal_12_CWT_df <- make_CWT_df(date_times = signal_12 %>% dplyr::select(date_time),
                             dt = dt_hours * 15,
                             cwt_result = signal_12_CWT)

signal_12_24_CWT_df <- make_CWT_df(date_times = signal_12_24 %>% dplyr::select(date_time),
                                   dt = dt_hours * 15,
                                   cwt_result = signal_12_24_CWT)

# 4. `plot_CWT_ggplot()` example

signal_12_CWT_plot <- plot_CWT_ggplot(cwt_df = signal_12_CWT_df,
                                   date = T,
                                   max_period = signal_12_CWT_df %>% dplyr::select(period) %>% max()) # signal_12_CWT_df %>% dplyr::select(period) %>% max() %>% round(digits = -2)

signal_12_24_CWT_plot <- plot_CWT_ggplot(cwt_df = signal_12_24_CWT_df,
                                         date = T,
                                         max_period = signal_12_24_CWT_df %>% dplyr::select(period) %>% max()) # signal_12_24_CWT_df %>% dplyr::select(period) %>% max() %>% round(digits = -2)
signal_12_CWT_plot
signal_12_24_CWT_plot

# save plots
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/plots/signal_12_CWT_plot.png"), plot = signal_12_CWT_plot + theme_bw(base_size = 5))
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/plots/signal_12_24_CWT_plot.png"), plot = signal_12_24_CWT_plot + theme_bw(base_size = 5))

# 4. compute cross wavelet analysis ####

signals_xwt <- compute_bivariate_wavelet_analysis(type = 'cross wavelet',
                                                  values1 = signal_12 %>% 
                                                    # dplyr::filter(date_time %>%
                                                    #                 between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(30),
                                                    #                         (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(30))) %>%
                                                    select(depth_m),
                                                  values2 = signal_12_24 %>% 
                                                    # dplyr::filter(date_time %>%
                                                                    # between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(30),
                                                                    #         (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(30))) %>%
                                                    select(depth_m),
                                                  dt = dt_hours * 15,
                                                  factor_smallest_scale = 4) 

signals_xwt_df <- make_CWT_df(date_times = signal_12 %>% 
                                # dplyr::filter(date_time %>%
                                #                 between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(30),
                                #                         (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(30))) %>% 
                                dplyr::select(date_time),
                                dt = dt_hours * 15,
                                cwt_result = signals_xwt)

signals_xwt_plot <- plot_CWT_ggplot(cwt_df = signals_xwt_df,
                                      date = T,
                                      max_period = signals_xwt_df %>% dplyr::select(period) %>% max()) # signal_12_CWT_df %>% dplyr::select(period) %>% max() %>% round(digits = -2)
signals_xwt_plot
