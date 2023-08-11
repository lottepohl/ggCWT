# Script to carry out cross wavelet analysis of data storage tag depth logs from two *Mustelus asterias* individuals from 2018/2019.
# Author: Lotte Pohl
# Date: 2023-08-10

library(dplyr)
library(scales)
library(tidyr)
library(ggplot2)
library(biwavelet)
# library(signal) # signals not filtered, does not really change anything

# load functions ####

base::paste0(base::getwd(), "/functions/load_save_data_rds.R") %>% base::source()
base::paste0(base::getwd(), "/functions/compute_CWT.R") %>% base::source()
base::paste0(base::getwd(), "/functions/make_wavelet_df.R") %>% base::source()
base::paste0(base::getwd(), "/functions/ggplot_wavelet.R") %>% base::source()
base::paste0(base::getwd(), "/functions/compute_bivariate_wavelet_analysis.R") %>% base::source()
base::paste0(base::getwd(), "/functions/calc_plot_fft.R") %>% base::source()


# plot language == English
Sys.setlocale("LC_TIME", "English")

## set plot theme ####

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


# generic signals ####
# sample rate same as the data from *Mustelus asterias*, i.e. one value every 2 min

# VERSTEHEN: warum f in days? Warum nicht hours?
f = 1/720 #because if 30 samples per hour are to be taken, these amount to 30 * 24 = 720 samples per day

signal_start_date <- "2018-07-01 01:00:00" %>% base::as.POSIXct(tz = "UTC")
signal_length_days <- 550

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

## 12h & 24h period: resting + DVM ####

signal_12_24 <- dplyr::tibble(t = seq(0,signal_length_days - f, f),
                           depth_m = (5 * (sin((2 * pi)*t) - 1)) + (5 * (sin((4 * pi)*t) - 1)),
                           date_time = seq(signal_start_date,
                                           signal_start_date + lubridate::days(signal_length_days),
                                           length.out = base::length(t)))

## combine generic signals ####

generic_signals <- signal_12 %>%
  dplyr::rename(depth_m_12 = depth_m) %>%
  dplyr::left_join(signal_24, by = c("t", "date_time")) %>%
  dplyr::rename(depth_m_24 = depth_m) %>%
  dplyr::left_join(signal_12_24, by = c("t", "date_time")) %>%
  dplyr::rename(depth_m_12_24 = depth_m) 

rm(signal_12, signal_24, signal_12_24)

## plot signals ####
generic_signals_plot_subset <- ggplot2::ggplot(data = generic_signals %>% dplyr::filter(date_time %>% between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(10), 
                                                                                                              (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(10)))) +
  geom_line(aes(x = date_time, y = depth_m_12), size = .75, colour = 'darkorange') + #, colour = "#357984"
  geom_line(aes(x = date_time, y = depth_m_24), size = .75, colour = 'darkgreen') +
  geom_line(aes(x = date_time, y = depth_m_12_24), size = .75, colour = 'lightblue') +
  labs(x = "", y = "Depth in m") +
  scale_x_datetime(
    date_minor_breaks = "1 day",
    date_breaks = "1 week",
    date_labels = "%b %d", #  
    expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

generic_signals_plot_subset 

# test: XWT of 12 and 24h signals

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

xwt_12_24_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                     values1 = generic_signals %>% dplyr::select(depth_m_12),
                                                     values2 = generic_signals %>% dplyr::select(depth_m_24),
                                                     dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = generic_signals %>% dplyr::select(date_time),
                  dt = 15 * dt_hours, signif_level = 0.95) %>%
  ggplot_wavelet(max_period = max(.$period)) # TODO: understand why argument for max_period is needed

xwt_12_24_plot #%>% View()

gridExtra::grid.arrange(cwt_12_plot, cwt_24_plot, xwt_12_24_plot, ncol = 1)

cwt_12_24_plot <- compute_CWT(values = generic_signals %>% dplyr::select(depth_m_12_24),
                           dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = generic_signals %>% dplyr::select(date_time),
                  dt = 15 * dt_hours) %>%
  ggplot_wavelet(max_period = max(.$period))

cwt_12_24_plot

xwt_12_12_24_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                     values1 = generic_signals %>% dplyr::select(depth_m_12),
                                                     values2 = generic_signals %>% dplyr::select(depth_m_12_24),
                                                     dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = generic_signals %>% dplyr::select(date_time),
                  dt = 15 * dt_hours) %>%
  ggplot_wavelet(max_period = max(.$period)) # TODO: understand why argument for max_period is needed

xwt_12_12_24_plot #%>% View()


# Mustelus asterias depthlogs ####

## load depthlogs ####
depthlog_female <- load_data(filestring = "depthlog_female", folder = base::paste0(base::getwd(), "/examples/data/")) %>%
  # dplyr::mutate(depth_m_sgolayfilt = depth_m %>% signal::sgolayfilt(p = 3, n = 9))
  dplyr::mutate(depth_m_raw = depth_m,
                depth_m = base::scale(depth_m))

depthlog_male <- load_data(filestring = "depthlog_male", folder = base::paste0(base::getwd(), "/examples/data/")) %>%
  # dplyr::mutate(depth_m_sgolayfilt = depth_m %>% signal::sgolayfilt(p = 3, n = 9))
  dplyr::mutate(depth_m_raw = depth_m,
                depth_m = base::scale(depth_m))

# plot depthlogs ####

depthlogs_plot <- ggplot2::ggplot() +
  geom_line(data = depthlog_female %>% 
              dplyr::filter(date_time %>% dplyr::between("2018-09-05 01:00:00" %>% base::as.POSIXct(tz = "UTC"),
                                                         "2018-11-28 01:00:00" %>% base::as.POSIXct(tz = "UTC"))),
            aes(x = date_time, y = -depth_m), size = .75, colour = 'darkorange') +
  # geom_line(data = depthlog_female %>% 
  #             dplyr::filter(date_time %>% dplyr::between("2018-09-05 01:00:00" %>% base::as.POSIXct(tz = "UTC"),
  #                                                        "2018-10-28 01:00:00" %>% base::as.POSIXct(tz = "UTC"))),
  #           aes(x = date_time, y = -depth_m_sgolayfilt), size = .75, colour = 'darkgreen') +
  geom_line(data = depthlog_male %>%
              dplyr::filter(date_time %>% dplyr::between("2018-09-05 01:00:00" %>% base::as.POSIXct(tz = "UTC"),
                                                         "2018-11-28 01:00:00" %>% base::as.POSIXct(tz = "UTC"))),
            aes(x = date_time, y = -depth_m), size = .75, colour = 'darkgreen') +
  labs(x = "", y = "Depth in m") +
  scale_x_datetime(
    date_minor_breaks = "1 day",
    date_breaks = "1 week",
    date_labels = "%b %d", #  
    expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

depthlogs_plot # when both depths are scaled they are different with respect to each other

# XWT female shark ####

## CWT 
cwt_female_plot <- compute_CWT(values = depthlog_female %>% dplyr::select(depth_m),
                               dt = 15 * dt_hours,
                               factor_smallest_scale = 8) %>%
  make_wavelet_df(date_times = depthlog_female %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.9999) %>%
  ggplot_wavelet(max_period = 150,
                 opacity = 0.8)

cwt_female_plot 

# save plot
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/cwt_female_plot.png"), plot = cwt_female_plot)


# XWT
xwt_female_12_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                         values1 = depthlog_female %>% dplyr::select(depth_m), #_sgolayfilt
                                                         values2 = generic_signals  %>%
                                                           dplyr::filter(date_time %>% 
                                                                           dplyr::between(min(depthlog_female$date_time),
                                                                                         max(depthlog_female$date_time))) %>% 
                                                           dplyr::select(depth_m_12),
                                                         dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = depthlog_female %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.95) %>%
  ggplot_wavelet(max_period = 70) # TODO: understand why argument for max_period is needed

xwt_female_12_plot

# save plot
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/xwt_female_12_plot.png"), plot = xwt_female_12_plot)


xwt_female_24_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                         values1 = depthlog_female %>% dplyr::select(depth_m),
                                                         values2 = generic_signals  %>%
                                                           dplyr::filter(date_time %>% 
                                                                           dplyr::between(min(depthlog_female$date_time),
                                                                                          max(depthlog_female$date_time))) %>% 
                                                           dplyr::select(depth_m_24),
                                                         dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = depthlog_female %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.95) %>%
  ggplot_wavelet(max_period = 70) # TODO: understand why argument for max_period is needed

xwt_female_24_plot

# save plot
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/xwt_female_24_plot.png"), plot = xwt_female_24_plot)


xwt_female_12_24_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                          values1 = depthlog_female %>% dplyr::select(depth_m),
                                                          values2 = generic_signals  %>%
                                                            dplyr::filter(date_time %>% 
                                                                            dplyr::between(min(depthlog_female$date_time),
                                                                                           max(depthlog_female$date_time))) %>% 
                                                            dplyr::select(depth_m_12_24),
                                                          dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = depthlog_female %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.95) %>%
  ggplot_wavelet(max_period = 70) # TODO: understand why argument for max_period is needed

xwt_female_12_24_plot

# save plot
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/xwt_female_12_24_plot.png"), plot = xwt_female_12_24_plot)


# XWT male shark ####

# CWT

cwt_male_plot <- compute_CWT(values = depthlog_male %>% dplyr::select(depth_m),
                               dt = 15 * dt_hours,
                               factor_smallest_scale = 8) %>%
  make_wavelet_df(date_times = depthlog_male %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.75) %>%
  ggplot_wavelet(max_period = 70)

cwt_male_plot 

ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/cwt_male_plot.png"), plot = cwt_male_plot)

# XWT

xwt_male_12_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                         values1 = depthlog_male %>% dplyr::select(depth_m),
                                                         values2 = generic_signals  %>%
                                                           dplyr::filter(date_time %>% 
                                                                           dplyr::between(min(depthlog_male$date_time),
                                                                                          max(depthlog_male$date_time))) %>% 
                                                           dplyr::select(depth_m_12),
                                                         dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = depthlog_male %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.95) %>%
  ggplot_wavelet(max_period = 70) # TODO: understand why argument for max_period is needed

xwt_male_12_plot

# save plot
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/xwt_male_12_plot.png"), plot = xwt_male_12_plot)


xwt_male_24_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                         values1 = depthlog_male %>% dplyr::select(depth_m),
                                                         values2 = generic_signals  %>%
                                                           dplyr::filter(date_time %>% 
                                                                           dplyr::between(min(depthlog_male$date_time),
                                                                                          max(depthlog_male$date_time))) %>% 
                                                           dplyr::select(depth_m_24),
                                                         dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = depthlog_male %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.95) %>%
  ggplot_wavelet(max_period = 70) # TODO: understand why argument for max_period is needed

xwt_male_24_plot

# save plot
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/xwt_male_24_plot.png"), plot = xwt_male_24_plot)


xwt_male_12_24_plot <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                       values1 = depthlog_male %>%
                                                         dplyr::filter(lubridate::date(date_time) %>% 
                                                                         dplyr::between(min(lubridate::date(depthlog_female$date_time)),
                                                                                        max(lubridate::date(depthlog_female$date_time)))) %>% 
                                                         dplyr::select(depth_m) ,
                                                       values2 = generic_signals  %>%
                                                         dplyr::filter(lubridate::date(date_time) %>% 
                                                                         dplyr::between(min(lubridate::date(depthlog_female$date_time)),
                                                                                        max(lubridate::date(depthlog_female$date_time)))) %>% 
                                                         dplyr::select(depth_m_12_24),
                                                       dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = depthlog_male%>%
                    dplyr::filter(lubridate::date(date_time) %>% 
                                    dplyr::between(min(lubridate::date(depthlog_female$date_time)),
                                                   max(lubridate::date(depthlog_female$date_time)))) %>% 
                    dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.95) %>%
  ggplot_wavelet(max_period = 70) # TODO: understand why argument for max_period is needed

xwt_male_12_24_plot

# save plot
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/xwt_male_12_24_plot.png"), plot = xwt_male_12_24_plot)


# gridExtra::grid.arrange(xwt_female_12_plot, xwt_male_12_plot,
#                         xwt_female_24_plot, xwt_male_24_plot)


# comparison XWT female & male ####

gridExtra::grid.arrange(xwt_female_12_24_plot + labs(title = 'female'), 
                        xwt_male_12_24_plot + labs(title = 'male'),
                        ncol = 2)

xwt_female_12_24_df <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                          values1 = depthlog_female %>% dplyr::select(depth_m),
                                                          values2 = generic_signals  %>%
                                                            dplyr::filter(date_time %>% 
                                                                            dplyr::between(min(depthlog_female$date_time),
                                                                                           max(depthlog_female$date_time))) %>% 
                                                            dplyr::select(depth_m_12_24),
                                                          dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = depthlog_female %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.95) %>%
  dplyr::filter(period <= 3000)

xwt_male_12_24_df <- compute_bivariate_wavelet_analysis(type = "cross wavelet",
                                                        values1 = depthlog_male %>% dplyr::select(depth_m),
                                                        values2 = generic_signals  %>%
                                                          dplyr::filter(date_time %>% 
                                                                          dplyr::between(min(depthlog_male$date_time),
                                                                                         max(depthlog_male$date_time))) %>% 
                                                          dplyr::select(depth_m_12_24),
                                                        dt = 15 * dt_hours) %>%
  make_wavelet_df(date_times = depthlog_male %>% dplyr::select(date_time),
                  dt = 15 * dt_hours,
                  signif_level = 0.95) %>%
  dplyr::filter(period <= 3000)

xwt_both_12_24_df <- xwt_male_12_24_df %>%
  dplyr::filter(date %>% 
                  dplyr::between(min(xwt_female_12_24_df$date),
                                 max(xwt_female_12_24_df$date))) %>%
  dplyr::mutate(power_log = power_log + (xwt_female_12_24_df$power_log),
                significance = significance + (xwt_female_12_24_df$significance),
                sig = sig * (xwt_female_12_24_df$sig))
                # sig = ifelse(significance >= stats::quantile(signif$significance, probs = 0.95), 1, 0))
  # dplyr::mutate(power_log = power_log + (xwt_female_12_24_df %>% dplyr::select(power_log)),
  #               significance = significance + (xwt_female_12_24_df %>% dplyr::select(significance)))
  
xwt_both_12_24_plot <- ggplot_wavelet(wavelet_df = xwt_both_12_24_df, max_period = 70)             
xwt_both_12_24_plot

# save plot
ggplot2::ggsave(filename = base::paste0(base::getwd(), "/examples/plots/xwt_both_12_24_plot.png"), plot = xwt_both_12_24_plot)

gridExtra::grid.arrange(xwt_female_12_24_plot,
                        xwt_male_12_24_plot,
                        xwt_both_12_24_plot)

# To do morgen:
# ploth mit both in 3 kategorien einteilen: beide, nur female, nur male,
# und das ist dann das endergebnis. dann alle plots schön abspeichern und in ne Präsi hauen
# und rausfinden wie ich die random funktion generiere, evtl. mit for-loops


