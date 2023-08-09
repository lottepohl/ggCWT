# script to make example figures for FFT and wavelet transform for thesis defense

rm(list = ls())

# WORKSPACE ####
library(biwavelet)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(sf)
library(ggspatial)
library(marmap)
library(ggrepel)
library(lubridate)
library(gridExtra)

# 
# library("rnaturalearth")
# library("rnaturalearthdata")

# dir_path <- "C:/Users/lotte.pohl/Documents/github_repos/MasterThesis_LottePohl"
dir_path <- "C:/users/lotte/Documents/Repos_Gitkraken/MasterThesis_LottePohl_MDA"
# path_envdata <- paste0(dir_path, "/00_data/environmental_layers/")
# path_boundaries <- paste0(dir_path, "/00_data/marine_boundaries/")
# path_maps <- "C:/Users/lotte/HiDrive/Master/MasterThesis/presentations/20230628_30_Thesis_defense/plots_maps/"

# paste0(dir_path, "/01_code/02_load_data/manuscript_figures/load_plots.R") %>% base::source()

paste0(dir_path, "/01_code/06_functions/functions.R") %>% source()
paste0(dir_path, "/01_code/04_analyses/FFT/calculate_fft_psd.R") %>% source()
# source(paste0(dir_path, "/01_code/02_load_data/load_environmental_data.R"))
# source(paste0(dir_path, "/01_code/02_load_data/load_human_activities.R"))
# source(paste0(dir_path, "/01_code/02_load_data/load_marine_boundaries.R"))
# source(paste0(dir_path, "/01_code/02_load_data/load_acoustic_detections.R"))
# source(paste0(dir_path, "/01_code/02_load_data/load_bathy.R"))
# paste0(dir_path, "/01_code/02_load_data/load_dst_geolocation_output.R") %>% base::source()
# paste0(dir_path, "/01_code/02_load_data/manuscript_figures/load_tables.R") %>% base::source()

## load data ####
paste0(dir_path, "/01_code/02_load_data/load_dst_summarystatistics.R") %>% base::source()
# paste0(dir_path, "/01_code/02_load_data/load_acoustic_detections.R") %>% base::source()
# to do: choose df's to load to reduce workspace size
paste0(dir_path, "/01_code/02_load_data/load_wavelet_results.R") %>% base::source()
# paste0(dir_path, "/01_code/02_load_data/load_autocorrelation_results.R") %>% base::source()
paste0(dir_path, "/01_code/02_load_data/load_depth_temp_logs.R") %>% base::source()
paste0(dir_path, "/01_code/02_load_data/load_fft_results.R") %>% base::source()
# paste0(dir_path, "/01_code/02_load_data/load_cpd_results.R") %>% base::source()
# paste0(dir_path, "/01_code/02_load_data/load_vertical_space_use_analysis.R") %>% base::source()
# paste0(dir_path, "/01_code/02_load_data/manuscript_figures/load_tables.R") %>% base::source()
# paste0(dir_path, "/01_code/02_load_data/manuscript_figures/load_models.R") %>% base::source()

## set path were all figures are saved ####
# plot_path <- paste0(dir_path, "/01_code/00_thesis_manuscript/figures/")

## theme ####

thesis_theme <- ggplot2::theme(
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
ggplot2::theme_set(thesis_theme)


# FFT ####

plot_fft <- function(fft_result, tag_serial_number_short, period_upperlim = 40, period_lowerlim = 0.05){
  periodogram <- ggplot(data = fft_result %>% dplyr::filter(period < period_upperlim & period > period_lowerlim)) + 
    geom_line(aes(x = period, y = spec), colour = "black") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(2, period_upperlim, by = 2)) +
    labs(y = "Spectral density", x = "Period in hours") #, title = paste0("tag ", tag_serial_number_short)
  return(periodogram)
}

calc_fft <- function(depth_log, sample_freq){
  # `depth_log` contains the columns `depth_m` and `date_time`.
  
  # prepare/format data
  depth_log <- depth_log %>% 
    # make time vector
    mutate(t = difftime(date_time, date_time[1], units = "hours") %>% as.numeric()) %>% 
    #dplyr::filter out NAs
    dplyr::filter(!is.na(depth_m))
  
  # define tmax and n
  tmax <- depth_log$t %>% max()
  n <- depth_log$t %>% length()
  
  # mit nullen auffuellen bis naechste n^2
  next_power <- 2^ceiling(log2(n))
  depth_log_padded <- c(depth_log$depth_m, rep(0, next_power - n))
  
  # compute the frequency spectrum
  spec <- abs(fft(depth_log_padded))^2
  
  # compute the frequency vector
  # freq <- seq(0, length(spec) - 1) * (1/sample_freq) / length(spec)
  freq <- seq(0, length(spec) - 1) * (sample_freq) / length(spec)
  
  # make result dataframe
  result_fft <- cbind(spec, freq) %>% as.data.frame() %>% 
    # calculate period
    mutate(period = 1 / freq)
  
  return(result_fft)
}

# plot periodogram 
plot_periodogram <- function(fft_result, tag_serial_number_short, period_upperlim = 40, period_lowerlim = 0.05, path = plot_path){
  # todo: get local rule or set values for period upper and lower lim
  periodogram <- ggplot(data = fft_result %>% dplyr::filter(period < period_upperlim & period > period_lowerlim)) + 
    geom_line(aes(x = period, y = spec), colour = "#357984", size = 1) + #theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    # theme_bw(base_size = 12) +
    scale_x_continuous(expand = c(0,0), breaks = seq(2, period_upperlim, by = 2)) +
    labs(y = "Spectral density", x = "Period in hours") #, title = paste0("tag ", tag_serial_number_short)
  
  # ggplot2::ggsave(filename = paste0(plot_path, "periodogram_", tag_serial_number_short, ".pdf"), plot = periodogram, width = 18, height = 12, units = "cm")
  # ggplot2::ggsave(filename = paste0(plot_path, "periodogram_", tag_serial_number_short, ".png"), plot = periodogram, width = 18, height = 12, units = "cm")
  return(periodogram)
}

fft_calc_plot <- function(depth_log, tag_serial_num_short, sample_frequency){
  fft_res <- calc_fft(depth_log = depth_log %>% dplyr::filter(tag_serial_number == paste0("1293", tag_serial_num_short)),
                      sample_freq = sample_frequency)
  
  pgram <- plot_periodogram(fft_result = fft_res, 
                            tag_serial_number_short = tag_serial_num_short)
  pgram %>% return()
}


# CWT ####


compute_wavelet_hrperiod <- function(depthlog = masterias_depth_temp, tag_serial_num = "1293308", parameter_name = "depth_m", dt = dt, factor_smallest_scale = 8){
  # downsample raw depthlog
  depthlog_downsampled <- depthlog %>% 
    ungroup() %>%
    dplyr::filter(tag_serial_number == tag_serial_num,
                  row_number() %% ((dt * 60) / 2) == 0) %>%
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

wavelet_output_compare_hrperiod <- function(depthlog = masterias_depth_temp, tag_serial_num = "1293308", parameter_name = "depth_m", dt = dt, wt_output = wt_output){
  
  # downsample raw depthlog
  depthlog_downsampled <- depthlog %>% 
    ungroup() %>%
    dplyr::filter(tag_serial_number == tag_serial_num,
                  row_number() %% ((dt * 60) / 2) == 0) %>%
    dplyr::select(c(parameter_name, date_time))# explanation: ... %% dt[hour] * 60[min] / 2[min] <- sample interval (every 2 min), e.g. for dt = 1: ...%% 30 (i.e. get every 30th val)
  
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


plot_wavelet_hrperiod_example <- function(wt_df = wt_df, type = c("power", "significance", "power_log"),
                                  date = TRUE, max_period = 72 #hours = 4 days
){
  # transformation function for the y axis
  my_trans <- scales::trans_new("log2_reverse", function(x) -log2(x), function(x) 2^-x)
  
  # wt_df <- wt_df %>% dplyr::filter(date > (min(wt_df$date) + lubridate::days(7))) # cut first week of data off to avoid looking at tagging effect
  
  # y axis labels
  y_breaks <- 2^floor(log2(wt_df$period)) %>% unique()
  y_breaks <- y_breaks[y_breaks <= max_period]
  # transform dates
  # wt_df$date_time <- wt_df$date_time %>% as.POSIXct(tz = "UTC") # not needed bc we have the posict objects
  
  # # x axis labels
  # ifelse(date %>% isTRUE(),
  #        x_breaks <- c(wt_df$date_time[1], wt_df$date_time[(1/5) * n_data], wt_df$date_time[(2/5) * n_data], wt_df$date_time[(3/5) * n_data],
  #                      wt_df$date_time[(4/5) * n_data], wt_df$date_time[(5/5) * n_data])
  # ,
  # x_breaks <- sprintf("%03d", seq(from = 0, to = n_data, by = 100)))
  # 
  #plot
  
  wt_df <- wt_df  %>% dplyr::filter(period <= max_period)
  
  # change max and min date to include max x axis label completely
  # max_date <- max(wt_df$date) + lubridate::days(10)
  # min_date <- min(wt_df$date) -lubridate::days(10)
  
  
  ifelse(date %>% isTRUE(),
         
         ifelse(type == "power_log",
                
                plot <- ggplot(data = wt_df) +
                  geom_tile(aes(x = date, y = period, fill = power_log),
                            position = "identity",
                            alpha = 0.7) + #0.5
                  # geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date, y = period, fill = power_log),
                  #           position = "identity") +
                  scale_y_continuous(trans = my_trans,
                                     breaks = y_breaks, 
                                     expand = c(0,0)) +
                  # scale_x_discrete(breaks = x_breaks) +
                  # scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                  scale_x_datetime(
                    date_minor_breaks = "1 day",
                    date_breaks = "1 week",
                    date_labels = "%b %d",
                    expand = c(0,0)
                    # ,limits = c(min_date, max_date)
                    ) +
                  scale_fill_viridis_c(direction = 1, option = "turbo") +
                  labs(x = "", y = "Period in hours", fill = "log2(Power)") +
                  theme(legend.position = "bottom", # "bottom",
                        legend.box = "horizontal",
                        legend.margin = margin(t = -15)) #+
                  # theme(axis.text.x = element_text(angle = 15, hjust = 0.5))
                ,
                
                ifelse(type == "significance",
                       
                       plot <- ggplot(data = wt_df) +
                         geom_tile(aes(x = date_posicxt, y = period, fill = significance),
                                   position = "identity",
                                   alpha = 0.65) +
                         geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date_posicxt, y = period, fill = significance),
                                   position = "identity") +
                         scale_y_continuous(trans = my_trans,
                                            breaks = y_breaks, 
                                            expand = c(0,0)) +
                         # scale_x_discrete(breaks = x_breaks) +
                         scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                         scale_fill_viridis_c(direction = 1, option = "turbo") +
                         labs(x = "date", y = "period in hours", fill = "significance") #+
                       # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                       ,
                       
                       plot <- ggplot(data = wt_df) +
                         geom_tile(aes(x = date_posicxt, y = period, fill = power),
                                   position = "identity",
                                   alpha = 0.65) +
                         geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date_posicxt, y = period, fill = power),
                                   position = "identity") +
                         scale_y_continuous(trans = my_trans,
                                            breaks = y_breaks, 
                                            expand = c(0,0)) +
                         # scale_x_discrete(breaks = x_breaks) +
                         scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                         scale_fill_viridis_c(direction = 1, option = "turbo") +
                         labs(x = "date", y = "period in hours", fill = "power") #+
                       # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                       
                )
         )
         ,
         
         ifelse(type == "power_log",
                
                plot <- ggplot(data = wt_df) +
                  geom_tile(aes(x = t, y = period, fill = power_log),
                            position = "identity",
                            alpha = 0.65) +
                  geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = power_log),
                            position = "identity") +
                  scale_y_continuous(trans = my_trans,
                                     breaks = y_breaks, 
                                     expand = c(0,0)) +
                  # scale_x_discrete(breaks = x_breaks) +
                  scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                  scale_fill_viridis_c(direction = 1, option = "turbo") +
                  labs(x = "date", y = "period in hours", fill = "log2(power)") #+
                # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                ,
                
                ifelse(type == "significance",
                       
                       plot <- ggplot(data = wt_df) +
                         geom_tile(aes(x = t, y = period, fill = significance),
                                   position = "identity",
                                   alpha = 0.65) +
                         geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = significance),
                                   position = "identity") +
                         scale_y_continuous(trans = my_trans,
                                            breaks = y_breaks, 
                                            expand = c(0,0)) +
                         # scale_x_discrete(breaks = x_breaks) +
                         scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                         scale_fill_viridis_c(direction = 1, option = "turbo") +
                         labs(x = "date", y = "period in hours", fill = "significance") #+
                       # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                       ,
                       
                       plot <- ggplot(data = wt_df) +
                         geom_tile(aes(x = t, y = period, fill = power),
                                   position = "identity",
                                   alpha = 0.65) +
                         geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = power),
                                   position = "identity") +
                         scale_y_continuous(trans = my_trans,
                                            breaks = y_breaks, 
                                            expand = c(0,0)) +
                         # scale_x_discrete(breaks = x_breaks) +
                         scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                         scale_fill_viridis_c(direction = 1, option = "turbo") +
                         labs(x = "date", y = "period in hours", fill = "power") #+
                       # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                )
         )
  )
  return(plot)
}

plot_wavelet_hrperiod <- function(wt_df = wt_df, type = c("power", "significance", "power_log"),
                                  date = TRUE, max_period = 72 #hours = 4 days
){
  # transformation function for the y axis
  my_trans <- scales::trans_new("log2_reverse", function(x) -log2(x), function(x) 2^-x)
  
  wt_df <- wt_df %>% dplyr::filter(date > (min(wt_df$date) + lubridate::days(7))) # cut first week of data off to avoid looking at tagging effect
  
  # y axis labels
  y_breaks <- 2^floor(log2(wt_df$period)) %>% unique()
  y_breaks <- y_breaks[y_breaks <= max_period]
  # transform dates
  # wt_df$date_time <- wt_df$date_time %>% as.POSIXct(tz = "UTC") # not needed bc we have the posict objects
  
  # # x axis labels
  # ifelse(date %>% isTRUE(),
  #        x_breaks <- c(wt_df$date_time[1], wt_df$date_time[(1/5) * n_data], wt_df$date_time[(2/5) * n_data], wt_df$date_time[(3/5) * n_data],
  #                      wt_df$date_time[(4/5) * n_data], wt_df$date_time[(5/5) * n_data])
  # ,
  # x_breaks <- sprintf("%03d", seq(from = 0, to = n_data, by = 100)))
  # 
  #plot
  
  wt_df <- wt_df  %>% dplyr::filter(period <= max_period)
  
  # change max and min date to include max x axis label completely
  max_date <- max(wt_df$date) + lubridate::days(10)
  min_date <- min(wt_df$date) -lubridate::days(10)
  
  
  ifelse(date %>% isTRUE(),
         
         ifelse(type == "power_log",
                
                plot <- ggplot(data = wt_df) +
                  geom_tile(aes(x = date, y = period, fill = power_log),
                            position = "identity",
                            alpha = 0.7) + #0.5
                  # geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date, y = period, fill = power_log),
                  #           position = "identity") +
                  scale_y_continuous(trans = my_trans,
                                     breaks = y_breaks, 
                                     expand = c(0,0)) +
                  # scale_x_discrete(breaks = x_breaks) +
                  # scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                  scale_x_datetime(
                    # date_minor_breaks = "1 month",
                    date_breaks = "1 month",
                    date_labels = "%b '%y",
                    expand = c(0,0),
                    limits = c(min_date, max_date)) +
                  scale_fill_viridis_c(direction = 1, option = "turbo") +
                  labs(x = "", y = "Period in hours", fill = "log2(Power)") +
                  theme(legend.position = "bottom", # "bottom",
                        legend.box = "horizontal",
                        legend.margin = margin(t = -15)) +
                  theme(axis.text.x = element_text(angle = 15, hjust = 0.5))
                ,
                
                ifelse(type == "significance",
                       
                       plot <- ggplot(data = wt_df) +
                         geom_tile(aes(x = date_posicxt, y = period, fill = significance),
                                   position = "identity",
                                   alpha = 0.65) +
                         geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date_posicxt, y = period, fill = significance),
                                   position = "identity") +
                         scale_y_continuous(trans = my_trans,
                                            breaks = y_breaks, 
                                            expand = c(0,0)) +
                         # scale_x_discrete(breaks = x_breaks) +
                         scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                         scale_fill_viridis_c(direction = 1, option = "turbo") +
                         labs(x = "date", y = "period in hours", fill = "significance") #+
                       # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                       ,
                       
                       plot <- ggplot(data = wt_df) +
                         geom_tile(aes(x = date_posicxt, y = period, fill = power),
                                   position = "identity",
                                   alpha = 0.65) +
                         geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date_posicxt, y = period, fill = power),
                                   position = "identity") +
                         scale_y_continuous(trans = my_trans,
                                            breaks = y_breaks, 
                                            expand = c(0,0)) +
                         # scale_x_discrete(breaks = x_breaks) +
                         scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                         scale_fill_viridis_c(direction = 1, option = "turbo") +
                         labs(x = "date", y = "period in hours", fill = "power") #+
                       # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                       
                )
         )
         ,
         
         ifelse(type == "power_log",
                
                plot <- ggplot(data = wt_df) +
                  geom_tile(aes(x = t, y = period, fill = power_log),
                            position = "identity",
                            alpha = 0.65) +
                  geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = power_log),
                            position = "identity") +
                  scale_y_continuous(trans = my_trans,
                                     breaks = y_breaks, 
                                     expand = c(0,0)) +
                  # scale_x_discrete(breaks = x_breaks) +
                  scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                  scale_fill_viridis_c(direction = 1, option = "turbo") +
                  labs(x = "date", y = "period in hours", fill = "log2(power)") #+
                # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                ,
                
                ifelse(type == "significance",
                       
                       plot <- ggplot(data = wt_df) +
                         geom_tile(aes(x = t, y = period, fill = significance),
                                   position = "identity",
                                   alpha = 0.65) +
                         geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = significance),
                                   position = "identity") +
                         scale_y_continuous(trans = my_trans,
                                            breaks = y_breaks, 
                                            expand = c(0,0)) +
                         # scale_x_discrete(breaks = x_breaks) +
                         scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                         scale_fill_viridis_c(direction = 1, option = "turbo") +
                         labs(x = "date", y = "period in hours", fill = "significance") #+
                       # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                       ,
                       
                       plot <- ggplot(data = wt_df) +
                         geom_tile(aes(x = t, y = period, fill = power),
                                   position = "identity",
                                   alpha = 0.65) +
                         geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = power),
                                   position = "identity") +
                         scale_y_continuous(trans = my_trans,
                                            breaks = y_breaks, 
                                            expand = c(0,0)) +
                         # scale_x_discrete(breaks = x_breaks) +
                         scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                         scale_fill_viridis_c(direction = 1, option = "turbo") +
                         labs(x = "date", y = "period in hours", fill = "power") #+
                       # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
                )
         )
  )
  return(plot)
}

# DEPTHLOG ####

plot_raw_and_sum_depth <- function(depthlog = masterias_depth_temp, summary_depthlog = long_dst_date, tag_serial_num = "1293308",
                           start_date_summer_chr, end_date_summer_chr, ymin_summer, ymax_summer, colour_summer,
                           start_date_winter_chr, end_date_winter_chr, ymin_winter, ymax_winter, colour_winter){
  data_raw <- depthlog %>% dplyr::filter(tag_serial_number == tag_serial_num, 
                                     row_number() %% 15 == 0)
  # , date > (min(wt_df$date) + lubridate::days(7))) # cut first week of data off to avoid looking at tagging effect
  
  data_sum <- summary_depthlog %>% dplyr::filter(tag_serial_number == tag_serial_num)
  
  max_date <- (max(data_raw$date_time %>% lubridate::date()) + lubridate::days(10)) %>% as.POSIXct()
  min_date <- (min(data_raw$date_time %>% lubridate::date()) -lubridate::days(10)) %>% as.POSIXct()
  
  plot <- ggplot2::ggplot() + 
    # annotate(geom = "rect", xmin = start_date_summer_chr %>% as.POSIXct(), xmax = end_date_summer_chr %>% as.POSIXct(),
    #          ymin = ymin_summer, ymax = ymax_summer, fill = colour_summer, alpha = 0.5) +
    # annotate(geom = "rect", xmin = start_date_winter_chr %>% as.POSIXct(), xmax = end_date_winter_chr %>% as.POSIXct(),
    #          ymin = ymin_winter, ymax = ymax_winter, fill = colour_winter, alpha = 0.5) +
    # depth line
    geom_line(data = data_raw, aes(x = date_time, y = -depth_m), size = 0.15, colour = "#478B96") + # darker: #003E4C, lighter: #478B96, mittel: #357984
    geom_line(data = data_sum,  aes(x = date, y = -depth_median_sgolay), size = 0.75, colour = "black") +
    # panel outline
    # geom_rect(xmin = start_date_summer_chr %>% as.POSIXct(), xmax = end_date_summer_chr %>% as.POSIXct(),
    #           ymin = ymin_summer, ymax = ymax_summer, colour = colour_summer, size = 0.25, fill = "transparent") +
    # geom_rect(xmin = start_date_winter_chr %>% as.POSIXct(), xmax = end_date_winter_chr %>% as.POSIXct(),
    #           ymin = ymin_winter, ymax = ymax_winter, colour = colour_winter, size = 2, fill = "transparent") +
    # settings
    scale_x_datetime(date_breaks = "1 month",
                     date_labels = "%b '%y",
                     expand = c(0,0),
                     limits = c(min_date, max_date)) +
    theme(axis.text.x = element_text(angle = 15, hjust = 0.5)) +
    scale_y_continuous(expand = c(0,0), limits = c(-(max(data_raw$depth_m)) -2,0)) +
    labs(x = "", y = "Depth in m") +
    theme(legend.position = "none",
          legend.box = "horizontal", legend.margin = margin(t = -15))
  
  plot
}

plot_depth_subset_summer <- function(depthlog = masterias_depth_temp, tag_serial_num, start_date_chr, end_date_chr){
  data <- depthlog %>% dplyr::filter(tag_serial_number == tag_serial_num,
                                     # lubridate::date(date_time) %>% between(as.POSIXct(start_date_chr), as.POSIXct(end_date_chr)),
                                     date_time %>% between(as.POSIXct(start_date_chr), as.POSIXct(end_date_chr)),
                                     row_number() %% 1 == 0)
  plot <- ggplot2::ggplot(data = data, mapping = aes(x = date_time, y = -depth_m)) + 
    # geom_point(size = 1.5) + # aes(colour = lubridate::hour(date_time)),%>% as.factor()
    geom_line(size = 0.15, colour = "#00718A") +
    scale_x_datetime(date_breaks = "1 day",
                     date_minor_breaks = "12 hours",
                     date_labels = "%d.%m.%y", #'%y
                     expand = c(0,0)
                     # ,limits = c(min_date, max_date)
    ) +
    # theme(axis.text.x = element_text(angle = 15, hjust = 0.25)) +
    scale_y_continuous(expand = c(0,0), limits = c(-21,0)) +
    theme(legend.position = "bottom", legend.box = "horizontal") +
    labs(x = "", y = "Depth in m", color = "Hour of the day") #+ #title = paste0("tag ", tag_serial_number_short), 
  
}

plot_depth_subset_winter <- function(depthlog = masterias_depth_temp, tag_serial_num, start_date_chr, end_date_chr){
  data <- depthlog %>% dplyr::filter(tag_serial_number == tag_serial_num,
                                     # lubridate::date(date_time) %>% between(as.POSIXct(start_date_chr), as.POSIXct(end_date_chr)),
                                     date_time %>% between(as.POSIXct(start_date_chr), as.POSIXct(end_date_chr)),
                                     row_number() %% 2 == 0)
  
  # dplyr::select(date_time) %>%
  # lubridate::date() %>% pull() %>% max()
  
  plot <- ggplot2::ggplot(data = data, mapping = aes(x = date_time, y = -depth_m)) + 
    # geom_point(size = 1.5) + # aes(colour = lubridate::hour(date_time)),%>% as.factor()
    geom_line(size = 0.3, colour = "#00718A") +
    scale_x_datetime(date_breaks = "1 week",
                     date_minor_breaks = "1 day",
                     date_labels = "%d.%m.%y", #'%y
                     expand = c(0,0)
                     # ,limits = c(min_date, max_date)
    ) +
    # theme(axis.text.x = element_text(angle = 15, hjust = 0.25)) +
    scale_y_continuous(expand = c(0,0), limits = c(-round(max(data$depth_m, na.rm = T)),0)) + #round to nearest multiple of 10 ceiling(max(data$depth_m)
    theme(legend.position = "bottom", legend.box = "horizontal") +
    labs(x = "", y = "Depth in m", color = "Hour of the day") #+ #title = paste0("tag ", tag_serial_number_short), 
  
}


# illustration T and f ####

signal_T_f <- tibble(t = seq(0, 100, 0.1),
                 depth_m =  sin((0.045 * pi)*t)) 

## plot signal ####
signal_T_f_plot <- ggplot() +
  geom_line(data = signal_T_f, aes(x = t, y = depth_m), size = .75, colour = "#357984" ) +
  labs(x = "time step", y = "y") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

signal_T_f_plot

## 12 and 24h illustration

signal_12_24 <- tibble(t = seq(0, 100, 0.1),
                     depth_m =  sin((0.15 * pi)*t)) 

## plot signal
signal_12_24_plot <- ggplot() +
  geom_line(data = signal_12_24, aes(x = t, y = depth_m), size = .75, colour = "black" ) +
  labs(x = "time step", y = "y") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

signal_12_24_plot

# 12 and 24h frequencies ####

f = 1/720 #because if 30 samples per hour are to be taken, these amount to 30 * 24 = 720 samples per day
sample_frequency_hours <- 30

signal <- tibble(t = seq(0,30 - f, f),
                 depth_m = 10 * (sin((4 * pi)*t) - 1),
                 tag_serial_number = "1293308",
                 date_time = seq("2020-01-01 01:00:00" %>% as.POSIXct(tz = "UTC"),
                                 ("2020-01-01 01:00:00" %>% as.POSIXct(tz = "UTC")) + lubridate::days(30),
                                 length.out = length(t))) %>%
  dplyr::mutate(depth_m = ifelse(date_time > ("2020-01-01 01:00:00" %>% as.POSIXct(tz = "UTC")) + lubridate::days(15),
                                 10 * (sin((2 * pi)*t) - 1),
                                 depth_m))


## plot signal ####
signal_plot <- ggplot() +
  geom_line(data = signal, aes(x = date_time, y = depth_m), size = .75, colour = "#357984" ) +
  labs(x = "", y = "Depth in m") +
  scale_x_datetime(
    date_minor_breaks = "1 day",
    date_breaks = "1 week",
    date_labels = "%b %d",
    expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

signal_plot

## calc & plot fft ####

signal_fft_plot <- fft_calc_plot(depth_log = signal, tag_serial_num_short = "308",
                            sample_frequency = sample_frequency_hours)
signal_fft_plot

# CWT ####

dt_signal <- (1/6)
factor_smallest_scale_signal <- 12

signal_cwt_res <- compute_wavelet_hrperiod(depthlog = signal, tag_serial_num = "1293308", parameter_name = "depth_m", 
                                           dt = dt_signal, factor_smallest_scale = factor_smallest_scale_signal)

signal_cwt_input_for_ggplot <- wavelet_output_compare_hrperiod(depthlog = signal, tag_serial_num = "1293308", parameter_name = "depth_m",
                                                               dt = dt_signal, wt_output = signal_cwt_res)

signal_cwt_plot <- plot_wavelet_hrperiod_example(wt_df = signal_cwt_input_for_ggplot, type = "power_log", date = T, max_period = 72)

signal_cwt_plot

### tag 308 ####

p_308_wavelet_depth_hr <- plot_wavelet_hrperiod(wt_df = wt_df_308_depth_hr, type = "power_log",
                                                date = TRUE, max_period = 512)

p_308_wavelet_depth_hr

### tag 321 ####

p_321_wavelet_depth_hr <- plot_wavelet_hrperiod(wt_df = wt_df_321_depth_hr, type = "power_log",
                                                date = TRUE, max_period = 512)
p_321_wavelet_depth_hr


# compare plots ####

gridExtra::grid.arrange(signal_plot, signal_cwt_plot, ncol = 1)


# depthlogs ####

p_308_depth <- plot_raw_and_sum_depth(depthlog = masterias_depth_temp, summary_depthlog = long_dst_date, tag_serial_num = "1293308",
                              start_date_summer_chr = "2018-09-21 14:00:00 UTC", end_date_summer_chr = "2018-09-24 14:00:00 UTC", ymin_summer = -20, ymax_summer = 0, colour_summer = "#F69A41",
                              start_date_winter_chr = "2019-02-08", end_date_winter_chr = "2019-03-07", ymin_winter = -75.52, ymax_winter = 0, colour_winter = "#508893") +
  theme(axis.text = element_text(family = "serif", size = 5.5))

p_308_depth

p_321_depth <- plot_raw_and_sum_depth(depthlog = masterias_depth_temp, summary_depthlog = long_dst_date, tag_serial_num = "1293321",
                                      start_date_summer_chr = "2018-09-21 14:00:00 UTC", end_date_summer_chr = "2018-09-24 14:00:00 UTC", ymin_summer = -20, ymax_summer = 0, colour_summer = "#F69A41",
                                      start_date_winter_chr = "2019-02-08", end_date_winter_chr = "2019-03-07", ymin_winter = -75.52, ymax_winter = 0, colour_winter = "#508893") +
  theme(axis.text = element_text(family = "serif", size = 5.5))

p_321_depth

## summer & winter subsets ####


p_308_depth_summer <- plot_depth_subset_summer(depthlog = masterias_depth_temp, tag_serial_num = "1293308", start_date_chr = "2018-09-21 14:00:00 UTC", end_date_chr = "2018-09-24 14:00:00 UTC") +
  theme(panel.border = element_rect(color = "#D5CF5E", fill = NA, size = 2))
p_308_depth_summer

p_308_depth_winter <- plot_depth_subset_winter(depthlog = masterias_depth_temp, tag_serial_num = "1293308", start_date_chr = "2019-02-08 14:00:00 UTC", end_date_chr = "2019-03-07 14:00:00 UTC") +
  theme(panel.border = element_rect(color = "#88C8D9", fill = NA, size = 2))
p_308_depth_winter

p_308_depth_winter_70m <- plot_depth_subset_winter(depthlog = masterias_depth_temp, tag_serial_num = "1293308", start_date_chr = "2018-10-10 14:00:00 UTC", end_date_chr = "2018-12-01 14:00:00 UTC") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
p_308_depth_winter_70m


p_321_depth_summer <- plot_depth_subset_summer(depthlog = masterias_depth_temp, tag_serial_num = "1293321", start_date_chr = "2018-09-15 14:00:00 UTC", end_date_chr = "2018-09-18 14:00:00 UTC") +
  theme(panel.border = element_rect(color = "#D5CF5E", fill = NA, size = 2))
p_321_depth_summer

# p_321_depth_winter <- plot_depth_subset_winter(depthlog = masterias_depth_temp, tag_serial_num = "1293321", start_date_chr = "2019-02-08", end_date_chr = "2019-02-27") +
p_321_depth_winter <- plot_depth_subset_winter(depthlog = masterias_depth_temp, tag_serial_num = "1293321", start_date_chr = "2018-12-21", end_date_chr = "2019-01-08") +
  theme(panel.border = element_rect(color = "#88C8D9", fill = NA, size = 2))
p_321_depth_winter #%>% ggplotly()

# save plots ####

plot_height <- 4
plot_width <- 15
plot_halfpagewidth <- 9
plot_quarterpagewidth <- 4.5

### T and F signal ####
save_data(data = signal_T_f_plot, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_T_f_plot)), ".png"), plot = signal_T_f_plot, width = plot_width, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_T_f_plot)), ".pdf"), plot = signal_T_f_plot, width = plot_width, height = plot_height, units = "cm")

save_data(data = signal_12_24_plot, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_12_24_plot)), ".png"), plot = signal_12_24_plot, width = plot_width, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_12_24_plot)), ".pdf"), plot = signal_12_24_plot, width = plot_width, height = plot_height, units = "cm")


### signal ####
save_data(data = signal_plot, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_plot)), ".png"), plot = signal_plot, width = plot_width, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_plot)), ".pdf"), plot = signal_plot, width = plot_width, height = plot_height, units = "cm")


### FFT ####
save_data(data = signal_fft_plot, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_fft_plot)), ".png"), plot = signal_fft_plot, width = plot_width, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_fft_plot)), ".pdf"), plot = signal_fft_plot, width = plot_width, height = plot_height, units = "cm")


### CWT ####
save_data(data = signal_cwt_plot, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_cwt_plot)), ".png"), plot = signal_cwt_plot, width = plot_width, height = 4.75, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(signal_cwt_plot)), ".pdf"), plot = signal_cwt_plot, width = plot_width, height = 4.75, units = "cm")

save_data(data = p_308_wavelet_depth_hr, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_wavelet_depth_hr)), ".png"), plot = p_308_wavelet_depth_hr, width = plot_halfpagewidth, height = 4.75, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_wavelet_depth_hr)), ".pdf"), plot = p_308_wavelet_depth_hr, width = plot_halfpagewidth, height = 4.75, units = "cm")

save_data(data = p_321_wavelet_depth_hr, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_321_wavelet_depth_hr)), ".png"), plot = p_321_wavelet_depth_hr, width = plot_halfpagewidth, height = 4.75, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_321_wavelet_depth_hr)), ".pdf"), plot = p_321_wavelet_depth_hr, width = plot_halfpagewidth, height = 4.75, units = "cm")

### depthlogs ####

# tag 308
save_data(data = p_308_depth, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_depth)), ".png"), plot = p_308_depth, width = plot_halfpagewidth, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_depth)), ".pdf"), plot = p_308_depth, width = plot_halfpagewidth, height = plot_height, units = "cm")


# tag 321
save_data(data = p_321_depth, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_321_depth)), ".png"), plot = p_321_depth, width = plot_halfpagewidth, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_321_depth)), ".pdf"), plot = p_321_depth, width = plot_halfpagewidth, height = plot_height, units = "cm")

#### summer winter #####

# tag 308
save_data(data = p_308_depth_summer, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_depth_summer)), ".png"), plot = p_308_depth_summer, width = plot_quarterpagewidth, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_depth_summer)), ".pdf"), plot = p_308_depth_summer, width = plot_quarterpagewidth, height = plot_height, units = "cm")

save_data(data = p_308_depth_winter, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_depth_winter)), ".png"), plot = p_308_depth_winter, width = plot_quarterpagewidth, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_depth_winter)), ".pdf"), plot = p_308_depth_winter, width = plot_quarterpagewidth, height = plot_height, units = "cm")

save_data(data = p_308_depth_winter_70m, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_depth_winter_70m)), ".png"), plot = p_308_depth_winter_70m, width = 12, height = 5.5, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_308_depth_winter_70m)), ".pdf"), plot = p_308_depth_winter_70m, width = 12, height = 5.5, units = "cm")


# tag 321
save_data(data = p_321_depth_summer, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_321_depth_summer)), ".png"), plot = p_321_depth_summer, width = plot_quarterpagewidth, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_321_depth_summer)), ".pdf"), plot = p_321_depth_summer, width = plot_quarterpagewidth, height = plot_height, units = "cm")

save_data(data = p_321_depth_winter, folder = path_maps)
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_321_depth_winter)), ".png"), plot = p_321_depth_winter, width = plot_quarterpagewidth, height = plot_height, units = "cm")
ggplot2::ggsave(filename = paste0(path_maps, deparse(substitute(p_321_depth_winter)), ".pdf"), plot = p_321_depth_winter, width = plot_quarterpagewidth, height = plot_height, units = "cm")



# old ####
# signal_24h <- signal %>% 
#   dplyr::filter(date_time > ("2020-01-01 12:00:00" %>% as.POSIXct(tz = "UTC")) + lubridate::days(15)) %>%
#   mutate(depth_m = 10 * (sin((2 * pi)*t) - 1))
# 
# signal_12h_24h <- signal %>%
#   left_join(signal_24h, by = c("date_time", "depth_m", "t", "tag_serial_number")) # funzt noch nicht
# 
# # NEXT STEP: INPUT PERIOD OF 24H FOR HALF THE SIGNAL
# 
# # plot signal
# ggplot() +
#   geom_line(data = signal_12h_24h, aes(x = date_time, y = depth_m)) +
#   labs(x = "Date", y = "Depth in m")