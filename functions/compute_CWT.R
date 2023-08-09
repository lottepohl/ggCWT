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
    ungroup() %>%
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


# 
# plot_wavelet_hrperiod_example <- function(wt_df = wt_df, type = c("power", "significance", "power_log"),
#                                           date = TRUE, max_period = 72 #hours = 4 days
# ){
#   # transformation function for the y axis
#   my_trans <- scales::trans_new("log2_reverse", function(x) -log2(x), function(x) 2^-x)
#   
#   # wt_df <- wt_df %>% dplyr::filter(date > (min(wt_df$date) + lubridate::days(7))) # cut first week of data off to avoid looking at tagging effect
#   
#   # y axis labels
#   y_breaks <- 2^floor(log2(wt_df$period)) %>% unique()
#   y_breaks <- y_breaks[y_breaks <= max_period]
#   # transform dates
#   # wt_df$date_time <- wt_df$date_time %>% as.POSIXct(tz = "UTC") # not needed bc we have the posict objects
#   
#   # # x axis labels
#   # ifelse(date %>% isTRUE(),
#   #        x_breaks <- c(wt_df$date_time[1], wt_df$date_time[(1/5) * n_data], wt_df$date_time[(2/5) * n_data], wt_df$date_time[(3/5) * n_data],
#   #                      wt_df$date_time[(4/5) * n_data], wt_df$date_time[(5/5) * n_data])
#   # ,
#   # x_breaks <- sprintf("%03d", seq(from = 0, to = n_data, by = 100)))
#   # 
#   #plot
#   
#   wt_df <- wt_df  %>% dplyr::filter(period <= max_period)
#   
#   # change max and min date to include max x axis label completely
#   # max_date <- max(wt_df$date) + lubridate::days(10)
#   # min_date <- min(wt_df$date) -lubridate::days(10)
#   
#   
#   ifelse(date %>% isTRUE(),
#          
#          ifelse(type == "power_log",
#                 
#                 plot <- ggplot(data = wt_df) +
#                   geom_tile(aes(x = date, y = period, fill = power_log),
#                             position = "identity",
#                             alpha = 0.7) + #0.5
#                   # geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date, y = period, fill = power_log),
#                   #           position = "identity") +
#                   scale_y_continuous(trans = my_trans,
#                                      breaks = y_breaks, 
#                                      expand = c(0,0)) +
#                   # scale_x_discrete(breaks = x_breaks) +
#                   # scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                   scale_x_datetime(
#                     date_minor_breaks = "1 day",
#                     date_breaks = "1 week",
#                     date_labels = "%b %d",
#                     expand = c(0,0)
#                     # ,limits = c(min_date, max_date)
#                   ) +
#                   scale_fill_viridis_c(direction = 1, option = "turbo") +
#                   labs(x = "", y = "Period in hours", fill = "log2(Power)") +
#                   theme(legend.position = "bottom", # "bottom",
#                         legend.box = "horizontal",
#                         legend.margin = margin(t = -15)) #+
#                 # theme(axis.text.x = element_text(angle = 15, hjust = 0.5))
#                 ,
#                 
#                 ifelse(type == "significance",
#                        
#                        plot <- ggplot(data = wt_df) +
#                          geom_tile(aes(x = date_posicxt, y = period, fill = significance),
#                                    position = "identity",
#                                    alpha = 0.65) +
#                          geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date_posicxt, y = period, fill = significance),
#                                    position = "identity") +
#                          scale_y_continuous(trans = my_trans,
#                                             breaks = y_breaks, 
#                                             expand = c(0,0)) +
#                          # scale_x_discrete(breaks = x_breaks) +
#                          scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                          scale_fill_viridis_c(direction = 1, option = "turbo") +
#                          labs(x = "date", y = "period in hours", fill = "significance") #+
#                        # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                        ,
#                        
#                        plot <- ggplot(data = wt_df) +
#                          geom_tile(aes(x = date_posicxt, y = period, fill = power),
#                                    position = "identity",
#                                    alpha = 0.65) +
#                          geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date_posicxt, y = period, fill = power),
#                                    position = "identity") +
#                          scale_y_continuous(trans = my_trans,
#                                             breaks = y_breaks, 
#                                             expand = c(0,0)) +
#                          # scale_x_discrete(breaks = x_breaks) +
#                          scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                          scale_fill_viridis_c(direction = 1, option = "turbo") +
#                          labs(x = "date", y = "period in hours", fill = "power") #+
#                        # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                        
#                 )
#          )
#          ,
#          
#          ifelse(type == "power_log",
#                 
#                 plot <- ggplot(data = wt_df) +
#                   geom_tile(aes(x = t, y = period, fill = power_log),
#                             position = "identity",
#                             alpha = 0.65) +
#                   geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = power_log),
#                             position = "identity") +
#                   scale_y_continuous(trans = my_trans,
#                                      breaks = y_breaks, 
#                                      expand = c(0,0)) +
#                   # scale_x_discrete(breaks = x_breaks) +
#                   scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                   scale_fill_viridis_c(direction = 1, option = "turbo") +
#                   labs(x = "date", y = "period in hours", fill = "log2(power)") #+
#                 # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                 ,
#                 
#                 ifelse(type == "significance",
#                        
#                        plot <- ggplot(data = wt_df) +
#                          geom_tile(aes(x = t, y = period, fill = significance),
#                                    position = "identity",
#                                    alpha = 0.65) +
#                          geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = significance),
#                                    position = "identity") +
#                          scale_y_continuous(trans = my_trans,
#                                             breaks = y_breaks, 
#                                             expand = c(0,0)) +
#                          # scale_x_discrete(breaks = x_breaks) +
#                          scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                          scale_fill_viridis_c(direction = 1, option = "turbo") +
#                          labs(x = "date", y = "period in hours", fill = "significance") #+
#                        # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                        ,
#                        
#                        plot <- ggplot(data = wt_df) +
#                          geom_tile(aes(x = t, y = period, fill = power),
#                                    position = "identity",
#                                    alpha = 0.65) +
#                          geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = power),
#                                    position = "identity") +
#                          scale_y_continuous(trans = my_trans,
#                                             breaks = y_breaks, 
#                                             expand = c(0,0)) +
#                          # scale_x_discrete(breaks = x_breaks) +
#                          scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                          scale_fill_viridis_c(direction = 1, option = "turbo") +
#                          labs(x = "date", y = "period in hours", fill = "power") #+
#                        # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                 )
#          )
#   )
#   return(plot)
# }
# 
# plot_wavelet_hrperiod <- function(wt_df = wt_df, type = c("power", "significance", "power_log"),
#                                   date = TRUE, max_period = 72 #hours = 4 days
# ){
#   # transformation function for the y axis
#   my_trans <- scales::trans_new("log2_reverse", function(x) -log2(x), function(x) 2^-x)
#   
#   wt_df <- wt_df %>% dplyr::filter(date > (min(wt_df$date) + lubridate::days(7))) # cut first week of data off to avoid looking at tagging effect
#   
#   # y axis labels
#   y_breaks <- 2^floor(log2(wt_df$period)) %>% unique()
#   y_breaks <- y_breaks[y_breaks <= max_period]
#   # transform dates
#   # wt_df$date_time <- wt_df$date_time %>% as.POSIXct(tz = "UTC") # not needed bc we have the posict objects
#   
#   # # x axis labels
#   # ifelse(date %>% isTRUE(),
#   #        x_breaks <- c(wt_df$date_time[1], wt_df$date_time[(1/5) * n_data], wt_df$date_time[(2/5) * n_data], wt_df$date_time[(3/5) * n_data],
#   #                      wt_df$date_time[(4/5) * n_data], wt_df$date_time[(5/5) * n_data])
#   # ,
#   # x_breaks <- sprintf("%03d", seq(from = 0, to = n_data, by = 100)))
#   # 
#   #plot
#   
#   wt_df <- wt_df  %>% dplyr::filter(period <= max_period)
#   
#   # change max and min date to include max x axis label completely
#   max_date <- max(wt_df$date) + lubridate::days(10)
#   min_date <- min(wt_df$date) -lubridate::days(10)
#   
#   
#   ifelse(date %>% isTRUE(),
#          
#          ifelse(type == "power_log",
#                 
#                 plot <- ggplot(data = wt_df) +
#                   geom_tile(aes(x = date, y = period, fill = power_log),
#                             position = "identity",
#                             alpha = 0.7) + #0.5
#                   # geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date, y = period, fill = power_log),
#                   #           position = "identity") +
#                   scale_y_continuous(trans = my_trans,
#                                      breaks = y_breaks, 
#                                      expand = c(0,0)) +
#                   # scale_x_discrete(breaks = x_breaks) +
#                   # scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                   scale_x_datetime(
#                     # date_minor_breaks = "1 month",
#                     date_breaks = "1 month",
#                     date_labels = "%b '%y",
#                     expand = c(0,0),
#                     limits = c(min_date, max_date)) +
#                   scale_fill_viridis_c(direction = 1, option = "turbo") +
#                   labs(x = "", y = "Period in hours", fill = "log2(Power)") +
#                   theme(legend.position = "bottom", # "bottom",
#                         legend.box = "horizontal",
#                         legend.margin = margin(t = -15)) +
#                   theme(axis.text.x = element_text(angle = 15, hjust = 0.5))
#                 ,
#                 
#                 ifelse(type == "significance",
#                        
#                        plot <- ggplot(data = wt_df) +
#                          geom_tile(aes(x = date_posicxt, y = period, fill = significance),
#                                    position = "identity",
#                                    alpha = 0.65) +
#                          geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date_posicxt, y = period, fill = significance),
#                                    position = "identity") +
#                          scale_y_continuous(trans = my_trans,
#                                             breaks = y_breaks, 
#                                             expand = c(0,0)) +
#                          # scale_x_discrete(breaks = x_breaks) +
#                          scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                          scale_fill_viridis_c(direction = 1, option = "turbo") +
#                          labs(x = "date", y = "period in hours", fill = "significance") #+
#                        # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                        ,
#                        
#                        plot <- ggplot(data = wt_df) +
#                          geom_tile(aes(x = date_posicxt, y = period, fill = power),
#                                    position = "identity",
#                                    alpha = 0.65) +
#                          geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = date_posicxt, y = period, fill = power),
#                                    position = "identity") +
#                          scale_y_continuous(trans = my_trans,
#                                             breaks = y_breaks, 
#                                             expand = c(0,0)) +
#                          # scale_x_discrete(breaks = x_breaks) +
#                          scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                          scale_fill_viridis_c(direction = 1, option = "turbo") +
#                          labs(x = "date", y = "period in hours", fill = "power") #+
#                        # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                        
#                 )
#          )
#          ,
#          
#          ifelse(type == "power_log",
#                 
#                 plot <- ggplot(data = wt_df) +
#                   geom_tile(aes(x = t, y = period, fill = power_log),
#                             position = "identity",
#                             alpha = 0.65) +
#                   geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = power_log),
#                             position = "identity") +
#                   scale_y_continuous(trans = my_trans,
#                                      breaks = y_breaks, 
#                                      expand = c(0,0)) +
#                   # scale_x_discrete(breaks = x_breaks) +
#                   scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                   scale_fill_viridis_c(direction = 1, option = "turbo") +
#                   labs(x = "date", y = "period in hours", fill = "log2(power)") #+
#                 # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                 ,
#                 
#                 ifelse(type == "significance",
#                        
#                        plot <- ggplot(data = wt_df) +
#                          geom_tile(aes(x = t, y = period, fill = significance),
#                                    position = "identity",
#                                    alpha = 0.65) +
#                          geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = significance),
#                                    position = "identity") +
#                          scale_y_continuous(trans = my_trans,
#                                             breaks = y_breaks, 
#                                             expand = c(0,0)) +
#                          # scale_x_discrete(breaks = x_breaks) +
#                          scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                          scale_fill_viridis_c(direction = 1, option = "turbo") +
#                          labs(x = "date", y = "period in hours", fill = "significance") #+
#                        # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                        ,
#                        
#                        plot <- ggplot(data = wt_df) +
#                          geom_tile(aes(x = t, y = period, fill = power),
#                                    position = "identity",
#                                    alpha = 0.65) +
#                          geom_tile(data = wt_df %>% dplyr::filter(sig == 1), aes(x = t, y = period, fill = power),
#                                    position = "identity") +
#                          scale_y_continuous(trans = my_trans,
#                                             breaks = y_breaks, 
#                                             expand = c(0,0)) +
#                          # scale_x_discrete(breaks = x_breaks) +
#                          scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
#                          scale_fill_viridis_c(direction = 1, option = "turbo") +
#                          labs(x = "date", y = "period in hours", fill = "power") #+
#                        # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
#                 )
#          )
#   )
#   return(plot)
# }
# 
