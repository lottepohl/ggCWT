# Script with function to plot a wavelet analysis output from `make_wavelet_df()` as a ggplot object
# Author: Lotte Pohl, Date: 2023-08-08

#' function to plot a wavelet analysis output from `make_wavelet_df()` as a ggplot object.
#' 
#' @author Lotte Pohl (lotte.pohl@@gmail.com)
#'
#' @param wavelet_df dataframe, for example an output from `make_wavelet_df()` that contains the periods and time stamps / dates of a wavelet result biwavelet object.
#' @param date logical operator stating if the time should be plotted as time steps or as dates.
#' @param max_period the maximum period that should be included in the ggplot object.

library(scales)
library(dplyr)
library(ggplot2)

# test
# wavelet_df <- signal_cwt_df
# date <- TRUE
# max_period <- 4200

# # ## test xwt
# wavelet_df <- signals_xwt_df
# date <- T
# # max_period <- wavelet_df %>% select(period) %>% max()
# max_period <- NULL

# test again
# wavelet_df <- xwt_both_12_24_df
# max_period <- wavelet_df %>% select(period) %>% max()

# TODO: fix bug with `max_period`: when piping `ggplot_wavelet()`, `max_period` should be able to be left emtpy and not return errors
# potentially has to do with the übergabe von dem dataframe an die funktion, während der Ausführung kann 
# wsl nicht direkt drauf zugegriffen werden, deswegen kann das max nicht rausgezogen werden

ggplot_wavelet <- function(wavelet_df, date = TRUE, max_period = NULL, opacity = 0.6){
  # transformation function for the y axis
  my_trans <- scales::trans_new("log2_reverse", function(x) -log2(x), function(x) 2^-x)
  
  # change or leave out to make more general
  wavelet_df <- wavelet_df %>% dplyr::filter(date %>% dplyr::between((min(wavelet_df$date, na.rm = T) + lubridate::days(7)),
                                                                     (max(wavelet_df$date, na.rm = T) - lubridate::days(7)))) # cut first and last week of data off to avoid looking at tagging effect
  # wavelet_df <- wavelet_df %>% dplyr::filter(date > (min(wavelet_df$date, na.rm = T) + lubridate::days(7))) 
  
  # y axis labels
  y_breaks <- 2^floor(log2(wavelet_df$period)) %>% unique()
  
  # filter out the undesired periods, if specified
  if(!base::is.null(max_period)){
    wavelet_df <- wavelet_df %>% dplyr::filter(period <= max_period)
    y_breaks <- y_breaks[y_breaks <= max_period]
  }
  
  # change max and min date to include max x axis label completely --> make customisable
  max_date <- max(wavelet_df$date, na.rm = T) + lubridate::days(10)
  min_date <- min(wavelet_df$date, na.rm = T) -lubridate::days(10)
  
  ifelse(date %>% base::isTRUE(),
                  # for now: plot only power_log
                  plot <- ggplot2::ggplot(data = wavelet_df) +
                  geom_tile(aes(x = date, y = period, fill = power_log),
                            position = "identity",
                            alpha = opacity) + 
                  # plot the significant periods fully opaque
                  geom_tile(data = wavelet_df %>% dplyr::filter(sig == 1),
                            aes(x = date, y = period, fill = power_log),
                            position = "identity") +
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
                  theme(legend.position = "bottom",
                        legend.box = "horizontal",
                        legend.margin = margin(t = -15)) +
                  theme(axis.text.x = element_text(angle = 15, hjust = 0.5))
         
         # if(!is.null(signif_level))
         
         , # if date != TRUE
         
                plot <- ggplot(data = wavelet_df) +
                  geom_tile(aes(x = t, y = period, fill = power_log),
                            position = "identity",
                            alpha = opacity) + 
                  # plot the significant periods fully opaque
                  geom_tile(data = wavelet_df %>% dplyr::filter(sig == 1), 
                            aes(x = t, y = period, fill = power_log),
                            position = "identity") +
                  scale_y_continuous(trans = my_trans,
                                     breaks = y_breaks, 
                                     expand = c(0,0)) +
                  # scale_x_discrete(breaks = x_breaks) +
                  scale_x_datetime(date_breaks = "6 weeks", date_labels = "%b %d", expand = c(0,0)) +
                  scale_fill_viridis_c(direction = 1, option = "turbo") +
                  labs(x = "time step", y = "period in hours", fill = "log2(power)") #+
                # theme(axis.text.x = element_text(angle = 60, hjust = 0.5))
  )
  
  return(plot)
}
