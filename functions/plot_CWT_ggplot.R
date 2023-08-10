# Script with function to plot a wavelet analysis output from `make_Ccwt_df()` as a ggplot object
# Author: Lotte Pohl, Date: 2023-08-08

#' function to plot a wavelet analysis output from `make_CWT_df()` as a ggplot object.
#' 
#' @author Lotte Pohl (lotte.pohl@@gmail.com)
#'
#' @param cwt_df dataframe, for example an output from `make_CWT_df()` that contains the periods and time stamps / dates of a CWT result biwavelet object.
#' @param date logical operator stating if the time should be plotted as time steps or as dates.
#' @param max_period the maximum period that should be included in the ggplot object.

library(scales)
library(dplyr)
library(ggplot2)

# test
# cwt_df <- signal_CWT_df
# date <- TRUE
# max_period <- 4200

# ## test xwt
# cwt_df <- signals_xwt_df
# date <- T
# max_period <- cwt_df %>% select(period) %>% max()


plot_CWT_ggplot <- function(cwt_df, date = TRUE, max_period){
  # transformation function for the y axis
  my_trans <- scales::trans_new("log2_reverse", function(x) -log2(x), function(x) 2^-x)
  
  # change or leave out to make more general
  cwt_df <- cwt_df %>% dplyr::filter(date > (min(cwt_df$date) + lubridate::days(7))) # cut first week of data off to avoid looking at tagging effect
  
  # y axis labels
  y_breaks <- 2^floor(log2(cwt_df$period)) %>% unique()
  y_breaks <- y_breaks[y_breaks <= max_period]

  # filter out the undesired periods
  cwt_df <- cwt_df %>% dplyr::filter(period <= max_period) 
  
  # change max and min date to include max x axis label completely --> make customisable
  max_date <- max(cwt_df$date) + lubridate::days(10)
  min_date <- min(cwt_df$date) -lubridate::days(10)
  
  ifelse(date %>% base::isTRUE(),
                  # for now: plot only power_log
                  plot <- ggplot2::ggplot(data = cwt_df) +
                  geom_tile(aes(x = date, y = period, fill = power_log),
                            position = "identity",
                            alpha = 0.6) + 
                  # plot the significant periods fully opaque
                  geom_tile(data = cwt_df %>% dplyr::filter(sig == 1), 
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
         
         , # if date != TRUE
         
                plot <- ggplot(data = cwt_df) +
                  geom_tile(aes(x = t, y = period, fill = power_log),
                            position = "identity",
                            alpha = 0.6) + 
                  # plot the significant periods fully opaque
                  geom_tile(data = cwt_df %>% dplyr::filter(sig == 1), 
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
