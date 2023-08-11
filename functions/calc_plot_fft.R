# compute and plot fft

# FFT ####

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

# 
# plot_fft <- function(fft_result, period_upperlim = 40, period_lowerlim = 0.05){
#   periodogram <- ggplot(data = fft_result %>% dplyr::filter(period < period_upperlim & period > period_lowerlim)) + 
#     geom_line(aes(x = period, y = spec), colour = "black") +
#     scale_y_continuous(expand = c(0,0)) +
#     scale_x_continuous(expand = c(0,0), breaks = seq(2, period_upperlim, by = 2)) +
#     labs(y = "Spectral density", x = "Period in hours")
#   return(periodogram)
# }

# plot periodogram 
plot_periodogram <- function(fft_result, period_upperlim = 40, period_lowerlim = 0.05){
  # todo: get local rule or set values for period upper and lower lim
  periodogram <- ggplot(data = fft_result %>% dplyr::filter(period < period_upperlim & period > period_lowerlim)) + 
    geom_line(aes(x = period, y = spec), colour = "black", size = 1) + #theme_bw() +
    scale_y_continuous(expand = c(0,0)) +
    # theme_bw(base_size = 12) +
    # scale_x_continuous(expand = c(0,0), breaks = seq(2, period_upperlim, by = 2)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(2, period_upperlim, by = 10)) +
    labs(y = "Spectral density", x = "Period in hours") 
  
  return(periodogram)
}

# fft_calc_plot <- function(depth_log, sample_frequency){
#   fft_res <- calc_fft(depth_log = depth_log,
#                       sample_freq = sample_frequency)
#   
#   pgram <- plot_periodogram(fft_result = fft_res)
#   pgram %>% return()
# }
