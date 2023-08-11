# old script to generate random signals


## rand freq for migration ####

x <- -2:4

# Generate the sum of sinusoidal components
# num_components <- 100  # Adjust the number of components as desired
signal_depth <- numeric(length(t))

for(i in 1:length(x)){
  print(x[i])
  component <- 10 * (sin(((2^x) * pi)*t) - 1)
  signal_depth <- signal_depth + component
}

signal_rand <- dplyr::tibble(t = seq(0,signal_length_days - f, f),
                             depth_m = signal_depth,
                             date_time = seq(signal_start_date,
                                             signal_start_date + lubridate::days(signal_length_days),
                                             length.out = base::length(t)))


# for (i in 1:num_components) {
#   freq <- runif(1, min = 0.0001, max = 1)  # Random frequency between 0.001 and 0.1
#   amplitude <- runif(1, min = 1, max = 10)  # Random amplitude between 1 and 10
#   phase <- runif(1, min = 0, max = 2 * pi)  # Random phase between 0 and 2*pi
#   
#   component1 <- amplitude * sin(2 * pi * freq * t + phase)
#   component2 <- amplitude * (sin(2 * pi * freq*t + phase) - amplitude)
#   signal <- signal + component
#   signal2 <- signal2 + component2
# }
# 
# # Create the dataframe
# signal_df <- tibble(
#   t = t,
#   # depth_m_3 = 10 * (sin((4 * pi) * t) - 1) + signal,
#   depth_m_2 = signal2,
#   depth_m = signal,
#   date_time = seq(
#     signal_start_date,
#     signal_start_date + lubridate::days(signal_length_days),
#     length.out = length(t)
#   )
# )

# head(signal_df)
# In this modified code, we generate a sum of num_components sinusoidal components with random frequencies, amplitudes, and phases. These components are added to the original sine wave. You can adjust the num_components and the ranges for frequency, amplitude, and phase according to your preferences. The resulting signal_df dataframe will contain the time, depth, and date-time columns for the generated signal.
# 
# Please note that generating a high number of random frequencies can create complex and potentially noisy signals. You may want to adjust the parameters and add additional logic to control the characteristics of the generated signal according to your specific requirements.

ggplot(data = signal_rand %>% dplyr::filter(date_time %>% between((signal_start_date + lubridate::days(signal_length_days / 2)) - lubridate::days(0), 
                                                                  (signal_start_date + lubridate::days(signal_length_days / 2)) + lubridate::days(1)))) +
  geom_line(aes(x = date_time, y = depth_m), colour = 'red') #+
# geom_line(aes(x = date_time, y = depth_m_2), colour = 'blue') #+
# geom_line(aes(x = date_time, y = depth_m_3), colour = 'green') 


rand_signal_cwt <- compute_CWT(values = signal_rand %>% dplyr::select(depth_m),
                               dt = dt_hours) %>%
  make_wavelet_df(date_times = signal_df %>% dplyr::select(date_time),
                  dt = dt_hours) %>%
  ggplot_wavelet(max_period = 500)

rand_signal_cwt

## FFT rand signal ####

fft_rand <- signal_rand %>% calc_fft(sample_freq = sampling_rate_hours) %>%
  plot_periodogram(period_upperlim = 500)

fft_rand
