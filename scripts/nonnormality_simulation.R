library(pwr)
library(parallel)

worst = list(allele = 'rs28519456_RC_altB', file_name = 'HepG2.r5') %>% get_nonnormal_allele()

worst_dens = worst %>% pull(activity) %>% density

worst_dens_fun = approxfun(x = worst_dens$x, y = worst_dens$y)

worst %>%
  ggplot(aes(x = activity)) + 
  geom_histogram(aes(y = ..density..), bins = 30) + 
  stat_function(fun = worst_dens_fun, color = 'dodgerblue2') + 
  ggtitle('rs28519456_RC_altB in sample HepG2.r5\nThe most highly non-normal allele in Tewhey et al., 2016')

worst_mean = mean(worst$activity)

pooled_SD = sqrt(((nrow(worst) - 1)*sd(worst$activity)**2 + (nrow(worst) - 1)*.926**2) / (nrow(worst)*2 - 2))

simulate_power = function(transcription_shift, sig_level = .05, n_sim = 10000){
  
  #simulate a t-test many times
  monte_carlo_p_values = sapply(1:n_sim, function(x){
    
    # randomly sample the most non-normal allele
    non_normal_samples = sample(worst$activity, size = nrow(worst), replace = TRUE) 
    
    # randomly sample a "typical" allele
    normal_samples = rnorm(nrow(worst), 
                           mean = worst_mean - transcription_shift, 
                           sd = .926) # using the average SD from Ulirsch and Tewhey together
    
    # compute the p-value of a t-test
    t.test(non_normal_samples, normal_samples)$p.value
  })
  
  # return how the fraction that pass the input significance level
  sum(monte_carlo_p_values < sig_level) / n_sim
}

sim_power = data_frame(transcription_shift = seq(-1.5, 1.5, length.out = 200),
                       theoretical_power = pwr.t.test(n = nrow(worst),
                                                      d = abs(transcription_shift) / pooled_SD,
                                                      type = 'two.sample',
                                                      alternative = 'two.sided')$power,
                       monte_carlo_power = mclapply(transcription_shift, simulate_power, mc.cores = 22) %>% unlist)

sim_power %>% 
  gather(power_type, power, -transcription_shift) %>% 
  ggplot(aes(transcription_shift, power, color = power_type)) + 
  geom_path() +
  scale_y_continuous(breaks = seq(0, 1, by = .1)) +
  ggtitle('Monte Carlow simulation comparing non-normal allele\nt-test results and theoretical results')
