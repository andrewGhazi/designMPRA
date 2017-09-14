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

one_simulation = function(x, transcription_shift, n_samp){
  # randomly sample the most non-normal allele
  non_normal_samples = sample(worst$activity, size = n_samp, replace = TRUE) 
  
  # randomly sample a "typical" allele
  normal_samples = rnorm(n_samp, 
                         mean = worst_mean - transcription_shift, 
                         sd = .926) # using the average SD from Ulirsch and Tewhey together
  
  # compute the p-value of a t-test
  t.test(non_normal_samples, normal_samples)$p.value
}

# used to calculate effect size
pooled_SD = sqrt(((nrow(worst) - 1)*sd(worst$activity)**2 + (nrow(worst) - 1)*.926**2) / (nrow(worst)*2 - 2))

simulate_power = function(transcription_shift, sig_level = .05, n_sim = 10000){
  
  #simulate a t-test many times
  monte_carlo_p_values = map_dbl(1:n_sim, one_simulation, transcription_shift = transcription_shift, n_samp = nrow(worst))
  
  # return the fraction that pass the input significance level
  sum(monte_carlo_p_values < sig_level) / n_sim
}

sim_power = data_frame(transcription_shift = seq(-1.1, 1.1, length.out = 50),
                       theoretical_power = pwr.t.test(n = nrow(worst),
                                                      d = abs(transcription_shift) / pooled_SD,
                                                      type = 'two.sample',
                                                      alternative = 'two.sided',
                                                      sig.level = 1e-5)$power,
                       monte_carlo_power = mclapply(transcription_shift, simulate_power, mc.cores = 20, sig_level = 1e-5, n_sim = 3e6) %>% unlist)

sim_power %>%
  mutate(diff = monte_carlo_power - theoretical_power) %>%
  ggplot(aes(transcription_shift, diff)) + geom_path()


### QQ analysis for Ulirsch
u_activities = UMPRA %>% 
  mutate_at(vars(contains('K562')), depthNormalize) %>% 
  mutate(dnaMean = (K562_minP_DNA1 + K562_minP_DNA2)/2) %>% 
  filter(dnaMean > .13) %>% 
  select(-K562_minP_DNA1, -K562_minP_DNA2) %>% 
  gather(sample, depthAdjCount, K562_CTRL_minP_RNA1:K562_GATA1_minP_RNA4) %>% 
  mutate(activity = log(depthAdjCount / dnaMean))

Ulirsch_lillie = u_activities %>%
  group_by(construct, type, sample) %>% 
  nest %>% 
  mutate(n = map_int(data, ~sum(is.finite(.x$activity)))) %>% 
  filter(n > 4) %>% 
  mutate(lillie_p = map_dbl(data, ~lillie.test(.x$activity[is.finite(.x$activity)])$p.value)) %>% 
  arrange(lillie_p)

qq_lines = Ulirsch_lillie %>% 
  .[1:20,] %>% 
  unnest %>%
  unite(construct_sample_type, construct, sample, type) %>% 
  mutate(construct_sample_type = construct_sample_type %>% gsub('K562_|minP_', '', .)) %>% 
  group_by(construct_sample_type) %>% 
  summarise(act25 = quantile(activity,.25),
            act75 = quantile(activity, .75),
            norm25 = qnorm(.25),
            norm75 = qnorm(.75),
            qq_slope = (act25 - act75) / (norm25 - norm75),
            qq_int = act25 - qq_slope * norm25)

Ulirsch_lillie %>% 
  .[1:20,] %>% 
  unnest %>%
  unite(construct_sample_type, construct, sample, type) %>% 
  mutate(construct_sample_type = construct_sample_type %>% gsub('K562_|minP_', '', .)) %>% 
  ggplot() + 
  stat_qq(aes(sample = activity)) +
  facet_wrap('construct_sample_type', scales = 'free') + 
  geom_abline(aes(slope = qq_slope, intercept = qq_int), 
              data = qq_lines, 
              color = 'grey60', lty = 2)
ggsave(filename = '~/designMPRA/outputs/plots/Ulirsch_most_nonnormal_qqplots.png')

