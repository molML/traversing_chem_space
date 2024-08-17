

df2 <- read_csv('figures/data/fig2.csv')

df2 = df2 %>%
  group_by(acquisition_method, bias, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "mean_total_sims", 'mean_total_hits_sims', 
                     'n_unique_scaffolds', 'enrichment', 'mean_tani_per_batch',
                     'mean_tani_batch_to_start_batch', 'mean_tani_all_mols_to_start_batch'), 
                   list(mean = mean, se = se), na.rm = TRUE)) %>% ungroup()

df2$architecture = factor(df2$architecture, levels=c('mlp', 'gcn'))
df2$bias = factor(df2$bias, levels=c('Innate', 'Small', 'Large'))


fig2a = ggplot(subset(df2, acquisition_method %in% c("Exploit (best predictions)")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = 'Similarity of acquired\nbatch to start data', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.45, by=0.05)) +
  default_plot_theme + theme(plot.margin = unit(c(0.5, 0, 0, 0.25), "cm"))

fig2b = ggplot(subset(df2, acquisition_method %in% c("BALD (least mutual information)")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.45, by=0.05)) +
  default_plot_theme + theme(plot.margin = unit(c(0.5, 0.25, 0, 0.25), "cm"))

fig2c = ggplot(subset(df2, acquisition_method %in% c("Exploit, no retraining")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = 'Similarity of acquired\nbatch to start data', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.45, by=0.05)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.5, 0.25), "cm"))

fig2d = ggplot(subset(df2, acquisition_method %in% c("Similarity")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.45, by=0.05)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.5, 0.25), "cm"))



fig2 = plot_grid(fig2a, fig2b, fig2c, fig2d,
                 ncol=2, 
                 label_size = 8, labels= c('a', 'b', 'c', 'd'))

pdf("figures/fig2.pdf", width = 89/25.4, height = 82/25.4)
print(fig2)
dev.off()

