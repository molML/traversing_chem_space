
#### Supplementary ####

df <- read_csv("figures/data/processed_results.csv")
df = subset(df, n_start <= 64)

df$seed = factor(df$seed, levels=unique(df$seed))
df$bias = stringr::str_to_title(gsub('random', 'innate', df$bias))
df$bias = factor(df$bias, levels=c('Innate', 'Small', 'Large'))
df$batch_size = factor(df$batch_size, levels=c(16, 32, 64))
df$n_start = factor(df$n_start, levels=c(2, 4, 8, 16, 32, 64))
df[df$retrain == "FALSE", ]$acquisition_method = 'Exploit, no retraining'
df$acquisition_method = gsub('bald', 'BALD (least mutual information)', df$acquisition_method)
df$acquisition_method = gsub('exploitation', 'Exploit (best predictions)', df$acquisition_method)
df$acquisition_method = gsub('exploration', 'Explore (most uncertain)', df$acquisition_method)
df$acquisition_method = gsub('random', 'Random', df$acquisition_method)
df$acquisition_method = gsub('similarity', 'Similarity', df$acquisition_method)





##### S Fig 1 #######

dfs1a = subset(df, train_cycle == 15 & architecture == 'mlp' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs1a$batch_size = as.numeric(as.character(dfs1a$batch_size))

sfig1a = ggplot(dfs1a, aes(x = batch_size, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0, 6.75), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.5), "cm"))


dfs1b = subset(df, train_cycle == 15 & architecture == 'gcn' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs1b$batch_size = as.numeric(as.character(dfs1b$batch_size))

sfig1b = ggplot(dfs1b, aes(x = batch_size, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1, linetype='dashed') +
  labs(y = '', x='') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0, 6.75), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"))


dfs1c = subset(df, train_cycle == 15 & architecture == 'mlp' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs1c$batch_size = as.numeric(as.character(dfs1c$batch_size))

sfig1c = ggplot(dfs1c, aes(x = batch_size, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1) +
  labs(y = 'Test set ROC AUC after\n screening 1000 molecules', x='Molecules acquired per cycle') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0.45, 0.65), expand=F) +
  scale_y_continuous(breaks = seq(0.45,0.65, by=0.05)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))


dfs1d = subset(df, train_cycle == 15 & architecture == 'gcn' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs1d$batch_size = as.numeric(as.character(dfs1d$batch_size))

sfig1d = ggplot(dfs1d, aes(x = batch_size, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1, linetype='dashed') +
  labs(y = '', x='Molecules acquired per cycle') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0.45, 0.65), expand=F) +
  scale_y_continuous(breaks = seq(0.45,0.65, by=0.05)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


sfig1 = plot_grid(sfig1a, sfig1b, sfig1c, sfig1d,
                  ncol=2, label_size = 8, labels= c('a', 'b', 'c', 'd'))


pdf('figures/sfig1.pdf', width = 88/25.4, height = 80/25.4)
print(sfig1)
dev.off()


##### S Fig 2 #######


umap_aldh1 <- read_csv("figures/data/UMAP_ALDH1.csv")
umap_aldh1 = subset(umap_aldh1, split == 's')
umap_aldh1 = umap_aldh1[order(umap_aldh1$y),]
umap_aldh1$size = umap_aldh1$y
umap_aldh1$size[umap_aldh1$size == 1] = 0.25
umap_aldh1$size[umap_aldh1$size == 0] = 1

p_umap_aldh1 = ggplot(umap_aldh1, aes(x=umap_x, y=umap_y, color = as.character(y)))+
  scale_color_manual(values=c("#cccccc", "#005f73")) +
  geom_point(size=umap_aldh1$size, alpha=0.5) +
  theme(
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'None',
    axis.line.y.left=element_line(color="#1e3648", size=0.5, linetype = 'dashed'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

umap_pkm2 <- read_csv("figures/data/UMAP_PKM2.csv")
umap_pkm2 = subset(umap_pkm2, split == 's')
umap_pkm2 = umap_pkm2[order(umap_pkm2$y),]
umap_pkm2$size = umap_pkm2$y
umap_pkm2$size[umap_pkm2$size == 1] = 0.25
umap_pkm2$size[umap_pkm2$size == 0] = 1

p_umap_pkm2 = ggplot(umap_pkm2, aes(x=umap_x, y=umap_y, color = as.character(y)))+
  scale_color_manual(values=c("#cccccc", "#005f73")) +
  geom_point(size=umap_aldh1$size, alpha=0.5) +
  theme(
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'None',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

umap_vdr <- read_csv("figures/data/UMAP_VDR.csv")
umap_vdr = subset(umap_vdr, split == 's')
umap_vdr = umap_vdr[order(umap_vdr$y),]
umap_vdr$size = umap_vdr$y
umap_vdr$size[umap_vdr$size == 1] = 0.25
umap_vdr$size[umap_vdr$size == 0] = 1

p_umap_vdr = ggplot(umap_vdr, aes(x=umap_x, y=umap_y, color = as.character(y)))+
  scale_color_manual(values=c("#cccccc", "#005f73")) +
  geom_point(size=umap_aldh1$size, alpha=0.5) +
  theme(
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line.y.left=element_line(color="#1e3648", size=0.5, linetype = 'dashed'),
    legend.position = 'None',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


sfig2 = plot_grid(p_umap_pkm2, p_umap_aldh1, p_umap_vdr, rel_widths = c(2, 2, 2),
                  labels = c('a', 'b', 'c'),
                  ncol=3, label_size=8)

pdf('figures/sfig2.pdf', width = 180/25.4, height = 60/25.4)
print(sfig2)
dev.off()


##### S Fig 3 #####

dfs3 <- read_csv("figures/data/similarity_histograms.csv")

sfig3a = ggplot(subset(dfs3, dataset == 'ALDH1'), aes(x=`Tanimoto similarity`, fill=kind)) +
  geom_density(alpha=0.4, linewidth = 0.45) +
  scale_fill_manual(values=c("#5c8095", "#dddddd")) +
  labs(x='Mean Tanimoto similarity', fill='', title='ALDH1') +
  scale_x_continuous(limits = c(0,0.3), expand = expansion(mult = c(0.01, 0.01)))+
  default_plot_theme + theme(legend.position = 'none')

sfig3b = ggplot(subset(dfs3, dataset == 'PKM2'), aes(x=`Tanimoto similarity`, fill=kind)) +
  geom_density(alpha=0.4, linewidth = 0.45) +
  scale_fill_manual(values=c("#5c8095", "#dddddd")) +
  labs(x='Mean Tanimoto similarity', fill='', title='PKM2') +
  scale_x_continuous(limits = c(0,0.3), expand = expansion(mult = c(0.01, 0.01)))+
  default_plot_theme + theme(legend.position = 'none')

sfig3c = ggplot(subset(dfs3, dataset == 'VDR'), aes(x=`Tanimoto similarity`, fill=kind)) +
  geom_density(alpha=0.4, linewidth = 0.45) +
  scale_fill_manual(values=c("#5c8095", "#dddddd")) +
  labs(x='Mean Tanimoto similarity', fill='', title='VDR') +
  scale_x_continuous(limits = c(0,0.3), expand = expansion(mult = c(0.01, 0.01)))+
  default_plot_theme + theme(legend.position = 'right')


sfig3 = plot_grid(sfig3a, sfig3b, sfig3c, plot_spacer(), ncol = 4, rel_widths = c(1, 1, 1.7, 0.2),
                  labels = c('a', 'b', 'c', ''), label_size=8)


pdf('figures/sfig3.pdf', width = 180/25.4, height = 45/25.4)
print(sfig3)
dev.off()


##### S Fig 4 #######

###### first row ######

dfs4ac = subset(df, train_cycle == 15 & batch_size == 64 & architecture == 'mlp' & dataset == 'PKM2' & n_start == 64)
mean_sim_per_method = subset(df, train_cycle == 0 & batch_size == 64 & architecture == 'mlp' & dataset == 'PKM2' & n_start == 64)  # $mean_total_sims

starting_bias = c()
for (i_row in 1:nrow(dfs4ac)){
  row = dfs4ac[i_row, ]
  mean_sim0 = subset(mean_sim_per_method, acquisition_method == row$acquisition_method & seed == row$seed & bias == row$bias)$mean_total_sims
  starting_bias = c(starting_bias, mean_sim0)
}
dfs4ac$starting_bias = starting_bias

figs4a = ggplot(dfs4ac, aes(x = starting_bias, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.5), "cm"))

figs4c = ggplot(dfs4ac, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1) +
  labs(y = 'Test set ROC AUC after\n screening 1000 molecules', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0.35, 0.75), expand=F) + 
  scale_y_continuous(breaks = seq(0.35,0.75, by=0.05)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))


dfs4bd = subset(df, train_cycle == 15 & batch_size == 64 & architecture == 'gcn' & dataset == 'PKM2' & n_start == 64)
mean_sim_per_method = subset(df, train_cycle == 0 & batch_size == 64 & architecture == 'gcn' & dataset == 'PKM2' & n_start == 64)  # $mean_total_sims

starting_bias = c()
for (i_row in 1:nrow(dfs4bd)){
  row = dfs4bd[i_row, ]
  mean_sim0 = subset(mean_sim_per_method, acquisition_method == row$acquisition_method & seed == row$seed & bias == row$bias)$mean_total_sims
  starting_bias = c(starting_bias, mean_sim0)
}
dfs4bd$starting_bias = starting_bias

figs4b = ggplot(dfs4bd, aes(x = starting_bias, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1, linetype='dashed') +
  labs(y = '', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"))

figs4d = ggplot(dfs4bd, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1, linetype='dashed') +
  labs(y = '', x='') +
  coord_cartesian(xlim = c(0.1, 0.525), ylim = c(0.35, 0.75), expand=F) + 
  scale_y_continuous(breaks = seq(0.35,0.75, by=0.05)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


###### second row ######

dfs4eg = subset(df, train_cycle == 15 & batch_size == 64 & architecture == 'mlp' & dataset == 'ALDH1' & n_start == 64)
mean_sim_per_method = subset(df, train_cycle == 0 & batch_size == 64 & architecture == 'mlp' & dataset == 'ALDH1' & n_start == 64)  # $mean_total_sims

starting_bias = c()
for (i_row in 1:nrow(dfs4eg)){
  row = dfs4eg[i_row, ]
  mean_sim0 = subset(mean_sim_per_method, acquisition_method == row$acquisition_method & seed == row$seed & bias == row$bias)$mean_total_sims
  starting_bias = c(starting_bias, mean_sim0)
}
dfs4eg$starting_bias = starting_bias


figs4e = ggplot(dfs4eg, aes(x = starting_bias, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.5), "cm"))

figs4g = ggplot(dfs4eg, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1) +
  labs(y = 'Test set ROC AUC after\n screening 1000 molecules', x='') +
  coord_cartesian(xlim = c(0.1, 0.525), ylim = c(0.35, 0.75), expand=F) + 
  scale_y_continuous(breaks = seq(0.35,0.75, by=0.05)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))


dfs4fh = subset(df, train_cycle == 15 & batch_size == 64 & architecture == 'gcn' & dataset == 'ALDH1' & n_start == 64)
mean_sim_per_method = subset(df, train_cycle == 0 & batch_size == 64 & architecture == 'gcn' & dataset == 'ALDH1' & n_start == 64)  # $mean_total_sims

starting_bias = c()
for (i_row in 1:nrow(dfs4fh)){
  row = dfs4fh[i_row, ]
  mean_sim0 = subset(mean_sim_per_method, acquisition_method == row$acquisition_method & seed == row$seed & bias == row$bias)$mean_total_sims
  starting_bias = c(starting_bias, mean_sim0)
}
dfs4fh$starting_bias = starting_bias

figs4f = ggplot(dfs4fh, aes(x = starting_bias, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1, linetype='dashed') +
  labs(y = '', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"))

figs4h = ggplot(dfs4fh, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1, linetype='dashed') +
  labs(y = '', x='') +
  coord_cartesian(xlim = c(0.1, 0.525), ylim = c(0.35, 0.75), expand=F) + 
  scale_y_continuous(breaks = seq(0.35,0.75, by=0.05)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


###### third row ######

dfs4ik= subset(df, train_cycle == 15 & batch_size == 64 & architecture == 'mlp' & dataset == 'VDR' & n_start == 64)
mean_sim_per_method = subset(df, train_cycle == 0 & batch_size == 64 & architecture == 'mlp' & dataset == 'VDR' & n_start == 64)  # $mean_total_sims

starting_bias = c()
for (i_row in 1:nrow(dfs4ik)){
  row = dfs4ik[i_row, ]
  mean_sim0 = subset(mean_sim_per_method, acquisition_method == row$acquisition_method & seed == row$seed & bias == row$bias)$mean_total_sims
  starting_bias = c(starting_bias, mean_sim0)
}
dfs4ik$starting_bias = starting_bias

figs4i = ggplot(dfs4ik, aes(x = starting_bias, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='Structural bias in start data') +
  coord_cartesian(xlim = c(0.12, 0.4), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7.5, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.5), "cm"))

figs4k = ggplot(dfs4ik, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1) +
  labs(y = 'Test set ROC AUC after\n screening 1000 molecules', x='Structural bias in start data') +
  coord_cartesian(xlim = c(0.12, 0.4), ylim = c(0.35, 0.75), expand=F) + 
  scale_y_continuous(breaks = seq(0.35,0.75, by=0.05)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))

dfs4jl = subset(df, train_cycle == 15 & batch_size == 64 & architecture == 'gcn' & dataset == 'VDR' & n_start == 64)
mean_sim_per_method = subset(df, train_cycle == 0 & batch_size == 64 & architecture == 'gcn' & dataset == 'VDR' & n_start == 64)  # $mean_total_sims

starting_bias = c()
for (i_row in 1:nrow(dfs4jl)){
  row = dfs4jl[i_row, ]
  mean_sim0 = subset(mean_sim_per_method, acquisition_method == row$acquisition_method & seed == row$seed & bias == row$bias)$mean_total_sims
  starting_bias = c(starting_bias, mean_sim0)
}
dfs4jl$starting_bias = starting_bias

figs4j = ggplot(dfs4jl, aes(x = starting_bias, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1, linetype='dashed') +
  labs(y = '', x='Structural bias in start data') +
  coord_cartesian(xlim = c(0.12, 0.4), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7.5, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"))

figs4l = ggplot(dfs4jl, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.45, alpha=0.1, linetype='dashed') +
  labs(y = '', x='Structural bias in start data') +
  coord_cartesian(xlim = c(0.12, 0.4), ylim = c(0.35, 0.75), expand=F) + 
  scale_y_continuous(breaks = seq(0.35,0.75, by=0.05)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


figs4 = plot_grid(figs4a, figs4b, figs4c,  figs4d, 
                  figs4e, figs4f, figs4g,  figs4h, 
                  figs4i, figs4j, figs4k,  figs4l, 
                  ncol=4, label_size = 8, 
                  labels= c('a', 'b', 'c', 'd',
                            'e', 'f', 'g', 'h', 
                            'i', 'j', 'k', 'l'))

pdf('figures/sfig4.pdf', width = 180/25.4, height = 130/25.4)
print(figs4)
dev.off()


##### S Fig 5 #######

dfs5 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64 & scrambledx == TRUE)


dfs5 = dfs5 %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc',
                     'test_balanced_accuracy', 'enrichment'),
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

sfig5a = ggplot(subset(dfs5, dataset == 'PKM2'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
# u r b l
sfig5b = ggplot(subset(dfs5, dataset == 'ALDH1'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfig5c = ggplot(subset(dfs5, dataset == 'VDR'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfig5d = ggplot(subset(dfs5, dataset == 'FEN1'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfig5e = ggplot(subset(dfs5, dataset == 'GBA'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig5f = ggplot(subset(dfs5, dataset == 'KAT2A'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig5g = ggplot(subset(dfs5, dataset == 'IDH1'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfig5h = ggplot(subset(dfs5, dataset == 'OPRK1'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

sfig5i = ggplot(subset(dfs5, dataset == 'ADRB2'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) +
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))


sfig5 = plot_grid(sfig5a, sfig5b, sfig5c,
                     sfig5d, sfig5e, sfig5f,
                     sfig5g, sfig5h, sfig5i,
                     labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' ,'i'),
                     ncol=3, label_size=8)

pdf('figures/sfig5.pdf', width = 120/25.4, height = 105/25.4)
print(sfig5)
dev.off()


##### S Fig 6 ######

dfs6 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64) %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

sfig6a = ggplot(subset(dfs6, dataset == 'PKM2' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean ROC AUC\non the test set', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

sfig6b = ggplot(subset(dfs6, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfig6c = ggplot(subset(dfs6, dataset == 'PKM2' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfig6d = ggplot(subset(dfs6, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean ROC AUC\non the test set', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfig6e = ggplot(subset(dfs6, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig6f = ggplot(subset(dfs6, dataset == 'ALDH1' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig6g = ggplot(subset(dfs6, dataset == 'VDR' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean ROC AUC\non the test set', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfig6h = ggplot(subset(dfs6, dataset == 'VDR' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

sfig6i = ggplot(subset(dfs6, dataset == 'VDR' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

sfig6 = plot_grid(sfig6a, sfig6b, sfig6c,
                  sfig6d, sfig6e, sfig6f,
                  sfig6g, sfig6h, sfig6i,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'), 
                  ncol=3, label_size=8)

pdf('figures/sfig6.pdf', width = 110/25.4, height = 90/25.4)
print(sfig6)
dev.off()


##### S Fig 7 #####

df7 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64) %>%
  group_by(architecture, dataset, acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

sfig7a = ggplot(subset(df7, architecture == 'mlp' & dataset == 'PKM2'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7b = ggplot(subset(df7, architecture == 'gcn' & dataset == 'PKM2'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7c = ggplot(subset(df7, architecture == 'rf' & dataset == 'PKM2'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7d = ggplot(subset(df7, architecture == 'mlp' & dataset == 'ALDH1'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7e = ggplot(subset(df7, architecture == 'gcn' & dataset == 'ALDH1'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7f = ggplot(subset(df7, architecture == 'rf' & dataset == 'ALDH1'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7g = ggplot(subset(df7, architecture == 'mlp' & dataset == 'VDR'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7h = ggplot(subset(df7, architecture == 'gcn' & dataset == 'VDR'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7i = ggplot(subset(df7, architecture == 'rf' & dataset == 'VDR'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig7 = plot_grid(sfig7a, sfig7b, sfig7c,
                  sfig7d, sfig7e, sfig7f,
                  sfig7g, sfig7h, sfig7i,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'),
                  ncol=3, label_size=8)

pdf('figures/sfig7.pdf', width = 110/25.4, height = 100/25.4)
print(sfig7)
dev.off()


##### S Fig 8 #####

dfs8 = subset(df, total_mols_screened == 1000 & batch_size == 64 & bias == 'Innate')
dfs8 = subset(dfs8, acquisition_method %in% c("Random", "Similarity", "BALD (least mutual information)"))

dfs8 = dfs8 %>%
  group_by(acquisition_method, n_start, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

sfig8_colours = c("#5c8095", "#bbbbbb", "#efc57b")

sfig8a = ggplot(subset(dfs8, dataset == 'PKM2' & architecture == 'mlp'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.45, alpha=0.5) + 
  geom_point(size=0.75) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=sfig8_colours) +
  scale_fill_manual(values=sfig8_colours) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
# u r b l
sfig8b = ggplot(subset(dfs8, dataset == 'PKM2' & architecture == 'gcn'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.45, alpha=0.5) + 
  geom_point(size=0.75) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=sfig8_colours) +
  scale_fill_manual(values=sfig8_colours) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfig8c = ggplot(subset(dfs8, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.45, alpha=0.5) + 
  geom_point(size=0.75) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=sfig8_colours) +
  scale_fill_manual(values=sfig8_colours) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfig8d = ggplot(subset(dfs8, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.45, alpha=0.5) + 
  geom_point(size=0.75) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=sfig8_colours) +
  scale_fill_manual(values=sfig8_colours) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig8e = ggplot(subset(dfs8, dataset == 'VDR' & architecture == 'mlp'),  aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.45, alpha=0.5) + 
  geom_point(size=0.75) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='molecules in the start set', color = 'Method', fill = 'Method')+
  scale_color_manual(values=sfig8_colours) +
  scale_fill_manual(values=sfig8_colours) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfig8f = ggplot(subset(dfs8, dataset == 'VDR' & architecture == 'gcn'),  aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.45, alpha=0.5) + 
  geom_point(size=0.75) +
  labs(y = '', x='molecules in the start set', color = 'Method', fill = 'Method')+
  scale_color_manual(values=sfig8_colours) +
  scale_fill_manual(values=sfig8_colours) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))


sfig8 = plot_grid(sfig8a, sfig8b, 
                  sfig8c, sfig8d,
                  sfig8e, sfig8f,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
                  ncol=2, label_size=8)

# 180 mm/ 88 mm
pdf('figures/sfig8.pdf', width = 88/25.4, height = 123/25.4)
print(sfig8)
dev.off()


##### S Fig 9 #######


dfs9 = df %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc',
                     'test_balanced_accuracy', 'enrichment'),
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

# first row
sfig9a = ggplot(subset(dfs9, dataset == 'FEN1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

sfig9b = ggplot(subset(dfs9, dataset == 'FEN1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfig9c = ggplot(subset(dfs9, dataset == 'GBA' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

sfig9d = ggplot(subset(dfs9, dataset == 'GBA' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

# second row
sfig9e = ggplot(subset(dfs9, dataset == 'KAT2A' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfig9f = ggplot(subset(dfs9, dataset == 'KAT2A' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfig9g = ggplot(subset(dfs9, dataset == 'IDH1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig9h = ggplot(subset(dfs9, dataset == 'IDH1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

# third row
sfig9i = ggplot(subset(dfs9, dataset == 'OPRK1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfig9j = ggplot(subset(dfs9, dataset == 'OPRK1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfig9k = ggplot(subset(dfs9, dataset == 'ADRB2' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

sfig9l = ggplot(subset(dfs9, dataset == 'ADRB2' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))


sfig9 = plot_grid(sfig9a, sfig9b, sfig9c, sfig9d,
                  sfig9e, sfig9f, sfig9g, sfig9h,
                  sfig9i, sfig9j, sfig9k, sfig9l,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'),
                  ncol=4, label_size=8)

pdf('figures/sfig9.pdf', width = 180/25.4, height = 123/25.4)
print(sfig9)
dev.off()


##### S Fig 10 #######

dfs10 = subset(df, batch_size == 64 & bias == 'Innate' & total_mols_screened==1000 & n_start==64) %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

dataset_yield =  data.frame(dataset = c('ALDH1', 'VDR', 'PKM2', 'FEN1', 'GBA', 'KAT2A', 'IDH1', 'OPRK1', 'ADRB2'),
                            yield = c(4986, 239, 223, 100, 55 ,54, 11, 9 ,5))
dataset_yield$label = paste0(dataset_yield$dataset, '\n(', dataset_yield$yield, ' hits)')
dfs10$yield = dataset_yield$yield[match(dfs10$dataset, dataset_yield$dataset)]
dfs10$yield_label = dataset_yield$label[match(dfs10$dataset, dataset_yield$dataset)]
x_axis_values = round(sort(unique(log(dfs10$yield))), 2)

sfig10a = ggplot(subset(dfs10, architecture == 'mlp'), aes(x = log(yield), y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture)) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='log(dataset yield)', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_smooth(method=lm, se=T, fullrange=TRUE, alpha=0.2, level=0.95, size=0.45)+
  coord_cartesian(xlim = c(x_axis_values[1], 8.6), ylim = c(0, 8), expand=F) +
  scale_y_continuous(breaks = seq(0,8, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))

sfig10b = ggplot(subset(dfs10, architecture == 'gcn'), aes(x = log(yield), y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture)) +
  labs(y = '', x='log(dataset yield)', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_smooth(method=lm, se=T, fullrange=TRUE, alpha=0.2, level=0.95, size=0.45)+
  coord_cartesian(xlim = c(x_axis_values[1], 8.6), ylim = c(0, 8), expand=F) +
  scale_y_continuous(breaks = seq(0,8, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0), "cm"))

sfig10 = plot_grid(sfig10a, sfig10b, 
                   labels = c('a', 'b'), 
                   ncol=2, label_size=8)

pdf('figures/sfig10.pdf', width = 88/25.4, height = 44/25.4)
print(sfig10)
dev.off()




