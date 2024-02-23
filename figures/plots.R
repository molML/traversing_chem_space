
library(readr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(data.table)

#### plot themes ####

default_plot_theme = theme(
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "plain"),
  axis.text.y = element_text(size=6, face="plain", colour = "#1e3648"),
  axis.text.x = element_text(size=6, face="plain", colour = "#1e3648"),
  axis.title.x = element_text(size=7, face="plain", colour = "#1e3648"),
  axis.title.y = element_text(size=7, face="plain", colour = "#1e3648"),
  axis.ticks = element_line(color="#1e3648"),
  axis.line.x.bottom=element_line(color="#1e3648", size=0.5),
  axis.line.y.left=element_line(color="#1e3648", size=0.5),
  legend.key = element_blank(),
  legend.position = 'None',
  legend.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

default_plot_theme_legend = theme(
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "plain"),
  axis.text.y = element_text(size=6, face="plain", colour = "#1e3648"),
  axis.text.x = element_text(size=6, face="plain", colour = "#1e3648"),
  axis.title.x = element_text(size=7, face="plain", colour = "#1e3648"),
  axis.title.y = element_text(size=7, face="plain", colour = "#1e3648"),
  axis.ticks = element_line(color="#1e3648"),
  axis.line.x.bottom=element_line(color="#1e3648", size=0.5),
  axis.line.y.left=element_line(color="#1e3648", size=0.5),
  legend.key = element_blank(),
  legend.position = 'right',
  legend.title = element_text(size=7),
  legend.background = element_blank(),
  legend.text = element_text(size=7),
  legend.spacing.y = unit(0., 'cm'),
  legend.key.size = unit(0.25, 'cm'),
  legend.key.width = unit(0.5, 'cm'),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())


####  Data prep  ####

se <- function(x, na.rm = FALSE) {sd(x, na.rm=na.rm) / sqrt(sum(1*(!is.na(x))))}
zip <- function(...) {mapply(list, ..., SIMPLIFY = FALSE)}

wd = "/Users/derekvantilborg/Dropbox/PycharmProjects/traversing_chemical_space"
setwd(wd)

# compute scaffold stuff in Python
df <- read_csv("processed_results_nosmiles.csv")
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

#### colours ####

custom_colours = c("#005f73", "#94d2bd", "#0a9396", "#ee9b00", "#bbbbbb", "#f8cd48")
# custom_colours = c("#005f73", "#6f9cbc", "#6d7889", "#ee9b00", "#bbbbbb", "#f8cd48")

# custom_colours = c("#6f9cbc", "#005f73", "#ee9b00", "#73a563", "#bbbbbb", "#f8cd48")
# custom_colours = c("#73a563", "#005f73", "#6f9cbc", "#ee9b00", "#bbbbbb", "#f8cd48")

acq_funcs = c("BALD (least mutual information)", 
              "Exploit (best predictions)", 
              "Exploit, no retraining", 
              "Explore (most uncertain)", 
              "Random",
              "Similarity")
acq_cols = list(custom_colours, acq_funcs)


#### Fig 2: acq function ####

df2 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64)

df2 = df2 %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig2a = ggplot(subset(df2, dataset == 'PKM2' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
                                                      # u r b l
fig2b = ggplot(subset(df2, dataset == 'PKM2' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

fig2c = ggplot(subset(df2, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

fig2d = ggplot(subset(df2, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

fig2e = ggplot(subset(df2, dataset == 'VDR' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

fig2f = ggplot(subset(df2, dataset == 'VDR' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

# "#005f73" "#6f9cbc" "#6d7889" "#ee9b00" "#bbbbbb" "#f8cd48"


fig2 = plot_grid(fig2a, fig2b, 
                 fig2c, fig2d,
                 fig2e, fig2f,
                 labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
                 ncol=2, label_size=8)
fig2

# 180 mm/ 88 mm
dev.print(pdf, 'al_v2_fig2.pdf', width = 88/25.4, height = 123/25.4)


#### Fig 3: Explored chemistry ####

###### Ridge plots ######

df_ridge_long <- read_csv("df_ridge_long.csv")

df_ridge = subset(df_ridge_long, dataset == 'ALDH1' & hit == 1)
df_ridge$train_cycle = factor(df_ridge$train_cycle, levels=unique(df_ridge$train_cycle))
df_ = reshape2::melt(df_ridge, id.vars = c("train_cycle", "seed", "acquisition_method", 'architecture', 'bias'), variable.name = "metric", measure.vars = 9:13)
df_$acquisition_method = gsub('random', 'Random', df_$acquisition_method)
df_$acquisition_method = gsub('similarity', 'Similarity', df_$acquisition_method)
df_$acquisition_method = gsub('exploitation no retrain', 'Exploit, no retraining', df_$acquisition_method)
df_$acquisition_method = gsub('exploitation', 'Exploit (best predictions)', df_$acquisition_method)
df_$acquisition_method = gsub('bald', 'BALD (least mutual information)', df_$acquisition_method)
df_$acquisition_method = gsub('exploration', 'Explore (most uncertain)', df_$acquisition_method)
df_ = subset(df_, acquisition_method %in% c("Random", 'BALD (least mutual information)', "Exploit (best predictions)"))

df__ = subset(df_, metric == 'TPSA' & bias == 'random' & architecture == 'mlp')
df__$acquisition_method = factor(df__$acquisition_method, levels = c("Random", 'BALD (least mutual information)', "Exploit (best predictions)"))

fig3a = ggplot(subset(df__, seed == 2 )) +
  geom_density_ridges(aes(x = value, y = train_cycle, fill = acquisition_method), alpha = 0.75) + 
  labs(x='TPSA of hits', y='Active learning cycle') +
  scale_fill_manual(values=acq_cols[[1]][match(levels(df__$acquisition_method), acq_cols[[2]])]) + default_plot_theme +
  theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))

df__ = subset(df_, metric == 'MolWt' & bias == 'random' & architecture == 'mlp')
df__$acquisition_method = factor(df__$acquisition_method, levels = c("Random", 'BALD (least mutual information)', "Exploit (best predictions)"))
fig3b = ggplot(subset(df__, seed == 2 )) +
  geom_density_ridges(aes(x = value, y = train_cycle, fill = acquisition_method), alpha = 0.75) + 
  labs(x='Molecular weight of hits', y='') +
  scale_fill_manual(values=acq_cols[[1]][match(levels(df__$acquisition_method), acq_cols[[2]])]) + default_plot_theme+
  theme(plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"))


###### Line plots #####

pattern_occurence <- read_csv("pattern_occurence_ALDH1.csv", 
                              col_types = cols(...1 = col_skip(), index = col_skip()))
pattern_occurence$total_patterns = pattern_occurence$hit_patterns + pattern_occurence$nonhit_patterns

dfx = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & architecture == 'mlp' & dataset == 'ALDH1' & n_start == 64)

margins_c = unit(c(0.25, 0.25, 0, 0.25), "cm")
margins_d = unit(c(0.25, 0.25, 0, 0), "cm")
margins_e = unit(c(0.25, 0.25, 0, 0.25), "cm")
margins_f = unit(c(0.25, 0.25, 0, 0.25), "cm")
margins_g = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
margins_h = unit(c(0.25, 0.25, 0.25, 0.25), "cm")

## c ##
dfxc = dfx %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig3c = ggplot(dfxc, aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = margins_c)

## d ##
tot_hit_patterns_found_per_method = dfx %>%
  group_by(acquisition_method, seed) %>%
  summarise(across(.cols = hit_unique_patterns, list(max = max))) %>%
  ungroup()

perc_of_tot_patterns = c()
for (i_row in 1:nrow(dfx)){
  row = dfx[i_row, ]
  tot_hit_patterns = subset(tot_hit_patterns_found_per_method, acquisition_method == row$acquisition_method & seed == row$seed)$hit_unique_patterns_max
  perc = row$hit_unique_patterns/tot_hit_patterns * 100
  perc_of_tot_patterns = c(perc_of_tot_patterns, perc)
}
dfx$perc_of_tot_patterns = perc_of_tot_patterns

dfxd = dfx %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig3d = ggplot(dfxd, aes(x = total_mols_screened, y=perc_of_tot_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = perc_of_tot_patterns_mean - perc_of_tot_patterns_se, ymax = perc_of_tot_patterns_mean + perc_of_tot_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nacquired hits (% of total)', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  default_plot_theme + theme(plot.margin = margins_d)


## e ##
pattern_enrichment_aldehyde = c()
for (x in zip(dfx$total_mols_screened, dfx['Aldehyde'][[1]])){
  ntot_mols = x[[1]]
  found = x[[2]]
  # break
  expected = subset(pattern_occurence, pattern == 'Aldehyde')$total_patterns*ntot_mols/100000
  enrichment = found/expected
  pattern_enrichment_aldehyde = c(pattern_enrichment_aldehyde, enrichment)
}
dfx$pattern_enrichment_aldehyde = pattern_enrichment_aldehyde

dfxe = dfx %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = pattern_enrichment_aldehyde, 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig3e = ggplot(dfxe, aes(x = total_mols_screened, y=pattern_enrichment_aldehyde_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = pattern_enrichment_aldehyde_mean - pattern_enrichment_aldehyde_se, ymax = pattern_enrichment_aldehyde_mean + pattern_enrichment_aldehyde_se),
              color=NA, alpha=0.1) + 
  labs(y = 'Aldehydes enrichment\nin acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  default_plot_theme + theme(plot.margin = margins_e)


## f ##
expected_nr_of_aldehydes_in_bacth = subset(pattern_occurence, pattern == 'Aldehyde')$total_patterns*64/100000
expected_aldehydes_per_hit = subset(pattern_occurence, pattern == 'Aldehyde')$hit_patterns/4986

aldehyde_hit_enrichment = c()
for (x in zip(dfx$hits_discovered, dfx['hit_Aldehyde'][[1]])){
  ntot_hits = x[[1]]
  found = x[[2]]

  expected = ntot_hits*subset(pattern_occurence, pattern == 'Aldehyde')$hit_patterns/4986
  enrichment = found/expected
  aldehyde_hit_enrichment = c(aldehyde_hit_enrichment, enrichment)
}
dfx$aldehyde_hit_enrichment = aldehyde_hit_enrichment

dfxf = dfx %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = aldehyde_hit_enrichment, 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig3f = ggplot(dfxf, aes(x = total_mols_screened, y=aldehyde_hit_enrichment_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = aldehyde_hit_enrichment_mean - aldehyde_hit_enrichment_se, ymax = aldehyde_hit_enrichment_mean + aldehyde_hit_enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Aldehyde enrichment\nin acquired hits', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  # scale_y_continuous(breaks = seq(0,20, by=10)) +
  default_plot_theme + theme(plot.margin = margins_f)


## g ##
pattern_enrichment_sulfonamide = c()
for (x in zip(dfx$total_mols_screened, dfx['Sulfonamide'][[1]])){
  ntot_mols = x[[1]]
  found = x[[2]]
  expected = subset(pattern_occurence, pattern == 'Sulfonamide')$total_patterns*ntot_mols/100000
  enrichment = found/expected
  pattern_enrichment_sulfonamide = c(pattern_enrichment_sulfonamide, enrichment)
}
dfx$pattern_enrichment_sulfonamide = pattern_enrichment_sulfonamide

dfxg = dfx %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = pattern_enrichment_sulfonamide, 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig3g = ggplot(dfxg, aes(x = total_mols_screened, y=pattern_enrichment_sulfonamide_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = pattern_enrichment_sulfonamide_mean - pattern_enrichment_sulfonamide_se, ymax = pattern_enrichment_sulfonamide_mean + pattern_enrichment_sulfonamide_se),
              color=NA, alpha=0.1) +
  labs(y = 'Sulfonamide enrichment\nin acquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(ylim = c(0.5, 4), xlim = c(64, 1010), expand=F) +
  scale_y_continuous(breaks = seq(0,3, by=1)) +
  default_plot_theme + theme(plot.margin = margins_g)


## h ##
expected_nr_of_aldehydes_in_bacth = subset(pattern_occurence, pattern == 'Sulfonamide')$total_patterns*64/100000
expected_aldehydes_per_hit = subset(pattern_occurence, pattern == 'Sulfonamide')$hit_patterns/4986

sulfonamide_hit_enrichment = c()
for (x in zip(dfx$hits_discovered, dfx['hit_Sulfonamide'][[1]])){
  ntot_hits = x[[1]]
  found = x[[2]]
  
  expected = ntot_hits*subset(pattern_occurence, pattern == 'Sulfonamide')$hit_patterns/4986
  enrichment = found/expected
  sulfonamide_hit_enrichment = c(sulfonamide_hit_enrichment, enrichment)
}
dfx$sulfonamide_hit_enrichment = sulfonamide_hit_enrichment

dfxh = dfx %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = sulfonamide_hit_enrichment, 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig3h = ggplot(dfxh, aes(x = total_mols_screened, y=sulfonamide_hit_enrichment_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = sulfonamide_hit_enrichment_mean - sulfonamide_hit_enrichment_se, ymax = sulfonamide_hit_enrichment_mean + sulfonamide_hit_enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Sulfonamide enrichment\nin acquired hits', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(ylim = c(0.5, 4), xlim = c(64, 1010), expand=F) +
  scale_y_continuous(breaks = seq(0,3, by=1)) +
  default_plot_theme + theme(plot.margin = margins_h)


## grid ##
fig3ab = plot_grid(fig3a, fig3b, 
                   ncol=2, 
                   rel_widths = c(1, 1), 
                   label_size = 8, labels= c('a', 'b'))

fig3ch = plot_grid(fig3c, fig3d,
                   fig3e, fig3f,
                   fig3g, fig3h,
                   labels = c('c', 'd', 'e', 'f', 'g', 'h'), 
                   ncol=2, label_size=8)

plot_grid(fig3ab, fig3ch, 
          ncol=2, 
          rel_widths = c(1, 1))


dev.print(pdf, 'al_v2_fig3_.pdf', width = 180/25.4, height = 123/25.4)


#### Fig 4: bias ####

df4abcd = subset(df,  batch_size == 64 & architecture == 'mlp' & dataset == 'ALDH1' & n_start == 64) # train_cycle >= 0 &
df4abcd = subset(df,  batch_size == 64 & dataset == 'ALDH1' & n_start == 64) # train_cycle >= 0 &

df4abcd = df4abcd %>%
  group_by(acquisition_method, bias, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "mean_total_sims", 'mean_total_hits_sims', 
                     'n_unique_scaffolds', 'enrichment', 'mean_tani_per_batch',
                     'mean_tani_batch_to_start_batch', 'mean_tani_all_mols_to_start_batch'), 
                   list(mean = mean, se = se), na.rm = TRUE)) %>% ungroup()

df4abcd$architecture = factor(df4abcd$architecture, levels=c('mlp', 'gcn'))

bias_colours = c("#94d2bd", "#0a9396", "#005f73")

fig4a = ggplot(subset(df4abcd, acquisition_method %in% c("Similarity")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = 'Similarity of acquired\nbatch to start data', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.45, by=0.05)) +
  default_plot_theme + theme(plot.margin = unit(c(0.5, 0, 0, 0.25), "cm"))

fig4b = ggplot(subset(df4abcd, acquisition_method %in% c("BALD (least mutual information)")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.45, by=0.05)) +
  default_plot_theme + theme(plot.margin = unit(c(0.5, 0.25, 0, 0.25), "cm"))

fig4c = ggplot(subset(df4abcd, acquisition_method %in% c("Exploit (best predictions)")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = 'Similarity of acquired\nbatch to start data', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.45, by=0.05)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.5, 0.25), "cm"))

fig4d = ggplot(subset(df4abcd, acquisition_method %in% c("Exploit, no retraining")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.45, by=0.05)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.5, 0.25), "cm"))


plot_grid(fig4a, fig4b, fig4c, fig4d,
          ncol=2, 
          label_size = 8, labels= c('a', 'b', 'c', 'd'))


dev.print(pdf, 'al_v2_fig4_v2.pdf', width = 88/25.4, height = 82/25.4)



#### Fig 5: start size ####


df5 = subset(df, total_mols_screened == 1000 & batch_size == 64 & bias == 'Innate')
df5 = subset(df5, acquisition_method %in% c("Random", "Similarity", "BALD (least mutual information)"))

df5 = df5 %>%
  group_by(acquisition_method, n_start, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()


fig5a = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'mlp'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
# u r b l
fig5b = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'gcn'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

fig5c = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

fig5d = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

fig5e = ggplot(subset(df5, dataset == 'VDR' & architecture == 'mlp'),  aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='molecules in the start set', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

fig5f = ggplot(subset(df5, dataset == 'VDR' & architecture == 'gcn'),  aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = '', x='molecules in the start set', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

# "#005f73" "#6f9cbc" "#6d7889" "#ee9b00" "#bbbbbb" "#f8cd48"

fig5 = plot_grid(fig5a, fig5b, 
                 fig5c, fig5d,
                 fig5e, fig5f,
                 labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
                 ncol=2, label_size=8)
fig5

# 180 mm/ 88 mm
dev.print(pdf, 'al_v2_fig5_v1.pdf', width = 88/25.4, height = 123/25.4)



####

df6 = subset(df, train_cycle <= 15 & batch_size == 64 & bias == 'Innate')
df6 = subset(df6, acquisition_method %in% c("Similarity", "BALD (least mutual information)"))

df6 = df6[c(2, 32, 34, 35, 36, 38, 41, 151)]
df6$id = paste(df6$architecture, df6$n_start, df6$seed, df6$dataset)

df6_sim = subset(df6, acquisition_method == "Similarity")
df6_bald = subset(df6, acquisition_method == "BALD (least mutual information)")

df6_sim = df6_sim[order(df6_sim$id),]
df6_bald = df6_bald[order(df6_bald$id),]

# difference in enrichment
df6_sim$enrichment = df6_bald$enrichment - df6_sim$enrichment 


df6 = df6_sim %>%
  group_by(acquisition_method, n_start, dataset, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

col_gradient = c("#bbbbbb", "#94d2bd", "#0a9396", "#005F73", "#f8cd48", "#ee9b00")


fig6a = ggplot(subset(df6, dataset == 'PKM2' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = 'Increase in enrichment\nover similarity search', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

fig6b = ggplot(subset(df6, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = '', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

fig6c = ggplot(subset(df6, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = 'Increase in enrichment\nover similarity search', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

fig6d = ggplot(subset(df6, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = '', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

fig6e = ggplot(subset(df6, dataset == 'VDR' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = 'Increase in enrichment\nover similarity search', x='active learning cycle', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

fig6f = ggplot(subset(df6, dataset == 'VDR' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = '', x='active learning cycle', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

fig6 = plot_grid(fig6a, fig6b, 
                 fig6c, fig6d,
                 fig6e, fig6f,
                 labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
                 ncol=2, label_size=8)
fig6

dev.print(pdf, 'al_v2_fig6_v2.pdf', width = 88/25.4, height = 123/25.4)


### Supp. Fig. 1: UMAP ####

df_UMAP_PKM2 <- read_csv("UMAP_PKM2.csv")
df_UMAP_PKM2 = df_UMAP_PKM2[order(df_UMAP_PKM2$y, decreasing=F),]
df_UMAP_PKM2$y = factor(df_UMAP_PKM2$y)
df_UMAP_PKM2$y_alpha = as.numeric(df_UMAP_PKM2$y )
df_UMAP_PKM2$y_size = df_UMAP_PKM2$y_alpha
df_UMAP_PKM2$y_size[df_UMAP_PKM2$y_size == 2] = 0.25
df_UMAP_ALDH1$y_alpha[df_UMAP_PKM2$y_alpha == 1] = 0.05

sfig1a = ggplot(data=df_UMAP_PKM2, aes(x=UMAP1, y=UMAP2, color=y))+
  geom_point(size=df_UMAP_PKM2$y_size, shape=16, alpha=df_UMAP_PKM2$y_alpha) +
  labs(title='PKM2')+
  scale_color_manual(values=c("#c3d6d9", "#005F73")) +
  default_plot_theme

df_UMAP_ALDH1 <- read_csv("UMAP_ALDH1.csv")
df_UMAP_ALDH1 = df_UMAP_ALDH1[order(df_UMAP_ALDH1$y, decreasing=F),]
df_UMAP_ALDH1$y = factor(df_UMAP_ALDH1$y)
df_UMAP_ALDH1$y_alpha = as.numeric(df_UMAP_ALDH1$y )
df_UMAP_ALDH1$y_size = df_UMAP_ALDH1$y_alpha
df_UMAP_ALDH1$y_size[df_UMAP_ALDH1$y_size == 2] = 0.25
df_UMAP_ALDH1$y_alpha[df_UMAP_ALDH1$y_alpha == 1] = 0.05

sfig1b = ggplot(data=df_UMAP_ALDH1, aes(x=UMAP1, y=UMAP2, color=y))+
  geom_point(size=df_UMAP_ALDH1$y_size, shape=16, alpha=df_UMAP_ALDH1$y_alpha) +
  labs(title='ALDH1')+
  scale_color_manual(values=c("#c3d6d9", "#005F73")) +
  default_plot_theme

df_UMAP_VDR <- read_csv("UMAP_VDR.csv")
df_UMAP_VDR = df_UMAP_VDR[order(df_UMAP_VDR$y, decreasing=F),]
df_UMAP_VDR$y = factor(df_UMAP_VDR$y)
df_UMAP_VDR$y_alpha = as.numeric(df_UMAP_VDR$y )
df_UMAP_VDR$y_size = df_UMAP_VDR$y_alpha
df_UMAP_VDR$y_size[df_UMAP_VDR$y_size == 2] = 0.25
df_UMAP_VDR$y_alpha[df_UMAP_VDR$y_alpha == 1] = 0.05

sfig1c = ggplot(data=df_UMAP_VDR, aes(x=UMAP1, y=UMAP2, color=y))+
  geom_point(size=df_UMAP_VDR$y_size, shape=16, alpha=df_UMAP_VDR$y_alpha) +
  labs(title='VDR')+
  scale_color_manual(values=c("#c3d6d9", "#005F73")) +
  default_plot_theme

plot_grid(sfig1a, sfig1b, sfig1c, 
          labels = c('a', 'b', 'c'), 
          ncol=3, label_size=8)

dev.print(pdf, 'al_sfig1_umap.pdf', width = 180/25.4, height = 60/25.4)


#### PCA ###

df_pca <- read_csv("pca.csv")

df_pca = df_pca[order(df_pca$exploitation_static, decreasing=F),]
df_pca$exploitation_static = factor(df_pca$exploitation_static)
df_pca$exploitation_static_alpha = as.numeric(df_pca$exploitation_static )
df_pca$exploitation_static_alpha[df_pca$exploitation_static_alpha == 1] = 0.05
p1 = ggplot(data=df_pca, aes(x=UMAP1, y=UMAP2, color=exploitation_static))+
  geom_point(size=0.75, shape=16, alpha=df_pca$exploitation_static_alpha) +
  labs(title='exploitation_static')+
  scale_color_manual(values=c("#c3d6d9", "#005F73", "#ee9b00")) +
  default_plot_theme

df_pca = df_pca[order(df_pca$similarity, decreasing=F),]
df_pca$similarity = factor(df_pca$similarity)
df_pca$similarity_alpha = as.numeric(df_pca$similarity )
df_pca$similarity_alpha[df_pca$similarity_alpha == 1] = 0.05
p2 = ggplot(data=df_pca, aes(x=UMAP1, y=UMAP2, color=similarity))+
  geom_point(size=0.75, shape=16, alpha=df_pca$similarity_alpha) +
  labs(title='similarity')+
  scale_color_manual(values=c("#c3d6d9", "#005F73", "#ee9b00")) +
  default_plot_theme

df_pca = df_pca[order(df_pca$exploitation, decreasing=F),]
df_pca$exploitation = factor(df_pca$exploitation)
df_pca$exploitation_alpha = as.numeric(df_pca$exploitation )
df_pca$exploitation_alpha[df_pca$exploitation_alpha == 1] = 0.05
p3 = ggplot(data=df_pca, aes(x=UMAP1, y=UMAP2, color=exploitation))+
  geom_point(size=0.75, shape=16, alpha=df_pca$exploitation_alpha) +
  labs(title='exploitation')+
  scale_color_manual(values=c("#c3d6d9", "#005F73", "#ee9b00")) +
  default_plot_theme

plot_grid(p2, p1, p3, 
          labels = c('a', 'b', 'c'), 
          ncol=3, label_size=8)

dev.print(pdf, 'al_v2_fig4_umap_.pdf', width = 180/25.4, height = 60/25.4)


#### Supplementary ####
##### S1. UMAP ######

umap_aldh1 <- read_csv("Dropbox/PycharmProjects/Active_Learning_Simulation/data/ALDH1/screen/umap.csv")
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

umap_pkm2 <- read_csv("Dropbox/PycharmProjects/Active_Learning_Simulation/data/PKM2/screen/umap.csv")
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

umap_vdr <- read_csv("Dropbox/PycharmProjects/Active_Learning_Simulation/data/VDR/screen/umap.csv")
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


fig1 = plot_grid(p_umap_pkm2, p_umap_aldh1, p_umap_vdr, rel_widths = c(2, 2, 2),
                 labels = c('a', 'b', 'c'),
                 ncol=3, label_size=8)
fig1


dev.print(pdf, 'al_Sup_fig1.pdf', width = 180/25.4, height = 60/25.4)


##### S2. ROC AUC #####

dfs2 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64) %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

figs2a = ggplot(subset(dfs2, dataset == 'PKM2' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean ROC AUC\non the test set', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

figs2b = ggplot(subset(dfs2, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

figs2c = ggplot(subset(dfs2, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean ROC AUC\non the test set', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

figs2d = ggplot(subset(dfs2, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

figs2e = ggplot(subset(dfs2, dataset == 'VDR' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean ROC AUC\non the test set', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

figs2f = ggplot(subset(dfs2, dataset == 'VDR' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

figs2 = plot_grid(figs2a, figs2b, 
                  figs2c, figs2d,
                  figs2e, figs2f,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
                  ncol=2, label_size=8)
figs2

# 180 mm/ 88 mm
dev.print(pdf, 'al_Sup_fig2.pdf', width = 88/25.4, height = 110/25.4)


##### S3. Batch size #####

dfs3a = subset(df, train_cycle == 15 & architecture == 'mlp' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs3a$batch_size = as.numeric(as.character(dfs3a$batch_size))

figs3a = ggplot(dfs3a, aes(x = batch_size, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0, 6.75), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.5), "cm"))


dfs3b = subset(df, train_cycle == 15 & architecture == 'gcn' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs3b$batch_size = as.numeric(as.character(dfs3b$batch_size))

figs3b = ggplot(dfs3b, aes(x = batch_size, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
  labs(y = '', x='') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0, 6.75), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"))


dfs3c = subset(df, train_cycle == 15 & architecture == 'mlp' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs3c$batch_size = as.numeric(as.character(dfs3c$batch_size))

figs3c = ggplot(dfs3c, aes(x = batch_size, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
  labs(y = 'Test set ROC AUC after\n screening 1000 molecules', x='Molecules acquired per cycle') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0.45, 0.65), expand=F) +
  scale_y_continuous(breaks = seq(0.45,0.65, by=0.05)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))


dfs3d = subset(df, train_cycle == 15 & architecture == 'gcn' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs3d$batch_size = as.numeric(as.character(dfs3d$batch_size))

figs3d = ggplot(dfs3d, aes(x = batch_size, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
  labs(y = '', x='Molecules acquired per cycle') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0.45, 0.65), expand=F) +
  scale_y_continuous(breaks = seq(0.45,0.65, by=0.05)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


plot_grid(figs3a, figs3b, figs3c, figs3d,
          ncol=2, 
          label_size = 8, labels= c('a', 'b', 'c', 'd'))

dev.print(pdf, 'al_Sup_fig3.pdf', width = 88/25.4, height = 80/25.4)


##### S4. Bias ######
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.5), "cm"))

figs4c = ggplot(dfs4ac, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
  labs(y = '', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"))

figs4d = ggplot(dfs4bd, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.5), "cm"))

figs4g = ggplot(dfs4eg, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
  labs(y = '', x='') +
  coord_cartesian(xlim = c(0.1, 0.53), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"))

figs4h = ggplot(dfs4fh, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='Structural bias in start data') +
  coord_cartesian(xlim = c(0.12, 0.4), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7.5, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.5), "cm"))

figs4k = ggplot(dfs4ik, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
  labs(y = '', x='Structural bias in start data') +
  coord_cartesian(xlim = c(0.12, 0.4), ylim = c(0, 7.5), expand=F) + 
  scale_y_continuous(breaks = seq(0,7.5, by=1)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"))

figs4l = ggplot(dfs4jl, aes(x = starting_bias, y = test_roc_auc, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
  labs(y = '', x='Structural bias in start data') +
  coord_cartesian(xlim = c(0.12, 0.4), ylim = c(0.35, 0.75), expand=F) + 
  scale_y_continuous(breaks = seq(0.35,0.75, by=0.05)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


###
plot_grid(figs4a, figs4b, figs4c,  figs4d, 
          figs4e, figs4f, figs4g,  figs4h, 
          figs4i, figs4j, figs4k,  figs4l, 
          ncol=4, 
          label_size = 8, labels= c('a', 'b', 'c', 'd',
                                    'e', 'f', 'g', 'h', 
                                    'i', 'j', 'k', 'l'))

dev.print(pdf, 'al_Sup_fig4.pdf', width = 180/25.4, height = 130/25.4)


##### S5. ####

dfx = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & architecture == 'mlp' & dataset == 'ALDH1' & n_start == 64)

dfxx = dfx %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

figsx = ggplot(dfxx, aes(x = total_mols_screened, y=hit_unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = hit_unique_patterns_mean - hit_unique_patterns_se, ymax = hit_unique_patterns_mean + hit_unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique substructures in\nacquired hits', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  # scale_y_continuous(breaks = seq(15,45, by=5)) +
  default_plot_theme + theme(plot.margin = margins_d)

dev.print(pdf, 'al_Sup_fig5.pdf', width = 44/25.4, height = 44/25.4)

#### End ####


## Video plot

# # # df$train_cycle
# cycle = 15
# df5 = subset(df, train_cycle == cycle & batch_size == 64 & bias == 'Innate')
# df5 = subset(df5, acquisition_method %in% c("Random", "Similarity", "BALD (least mutual information)"))
# 
# df5 = df5 %>%
#   group_by(acquisition_method, n_start, dataset, total_mols_screened, train_cycle, architecture) %>%
#   summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
#                      'test_balanced_accuracy', 'enrichment'), 
#                    list(mean = mean, sd = sd, se = se))) %>% ungroup()
# 
# 
# fig5a = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'mlp'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
#   geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
#   geom_point(size=1) +
#   labs(y = paste('Mean enrichment after', cycle, 'cycles'), x='', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0, 6), expand=F) +
#   scale_y_continuous(breaks = seq(0,6, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
# # u r b l
# fig5b = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'gcn'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
#   geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
#   geom_point(size=1) +
#   labs(y = '', x='', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0, 6), expand=F) +
#   scale_y_continuous(breaks = seq(0,6, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))
# 
# fig5c = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
#   geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
#   geom_point(size=1) +
#   labs(y = paste('Mean enrichment after', cycle, 'cycles'), x='', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0, 6), expand=F) +
#   scale_y_continuous(breaks = seq(0,6, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))
# 
# fig5d = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
#   geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
#   geom_point(size=1) +
#   labs(y = '', x='', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0, 6), expand=F) +
#   scale_y_continuous(breaks = seq(0,6, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))
# 
# fig5e = ggplot(subset(df5, dataset == 'VDR' & architecture == 'mlp'),  aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
#   geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
#   geom_point(size=1) +
#   labs(y = paste('Mean enrichment after', cycle, 'cycles'), x='molecules in the start set', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0, 6), expand=F) +
#   scale_y_continuous(breaks = seq(0,6, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))
# 
# fig5f = ggplot(subset(df5, dataset == 'VDR' & architecture == 'gcn'),  aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
#   geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
#   geom_point(size=1) +
#   labs(y = '', x='molecules in the start set', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
#   coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0, 6), expand=F) +
#   scale_y_continuous(breaks = seq(0,6, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))
# 
# # "#005f73" "#6f9cbc" "#6d7889" "#ee9b00" "#bbbbbb" "#f8cd48"
# 
# fig5 = plot_grid(fig5a, fig5b, 
#                  fig5c, fig5d,
#                  fig5e, fig5f,
#                  labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
#                  ncol=2, label_size=8)
# fig5
