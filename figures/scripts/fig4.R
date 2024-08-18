


df4_umap <- read_csv("figures/data/fig4_abc.csv")

df4_umap = df4_umap[order(df4_umap$exploitation, decreasing=F),]
df4_umap$exploitation = factor(df4_umap$exploitation)
df4_umap$exploitation_alpha = as.numeric(df4_umap$exploitation )
df4_umap$exploitation_alpha[df4_umap$exploitation_alpha == 1] = 0.05
fig4a = ggplot(data=df4_umap, aes(x=UMAP1, y=UMAP2, color=exploitation))+
  geom_point(size=0.75, shape=16, alpha=df4_umap$exploitation_alpha) +
  labs(title='exploitation')+
  scale_color_manual(values=umap_cols) +
  default_plot_theme

table(df4_umap$exploitation_alpha)

df4_umap = df4_umap[order(df4_umap$exploitation_static, decreasing=F),]
df4_umap$exploitation_static = factor(df4_umap$exploitation_static)
df4_umap$exploitation_static_alpha = as.numeric(df4_umap$exploitation_static )
df4_umap$exploitation_static_alpha[df4_umap$exploitation_static_alpha == 1] = 0.05
fig4b = ggplot(data=df4_umap, aes(x=UMAP1, y=UMAP2, color=exploitation_static))+
  geom_point(size=0.75, shape=16, alpha=df4_umap$exploitation_static_alpha) +
  labs(title='exploitation_static')+
  scale_color_manual(values=umap_cols) +
  default_plot_theme

df4_umap = df4_umap[order(df4_umap$similarity, decreasing=F),]
df4_umap$similarity = factor(df4_umap$similarity)
df4_umap$similarity_alpha = as.numeric(df4_umap$similarity )
df4_umap$similarity_alpha[df4_umap$similarity_alpha == 1] = 0.05
fig4c = ggplot(data=df4_umap, aes(x=UMAP1, y=UMAP2, color=similarity))+
  geom_point(size=0.75, shape=16, alpha=df4_umap$similarity_alpha) +
  labs(title='similarity')+
  scale_color_manual(values=umap_cols) +
  default_plot_theme


fig4abc = plot_grid(fig4a, fig4b, fig4c,
                    labels = c('a', 'b', 'c'),
                    ncol=3, label_size=8)

pdf('figures/fig4_umaps_.pdf', width = 180/25.4, height = 60/25.4)
print(fig4abc)
dev.off()

######

df4_ridge <- read_csv("figures/data/fig4_de.csv")

df4_ridge$acquisition_method = factor(df4_ridge$acquisition_method, levels = c("Random", 'BALD (least mutual information)', "Exploit (best predictions)"))
df4_ridge$train_cycle = factor(df4_ridge$train_cycle, levels=unique(df4_ridge$train_cycle))

fig4d = ggplot(subset(df4_ridge, metric == 'TPSA')) +
  geom_density_ridges(aes(x = value, y = train_cycle, fill = acquisition_method), alpha = 0.75, linewidth=0.35) + 
  labs(x='TPSA of hits', y='Active learning cycle') +
  scale_fill_manual(values=acq_cols[[1]][match(levels(df4_ridge$acquisition_method), acq_cols[[2]])]) + default_plot_theme +
  theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))


fig4e = ggplot(subset(df4_ridge, metric == 'MolWt')) +
  geom_density_ridges(aes(x = value, y = train_cycle, fill = acquisition_method), alpha = 0.75, linewidth=0.35) + 
  labs(x='Molecular weight of hits', y='') +
  scale_fill_manual(values=acq_cols[[1]][match(levels(df4_ridge$acquisition_method), acq_cols[[2]])]) + default_plot_theme+
  theme(plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"))

######

df4_fghijk <- read_csv("figures/data/fig4_fghijk.csv")

## f ##

df4f = df4_fghijk %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric),
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4f = ggplot(df4f, aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) +
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) +
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

## g ##

df4g = df4_fghijk %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric),
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4g = ggplot(df4g, aes(x = total_mols_screened, y=perc_of_tot_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = perc_of_tot_patterns_mean - perc_of_tot_patterns_se, ymax = perc_of_tot_patterns_mean + perc_of_tot_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nacquired hits (% of total)', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) +
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))

## h ##

df4h = df4_fghijk %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = pattern_enrichment_aldehyde,
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4h = ggplot(df4h, aes(x = total_mols_screened, y=pattern_enrichment_aldehyde_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = pattern_enrichment_aldehyde_mean - pattern_enrichment_aldehyde_se, ymax = pattern_enrichment_aldehyde_mean + pattern_enrichment_aldehyde_se),
              color=NA, alpha=0.1) +
  labs(y = 'Aldehydes enrichment\nin acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) +
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))


## i ## 

df4i = df4_fghijk %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = aldehyde_hit_enrichment,
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4i = ggplot(df4i, aes(x = total_mols_screened, y=aldehyde_hit_enrichment_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = aldehyde_hit_enrichment_mean - aldehyde_hit_enrichment_se, ymax = aldehyde_hit_enrichment_mean + aldehyde_hit_enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Aldehyde enrichment\nin acquired hits', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) +
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  # scale_y_continuous(breaks = seq(0,20, by=10)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

## j ## 

df4j = df4_fghijk %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = pattern_enrichment_sulfonamide,
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4j = ggplot(df4j, aes(x = total_mols_screened, y=pattern_enrichment_sulfonamide_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = pattern_enrichment_sulfonamide_mean - pattern_enrichment_sulfonamide_se, ymax = pattern_enrichment_sulfonamide_mean + pattern_enrichment_sulfonamide_se),
              color=NA, alpha=0.1) +
  labs(y = 'Sulfonamide enrichment\nin acquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) +
  coord_cartesian(ylim = c(0.5, 4), xlim = c(64, 1010), expand=F) +
  scale_y_continuous(breaks = seq(0,3, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))

## k ## 

df4k = df4_fghijk %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = sulfonamide_hit_enrichment,
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4k = ggplot(df4k, aes(x = total_mols_screened, y=sulfonamide_hit_enrichment_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = sulfonamide_hit_enrichment_mean - sulfonamide_hit_enrichment_se, ymax = sulfonamide_hit_enrichment_mean + sulfonamide_hit_enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Sulfonamide enrichment\nin acquired hits', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.45) +
  coord_cartesian(ylim = c(0.5, 4), xlim = c(64, 1010), expand=F) +
  scale_y_continuous(breaks = seq(0,3, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


cor.test(x,y, method="kendall")

## plotting ##

## grid ##
fig4de = plot_grid(fig4d, fig4e, 
                   ncol=2, 
                   rel_widths = c(1, 1), 
                   label_size = 8, labels= c('a', 'b'))

fig4fk = plot_grid(fig4f, fig4g,
                   fig4h, fig4i,
                   fig4j, fig4k,
                   labels = c('c', 'd', 'e', 'f', 'g', 'h'), 
                   ncol=2, label_size=8)

fig4 = plot_grid(fig4de, fig4fk, ncol=2, rel_widths = c(1, 1))

pdf('figures/fig4.pdf', width = 180/25.4, height = 123/25.4)
print(fig4)
dev.off()


