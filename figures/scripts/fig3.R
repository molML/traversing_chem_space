
df3 <- read_csv('figures/data/fig3.csv')


###### plots ######

df3 = df3 %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig3a = ggplot(subset(df3, dataset == 'PKM2' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
# u r b l
fig3b = ggplot(subset(df3, dataset == 'PKM2' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

fig3c = ggplot(subset(df3, dataset == 'PKM2' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

fig3d = ggplot(subset(df3, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

fig3e = ggplot(subset(df3, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

fig3f = ggplot(subset(df3, dataset == 'ALDH1' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

fig3g = ggplot(subset(df3, dataset == 'VDR' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

fig3h = ggplot(subset(df3, dataset == 'VDR' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

fig3i = ggplot(subset(df3, dataset == 'VDR' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

# "#005f73" "#6f9cbc" "#6d7889" "#ee9b00" "#bbbbbb" "#f8cd48"


fig3 = plot_grid(fig3a, fig3b, fig3c, 
                 fig3d, fig3e, fig3f,
                 fig3g, fig3h, fig3i,
                 labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' ,'i'), 
                 ncol=3, label_size=8)

# 180 mm/ 88 mm
pdf('figures/fig3.pdf', width = 120/25.4, height = 105/25.4)
print(fig3)
dev.off()




