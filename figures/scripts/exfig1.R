

df5 <- read_csv("figures/data/exfig1.csv")

df5$n_start = factor(df5$n_start, levels=sort(unique(df5$n_start)))


df5 = df5 %>%
  group_by(acquisition_method, n_start, dataset, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()


fig5a = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) +
  labs(y = 'Increase in enrichment\nover similarity search', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

fig5b = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) +
  labs(y = '', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

fig5b = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) +
  labs(y = '', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

fig5c = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) + 
  labs(y = 'Increase in enrichment\nover similarity search', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

fig5d = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) + 
  labs(y = '', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

fig5e = ggplot(subset(df5, dataset == 'VDR' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) + 
  labs(y = 'Increase in enrichment\nover similarity search', x='active learning cycle', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

fig5f = ggplot(subset(df5, dataset == 'VDR' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) + 
  labs(y = '', x='active learning cycle', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

fig5 = plot_grid(fig5a, fig5b, 
                 fig5c, fig5d,
                 fig5e, fig5f,
                 labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
                 ncol=2, label_size=8)

fig5 = plot_grid(plot_spacer(), fig5, plot_spacer(),
                 labels = c('', '', ''),
                 ncol=3, rel_widths = c(1, 2, 1))



pdf('figures/exfig1.pdf', width = 180/25.4, height = 123/25.4)
print(fig5)
dev.off()
