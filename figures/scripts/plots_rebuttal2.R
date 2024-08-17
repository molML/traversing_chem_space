

#### General ####

library(readr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(data.table)

##### plot themes #####

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


#####  Data prep  #####

se <- function(x, na.rm = FALSE) {sd(x, na.rm=na.rm) / sqrt(sum(1*(!is.na(x))))}
zip <- function(...) {mapply(list, ..., SIMPLIFY = FALSE)}

wd = "/Users/derekvantilborg/Dropbox/PycharmProjects/traversing_chemical_space"
# setwd(wd)

# compute scaffold stuff in Python
df <- read_csv("processed_results_nosmiles_rebuttal.csv")
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


##### colours #####

custom_colours = c("#005f73", "#94d2bd", "#0a9396", "#ee9b00", "#bbbbbb", "#f8cd48")

acq_funcs = c("BALD (least mutual information)", 
              "Exploit (best predictions)", 
              "Exploit, no retraining", 
              "Explore (most uncertain)", 
              "Random",
              "Similarity")
acq_funcs2 = c("BALD (least mutual information)", 
              "Exploit (best predictions)", 
              "Exploit, no retraining", 
              "Random",
              "Similarity")
acq_cols = list(custom_colours, acq_funcs)



#### Fig 2: Structural diversity ####

df2 = subset(df,  batch_size == 64 & dataset == 'ALDH1' & n_start == 64) # train_cycle >= 0 &

df2 = df2 %>%
  group_by(acquisition_method, bias, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "mean_total_sims", 'mean_total_hits_sims', 
                     'n_unique_scaffolds', 'enrichment', 'mean_tani_per_batch',
                     'mean_tani_batch_to_start_batch', 'mean_tani_all_mols_to_start_batch'), 
                   list(mean = mean, se = se), na.rm = TRUE)) %>% ungroup()

df2$architecture = factor(df2$architecture, levels=c('mlp', 'gcn'))

bias_colours = c("#94d2bd", "#0a9396", "#005f73")

fig2a = ggplot(subset(df2, acquisition_method %in% c("Similarity")), 
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

fig2b = ggplot(subset(df2, acquisition_method %in% c("BALD (least mutual information)")), 
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

fig2c = ggplot(subset(df2, acquisition_method %in% c("Exploit (best predictions)")), 
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

fig2d = ggplot(subset(df2, acquisition_method %in% c("Exploit, no retraining")), 
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


fig2 = plot_grid(fig2a, fig2b, fig2c, fig2d,
          ncol=2, 
          label_size = 8, labels= c('a', 'b', 'c', 'd'))
fig2

# dev.print(pdf, 'al_v4_fig2.pdf', width = 88/25.4, height = 82/25.4)
rm(fig2a, fig2b, fig2c, fig2d)



#### Fig 3: acq function ####

df3 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64)

###### Statistics Acq functions #######

acq_last = subset(df3, total_mols_screened == 1000)
acq_statistics = list()
for (dataset_ in unique(acq_last$dataset)){      
  for (algo in c('mlp', 'gcn', 'rf')){
    for (acq1 in acq_funcs){
      acq_statistics[['datasets']] = c(acq_statistics[['datasets']] , dataset_)
      acq_statistics[['algos']] = c(acq_statistics[['algos']] , algo)
      acq_statistics[['acqs']] = c(acq_statistics[['acqs']] , acq1)
      
      a_ = subset(acq_last, dataset == dataset_ & architecture == algo & acquisition_method == acq1)$enrichment
      b_ = subset(acq_last, dataset == dataset_ & architecture == algo & acquisition_method == "Random")$enrichment
      
      if (length(a_) + length(b_) > 0 & length(a_) == length(b_)){
        
        p = wilcox.test(a_, b_, paired=T)$p.value
        
        if (is.na(p)){
          p=1
        }
        if (p < 0.05){
          acq_statistics[["wilcox"]] = c(acq_statistics[["wilcox"]], '*')
        } else {
          acq_statistics[["wilcox"]] = c(acq_statistics[["wilcox"]], '')
        }
      } else {
        acq_statistics[["wilcox"]] = c(acq_statistics[["wilcox"]], '')
      }
    }
  }
}
acq_statistics = data.frame(acq_statistics)
rm(acq_last)


# # acq_funcs
# a_ = subset(acq_last, dataset == 'VDR' & architecture == 'mlp' & acquisition_method == "BALD (least mutual information)")$enrichment
# b_ = subset(acq_last, dataset == 'VDR' & architecture == 'gcn' & acquisition_method == "BALD (least mutual information)")$enrichment
# p = wilcox.test(a_, b_, paired=T)$p.value
# if (p<0.05){print('*')}
# 
# 

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
  geom_line(size=0.6) + 
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
  geom_line(size=0.6) + 
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
  geom_line(size=0.6) + 
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
  geom_line(size=0.6) + 
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
  geom_line(size=0.6) + 
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
  geom_line(size=0.6) + 
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
  geom_line(size=0.6) + 
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
  geom_line(size=0.6) + 
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
  geom_line(size=0.6) + 
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
fig3

# 180 mm/ 88 mm
dev.print(pdf, 'al_v4_fig3.pdf', width = 110/25.4, height = 100/25.4)
rm(fig3a, fig3b, fig3c, fig3d, fig3e, fig3f)


#### Fig 4: Explored chemistry ####

###### U-Maps ######

df4_umap <- read_csv("pca.csv")

df4_umap = df4_umap[order(df4_umap$exploitation, decreasing=F),]
df4_umap$exploitation = factor(df4_umap$exploitation)
df4_umap$exploitation_alpha = as.numeric(df4_umap$exploitation )
df4_umap$exploitation_alpha[df4_umap$exploitation_alpha == 1] = 0.05
fig4a = ggplot(data=df4_umap, aes(x=UMAP1, y=UMAP2, color=exploitation))+
  geom_point(size=0.75, shape=16, alpha=df4_umap$exploitation_alpha) +
  labs(title='exploitation')+
  scale_color_manual(values=c("#c3d6d9", "#005F73", "#ee9b00")) +
  default_plot_theme

df4_umap = df4_umap[order(df4_umap$exploitation_static, decreasing=F),]
df4_umap$exploitation_static = factor(df4_umap$exploitation_static)
df4_umap$exploitation_static_alpha = as.numeric(df4_umap$exploitation_static )
df4_umap$exploitation_static_alpha[df4_umap$exploitation_static_alpha == 1] = 0.05
fig4b = ggplot(data=df4_umap, aes(x=UMAP1, y=UMAP2, color=exploitation_static))+
  geom_point(size=0.75, shape=16, alpha=df4_umap$exploitation_static_alpha) +
  labs(title='exploitation_static')+
  scale_color_manual(values=c("#c3d6d9", "#005F73", "#ee9b00")) +
  default_plot_theme

df4_umap = df4_umap[order(df4_umap$similarity, decreasing=F),]
df4_umap$similarity = factor(df4_umap$similarity)
df4_umap$similarity_alpha = as.numeric(df4_umap$similarity )
df4_umap$similarity_alpha[df4_umap$similarity_alpha == 1] = 0.05
fig4c = ggplot(data=df4_umap, aes(x=UMAP1, y=UMAP2, color=similarity))+
  geom_point(size=0.75, shape=16, alpha=df4_umap$similarity_alpha) +
  labs(title='similarity')+
  scale_color_manual(values=c("#c3d6d9", "#005F73", "#ee9b00")) +
  default_plot_theme


fig4abc = plot_grid(fig4a, fig4b, fig4c, 
          labels = c('a', 'b', 'c'), 
          ncol=3, label_size=8)
fig4abc

dev.print(pdf, 'al_v4_fig4_umaps.pdf', width = 180/25.4, height = 60/25.4)
rm(fig4a, fig4b, fig4c)

###### Ridge plots ######


df4_ridge = subset(read_csv("df_ridge_long.csv"), dataset == 'ALDH1' & hit == 1)
df4_ridge$train_cycle = factor(df4_ridge$train_cycle, levels=unique(df4_ridge$train_cycle))
df4_ridge = reshape2::melt(df4_ridge, id.vars = c("train_cycle", "seed", "acquisition_method", 'architecture', 'bias'), variable.name = "metric", measure.vars = 9:13)
df4_ridge$acquisition_method = gsub('random', 'Random', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('similarity', 'Similarity', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('exploitation no retrain', 'Exploit, no retraining', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('exploitation', 'Exploit (best predictions)', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('bald', 'BALD (least mutual information)', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('exploration', 'Explore (most uncertain)', df4_ridge$acquisition_method)
df4_ridge = subset(df4_ridge, acquisition_method %in% c("Random", 'BALD (least mutual information)', "Exploit (best predictions)"))

df4_ridge_d = subset(df4_ridge, metric == 'TPSA' & bias == 'random' & architecture == 'mlp')
df4_ridge_d$acquisition_method = factor(df4_ridge_d$acquisition_method, levels = c("Random", 'BALD (least mutual information)', "Exploit (best predictions)"))

fig4d = ggplot(subset(df4_ridge_d, seed == 2 )) +
  geom_density_ridges(aes(x = value, y = train_cycle, fill = acquisition_method), alpha = 0.75) + 
  labs(x='TPSA of hits', y='Active learning cycle') +
  scale_fill_manual(values=acq_cols[[1]][match(levels(df4_ridge_d$acquisition_method), acq_cols[[2]])]) + default_plot_theme +
  theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))

df4_ridge_e = subset(df4_ridge, metric == 'MolWt' & bias == 'random' & architecture == 'mlp')
df4_ridge_e$acquisition_method = factor(df4_ridge_e$acquisition_method, levels = c("Random", 'BALD (least mutual information)', "Exploit (best predictions)"))
fig4e = ggplot(subset(df4_ridge_e, seed == 2 )) +
  geom_density_ridges(aes(x = value, y = train_cycle, fill = acquisition_method), alpha = 0.75) + 
  labs(x='Molecular weight of hits', y='') +
  scale_fill_manual(values=acq_cols[[1]][match(levels(df4_ridge_e$acquisition_method), acq_cols[[2]])]) + default_plot_theme+
  theme(plot.margin = unit(c(0.25, 0, 0.25, 0), "cm"))

remove(df4_ridge)


###### Line plots #####

# margins_h = unit(c(0.25, 0.25, 0, 0.25), "cm")
# margins_i = unit(c(0.25, 0.25, 0, 0.25), "cm")
# margins_j = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
# margins_k = unit(c(0.25, 0.25, 0.25, 0.25), "cm")

pattern_occurence <- read_csv("pattern_occurence_ALDH1.csv", 
                              col_types = cols(...1 = col_skip(), index = col_skip()))
pattern_occurence$total_patterns = pattern_occurence$hit_patterns + pattern_occurence$nonhit_patterns

df4 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & architecture == 'mlp' & dataset == 'ALDH1' & n_start == 64)
df4_statistics = list()

# statistics for f
for (acq in acq_funcs){
  df4_statistics[['subplot']] = c(df4_statistics[['subplot']], 'f')
  df4_statistics[['acq']] = c(df4_statistics[['acq']], acq)

  p_rand = wilcox.test(subset(df4, total_mols_screened==1000 & acquisition_method == acq)$unique_patterns,
              subset(df4, total_mols_screened==1000 & acquisition_method == 'Random')$unique_patterns,
              paired=T)$p.value
  if (is.na(p_rand)){
    p_rand = 1
  }
  if (p_rand < 0.05){
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '*')
  } else {
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '')
  }
}

## f ##
df4f = df4 %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4f = ggplot(df4f, aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))
rm(df4f)

## g ##
tot_hit_patterns_found_per_method = df4 %>%
  group_by(acquisition_method, seed) %>%
  summarise(across(.cols = hit_unique_patterns, list(max = max))) %>%
  ungroup()

perc_of_tot_patterns = c()
for (i_row in 1:nrow(df4)){
  row = df4[i_row, ]
  tot_hit_patterns = subset(tot_hit_patterns_found_per_method, acquisition_method == row$acquisition_method & seed == row$seed)$hit_unique_patterns_max
  perc = row$hit_unique_patterns/tot_hit_patterns * 100
  perc_of_tot_patterns = c(perc_of_tot_patterns, perc)
}
df4$perc_of_tot_patterns = perc_of_tot_patterns
rm(perc_of_tot_patterns)

df4g = df4 %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4g = ggplot(df4g, aes(x = total_mols_screened, y=perc_of_tot_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = perc_of_tot_patterns_mean - perc_of_tot_patterns_se, ymax = perc_of_tot_patterns_mean + perc_of_tot_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nacquired hits (% of total)', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
rm(df4g)


## h ##
pattern_enrichment_aldehyde = c()
for (x in zip(df4$total_mols_screened, df4['Aldehyde'][[1]])){
  ntot_mols = x[[1]]
  found = x[[2]]
  # break
  expected = subset(pattern_occurence, pattern == 'Aldehyde')$total_patterns*ntot_mols/100000
  enrichment = found/expected
  pattern_enrichment_aldehyde = c(pattern_enrichment_aldehyde, enrichment)
}
df4$pattern_enrichment_aldehyde = pattern_enrichment_aldehyde
rm(pattern_enrichment_aldehyde)

# statistics
# df4_stats$pattern_enrichment_aldehyde = subset(df4, total_mols_screened == 1000)$pattern_enrichment_aldehyde
for (acq in acq_funcs){
  df4_statistics[['subplot']] = c(df4_statistics[['subplot']], 'h')
  df4_statistics[['acq']] = c(df4_statistics[['acq']], acq)
  
  p_rand = wilcox.test(subset(df4, total_mols_screened==1000 & acquisition_method == acq)$pattern_enrichment_aldehyde,
                       subset(df4, total_mols_screened==1000 & acquisition_method == 'Random')$pattern_enrichment_aldehyde,
                       paired=T)$p.value
  if (is.na(p_rand)){
    p_rand = 1
  }
  if (p_rand < 0.05){
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '*')
  } else {
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '')
  }
}

df4h = df4 %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = pattern_enrichment_aldehyde, 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4h = ggplot(df4h, aes(x = total_mols_screened, y=pattern_enrichment_aldehyde_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = pattern_enrichment_aldehyde_mean - pattern_enrichment_aldehyde_se, ymax = pattern_enrichment_aldehyde_mean + pattern_enrichment_aldehyde_se),
              color=NA, alpha=0.1) + 
  labs(y = 'Aldehydes enrichment\nin acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))
rm(df4h)


## f ##
expected_nr_of_aldehydes_in_bacth = subset(pattern_occurence, pattern == 'Aldehyde')$total_patterns*64/100000
expected_aldehydes_per_hit = subset(pattern_occurence, pattern == 'Aldehyde')$hit_patterns/4986
print(paste0(expected_nr_of_aldehydes_in_bacth, "% of hits, ", expected_aldehydes_per_hit, "% of non-hits"))
rm(expected_nr_of_aldehydes_in_bacth, expected_aldehydes_per_hit)

aldehyde_hit_enrichment = c()
for (x in zip(df4$hits_discovered, df4['hit_Aldehyde'][[1]])){
  ntot_hits = x[[1]]
  found = x[[2]]

  expected = ntot_hits*subset(pattern_occurence, pattern == 'Aldehyde')$hit_patterns/4986
  enrichment = found/expected
  aldehyde_hit_enrichment = c(aldehyde_hit_enrichment, enrichment)
}
df4$aldehyde_hit_enrichment = aldehyde_hit_enrichment
rm(aldehyde_hit_enrichment)

# statistics
for (acq in acq_funcs){
  df4_statistics[['subplot']] = c(df4_statistics[['subplot']], 'i')
  df4_statistics[['acq']] = c(df4_statistics[['acq']], acq)
  
  p_rand = wilcox.test(subset(df4, total_mols_screened == 1000 & acquisition_method == acq)$aldehyde_hit_enrichment,
                       subset(df4, total_mols_screened == 1000 & acquisition_method == 'Random')$aldehyde_hit_enrichment,
                       paired=T)$p.value
  if (is.na(p_rand)){
    p_rand = 1
  }
  if (p_rand < 0.05){
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '*')
  } else {
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '')
  }
}

df4i = df4 %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = aldehyde_hit_enrichment, 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4i = ggplot(df4i, aes(x = total_mols_screened, y=aldehyde_hit_enrichment_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = aldehyde_hit_enrichment_mean - aldehyde_hit_enrichment_se, ymax = aldehyde_hit_enrichment_mean + aldehyde_hit_enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Aldehyde enrichment\nin acquired hits', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), expand=F) +
  # scale_y_continuous(breaks = seq(0,20, by=10)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))
rm(df4i)


## j ##
pattern_enrichment_sulfonamide = c()
for (x in zip(df4$total_mols_screened, df4['Sulfonamide'][[1]])){
  ntot_mols = x[[1]]
  found = x[[2]]
  expected = subset(pattern_occurence, pattern == 'Sulfonamide')$total_patterns*ntot_mols/100000
  enrichment = found/expected
  pattern_enrichment_sulfonamide = c(pattern_enrichment_sulfonamide, enrichment)
}
df4$pattern_enrichment_sulfonamide = pattern_enrichment_sulfonamide
rm(pattern_enrichment_sulfonamide)


# statistics
for (acq in acq_funcs){
  df4_statistics[['subplot']] = c(df4_statistics[['subplot']], 'j')
  df4_statistics[['acq']] = c(df4_statistics[['acq']], acq)
  
  p_rand = wilcox.test(subset(df4, total_mols_screened == 1000 & acquisition_method == acq)$pattern_enrichment_sulfonamide,
                       subset(df4, total_mols_screened == 1000 & acquisition_method == 'Random')$pattern_enrichment_sulfonamide,
                       paired=T)$p.value
  if (is.na(p_rand)){
    p_rand = 1
  }
  if (p_rand < 0.05){
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '*')
  } else {
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '')
  }
}

df4j = df4 %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = pattern_enrichment_sulfonamide, 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4j = ggplot(df4j, aes(x = total_mols_screened, y=pattern_enrichment_sulfonamide_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = pattern_enrichment_sulfonamide_mean - pattern_enrichment_sulfonamide_se, ymax = pattern_enrichment_sulfonamide_mean + pattern_enrichment_sulfonamide_se),
              color=NA, alpha=0.1) +
  labs(y = 'Sulfonamide enrichment\nin acquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(ylim = c(0.5, 4), xlim = c(64, 1010), expand=F) +
  scale_y_continuous(breaks = seq(0,3, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
rm(df4j)


## k ##
expected_nr_of_aldehydes_in_bacth = subset(pattern_occurence, pattern == 'Sulfonamide')$total_patterns*64/100000
expected_aldehydes_per_hit = subset(pattern_occurence, pattern == 'Sulfonamide')$hit_patterns/4986
print(paste0(expected_nr_of_aldehydes_in_bacth, "% of hits, ", expected_aldehydes_per_hit, "% of non-hits"))
rm(expected_nr_of_aldehydes_in_bacth, expected_aldehydes_per_hit)

sulfonamide_hit_enrichment = c()
for (x in zip(df4$hits_discovered, df4['hit_Sulfonamide'][[1]])){
  ntot_hits = x[[1]]
  found = x[[2]]
  
  expected = ntot_hits*subset(pattern_occurence, pattern == 'Sulfonamide')$hit_patterns/4986
  enrichment = found/expected
  sulfonamide_hit_enrichment = c(sulfonamide_hit_enrichment, enrichment)
}
df4$sulfonamide_hit_enrichment = sulfonamide_hit_enrichment
rm(sulfonamide_hit_enrichment)


# statistics
for (acq in acq_funcs){
  df4_statistics[['subplot']] = c(df4_statistics[['subplot']], 'k')
  df4_statistics[['acq']] = c(df4_statistics[['acq']], acq)
  
  p_rand = wilcox.test(subset(df4, total_mols_screened == 1000 & acquisition_method == acq)$sulfonamide_hit_enrichment,
                       subset(df4, total_mols_screened == 1000 & acquisition_method == 'Random')$sulfonamide_hit_enrichment,
                       paired=T)$p.value
  if (is.na(p_rand)){
    p_rand = 1
  }
  if (p_rand < 0.05){
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '*')
  } else {
    df4_statistics[['wilx_rand']] = c(df4_statistics[['wilx_rand']], '')
  }
}
df4_statistics = data.frame(df4_statistics)

df4k = df4 %>%
  group_by(acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = sulfonamide_hit_enrichment, 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig4k = ggplot(df4k, aes(x = total_mols_screened, y=sulfonamide_hit_enrichment_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = sulfonamide_hit_enrichment_mean - sulfonamide_hit_enrichment_se, ymax = sulfonamide_hit_enrichment_mean + sulfonamide_hit_enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Sulfonamide enrichment\nin acquired hits', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(ylim = c(0.5, 4), xlim = c(64, 1010), expand=F) +
  scale_y_continuous(breaks = seq(0,3, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
rm(df4k)


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
rm(fig4d, fig4e, fig4f, fig4g, fig4h, fig4i, fig4j, fig4k, fig4de, fig4fk)

dev.print(pdf, 'al_v4_fig4_.pdf', width = 180/25.4, height = 123/25.4)


#### Fig 5: Start set size ####

df5 = subset(df, train_cycle <= 15 & batch_size == 64 & bias == 'Innate')
df5 = subset(df5, acquisition_method %in% c("Similarity", "BALD (least mutual information)"))
df5 = subset(df5, dataset %in% c("PKM2" ,"ALDH1", "VDR"))

df5 = df5[c(2, 32, 34, 35, 36, 38, 41, 151)]
df5$id = paste(df5$architecture, df5$n_start, df5$seed, df5$dataset)

df5_sim = subset(df5, acquisition_method == "Similarity")
df5_bald = subset(df5, acquisition_method == "BALD (least mutual information)")


df5_bald = df5_bald[df5_bald$id %in% df5_sim$id, ]
df5_sim = df5_sim[df5_sim$id %in% df5_bald$id, ]

df5_sim = df5_sim[order(df5_sim$id),]
df5_bald = df5_bald[order(df5_bald$id),]

# difference in enrichment
df5_sim$enrichment = df5_bald$enrichment - df5_sim$enrichment 


###### Statistics ######

# Statistical tests between 64 and the other start 
df5_statistics = list()
for (dataset_ in c("PKM2" ,"ALDH1", "VDR")){
  for (algo in c('mlp', 'gcn')){

    for (start_size in c(2, 4, 8, 16, 32, 64)){
      df5_statistics[['datasets']] = c(df5_statistics[['datasets']] , dataset_)
      df5_statistics[['algos']] = c(df5_statistics[['algos']] , algo)
      df5_statistics[['start_size']] = c(df5_statistics[['start_size']] , start_size)
        p = wilcox.test(subset(df5_sim, train_cycle == 15 & dataset == dataset_ & architecture == algo & n_start == start_size)$enrichment,
                        subset(df5_sim, train_cycle == 15 & dataset == dataset_ & architecture == algo & n_start == 64)$enrichment,
                        paired=T)$p.value
        if (is.na(p)){
          p=1
        }
        if (p < 0.05){
          df5_statistics[['wilx_64']] = c(df5_statistics[['wilx_64']], '*')
        } else {
          df5_statistics[['wilx_64']] = c(df5_statistics[['wilx_64']], '')
        }
    }
  }
}
df5_statistics = data.frame(df5_statistics)


####### plots ########

df5 = df5_sim %>%
  group_by(acquisition_method, n_start, dataset, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

col_gradient = c("#bbbbbb", "#94d2bd", "#0a9396", "#005F73", "#f8cd48", "#ee9b00")


fig5a = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = 'Increase in enrichment\nover similarity search', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

fig5b = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
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

fig5b = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
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

fig5c = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = 'Increase in enrichment\nover similarity search', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

fig5d = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
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

fig5e = ggplot(subset(df5, dataset == 'VDR' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.6) +
  labs(y = 'Increase in enrichment\nover similarity search', x='active learning cycle', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

fig5f = ggplot(subset(df5, dataset == 'VDR' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
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

fig5 = plot_grid(fig5a, fig5b, 
                 fig5c, fig5d,
                 fig5e, fig5f,
                 labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
                 ncol=2, label_size=8)

fig5

dev.print(pdf, 'al_v4_fig5.pdf', width = 88/25.4, height = 123/25.4)
rm(fig5a, fig5b, fig5c, fig5d, fig5e, fig5f, df5_sim, df5_bald)



#### Supplementary ####


##### Random ####



# compute scaffold stuff in Python
df_rand <- read_csv("processed_results_nosmiles_rebuttal2.csv")
df_rand = subset(df_rand, n_start <= 64)

df_rand$seed = factor(df_rand$seed, levels=unique(df_rand$seed))
df_rand$bias = stringr::str_to_title(gsub('random', 'innate', df_rand$bias))
df_rand$bias = factor(df_rand$bias, levels=c('Innate', 'Small', 'Large'))
df_rand$batch_size = factor(df_rand$batch_size, levels=c(16, 32, 64))
df_rand$n_start = factor(df_rand$n_start, levels=c(2, 4, 8, 16, 32, 64))
df_rand[df_rand$retrain == "FALSE", ]$acquisition_method = 'Exploit, no retraining'
df_rand$acquisition_method = gsub('bald', 'BALD (least mutual information)', df_rand$acquisition_method)
df_rand$acquisition_method = gsub('exploitation', 'Exploit (best predictions)', df_rand$acquisition_method)
df_rand$acquisition_method = gsub('exploration', 'Explore (most uncertain)', df_rand$acquisition_method)
df_rand$acquisition_method = gsub('random', 'Random', df_rand$acquisition_method)
df_rand$acquisition_method = gsub('similarity', 'Similarity', df_rand$acquisition_method)


df_rand = subset(df_rand, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64)

df_rand = df_rand %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()
# table(df_rand$dataset)


sfigranda = ggplot(subset(df_rand, dataset == 'PKM2'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
# u r b l
sfigrandb = ggplot(subset(df_rand, dataset == 'ALDH1'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfigrandc = ggplot(subset(df_rand, dataset == 'VDR'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfigrandd = ggplot(subset(df_rand, dataset == 'FEN1'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfigrande = ggplot(subset(df_rand, dataset == 'GBA'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfigrandf = ggplot(subset(df_rand, dataset == 'KAT2A'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfigrandg = ggplot(subset(df_rand, dataset == 'IDH1'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfigrandh = ggplot(subset(df_rand, dataset == 'OPRK1'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

sfigrandi = ggplot(subset(df_rand, dataset == 'ADRB2'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

# "#005f73" "#6f9cbc" "#6d7889" "#ee9b00" "#bbbbbb" "#f8cd48"


sfigrand = plot_grid(sfigranda, sfigrandb, sfigrandc, 
                     sfigrandd, sfigrande, sfigrandf,
                     sfigrandg, sfigrandh, sfigrandi,
                 labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' ,'i'), 
                 ncol=3, label_size=8)
sfigrand

# 180 mm/ 88 mm
dev.print(pdf, 'al_v4_sfigrand.pdf', width = 110/25.4, height = 100/25.4)
rm(sfigranda, sfigrandab, sfigrandac, sfigrandad, sfigrandae, sfigrandaf,
   sfigrandg, sfigrandh, sfigrandi)




##### S1. Batch size #####

dfs1a = subset(df, train_cycle == 15 & architecture == 'mlp' & dataset == 'ALDH1' & bias == 'Innate' & n_start == 64)
dfs1a$batch_size = as.numeric(as.character(dfs1a$batch_size))

sfig1a = ggplot(dfs1a, aes(x = batch_size, y = enrichment, color = acquisition_method, fill = acquisition_method))+
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1) +
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
  geom_smooth(method=lm, level=0.95, size=0.6, alpha=0.1, linetype='dashed') +
  labs(y = '', x='Molecules acquired per cycle') +
  coord_cartesian(xlim = c(12, 68), ylim = c(0.45, 0.65), expand=F) +
  scale_y_continuous(breaks = seq(0.45,0.65, by=0.05)) +
  scale_x_continuous(breaks = c(16, 32, 64)) +
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))


sfig1 = plot_grid(sfig1a, sfig1b, sfig1c, sfig1d,
                  ncol=2, label_size = 8, labels= c('a', 'b', 'c', 'd'))
sfig1

dev.print(pdf, 'al_v4_sfig1.pdf', width = 88/25.4, height = 80/25.4)
rm(sfig1a, dfs1a, sfig1b, dfs1b, sfig1c, dfs1c, sfig1d, dfs1d)


##### S2. UMAP ######

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


sfig2 = plot_grid(p_umap_pkm2, p_umap_aldh1, p_umap_vdr, rel_widths = c(2, 2, 2),
                 labels = c('a', 'b', 'c'),
                 ncol=3, label_size=8)
sfig2

dev.print(pdf, 'al_v4_sfig2.pdf', width = 180/25.4, height = 60/25.4)
rm(umap_vdr, p_umap_pkm2, umap_aldh1, p_umap_aldh1, umap_vdr, p_umap_vdr)


##### S3. Hit similarity distribution #####

dfs3 <- read_csv("rebuttal_data.csv")

sfig3a = ggplot(subset(dfs3, dataset == 'ALDH1'), aes(x=`Tanimoto similarity`, fill=kind)) +
  geom_density(alpha=0.4) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb")) +
  labs(x='Mean Tanimoto similarity', fill='', title='ALDH1') +
  scale_x_continuous(limits = c(0,0.3), expand = expansion(mult = c(0.01, 0.01)))+
  default_plot_theme + theme(legend.position = 'none')

sfig3b = ggplot(subset(dfs3, dataset == 'PKM2'), aes(x=`Tanimoto similarity`, fill=kind)) +
  geom_density(alpha=0.4) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb")) +
  labs(x='Mean Tanimoto similarity', fill='', title='PKM2') +
  scale_x_continuous(limits = c(0,0.3), expand = expansion(mult = c(0.01, 0.01)))+
  default_plot_theme+ theme(legend.position = 'none')

sfig3c = ggplot(subset(dfs3, dataset == 'VDR'), aes(x=`Tanimoto similarity`, fill=kind)) +
  geom_density(alpha=0.4) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb")) +
  labs(x='Mean Tanimoto similarity', fill='', title='VDR') +
  scale_x_continuous(limits = c(0,0.3), expand = expansion(mult = c(0.01, 0.01)))+
  default_plot_theme_legend

sfig3 = plot_grid(sfig3a, sfig3b, sfig3c, ncol = 3, rel_widths = c(1, 1, 1.6))
sfig3

dev.print(pdf, 'al_v4_sfig3.pdf', width = 7, height = 2) # width = 7.205, height = 4
rm(sfig3a, sfig3b, sfig3c, dfs3)


##### S4. Structural Diversity ######
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


figs4 = plot_grid(figs4a, figs4b, figs4c,  figs4d, 
                  figs4e, figs4f, figs4g,  figs4h, 
                  figs4i, figs4j, figs4k,  figs4l, 
                  ncol=4, label_size = 8, 
                  labels= c('a', 'b', 'c', 'd',
                            'e', 'f', 'g', 'h', 
                            'i', 'j', 'k', 'l'))

dev.print(pdf, 'al_v4_sfig4.pdf', width = 180/25.4, height = 130/25.4)
rm(figs4a, figs4b, figs4c,  figs4d, figs4e, figs4f, figs4g,  figs4h, 
   figs4i, figs4j, figs4k,  figs4l, dfs4ac, dfs4bd, dfs4eg, dfs4fh, 
   dfs4ik, dfs4jl)


# ##### S5. Random Forest #####
# 
# acq_rf_last = subset(df, train_cycle >= 0 & architecture == 'rf' & total_mols_screened == 1000 & batch_size == 64 & bias == 'Innate' & n_start == 64)
# acq_rf_statistics = list()
# for (dataset_ in unique(acq_rf_last$dataset)){
#   for (acq1 in acq_funcs){
#     acq_rf_statistics[['datasets']] = c(acq_rf_statistics[['datasets']] , dataset_)
#     acq_rf_statistics[['algos']] = c(acq_rf_statistics[['algos']] , 'rf')
#     acq_rf_statistics[['acqs']] = c(acq_rf_statistics[['acqs']] , acq1)
#     
#     a_ = subset(acq_rf_last, dataset == dataset_ & acquisition_method == acq1)$test_roc_auc
#     b_ = subset(acq_rf_last, dataset == dataset_ & acquisition_method == "Random")$test_roc_auc
#     
#     if (length(a_) + length(b_) > 0 & length(a_) == length(b_)){
#       
#       p = wilcox.test(a_, b_, paired=T)$p.value
#       
#       if (is.na(p)){
#         p=1
#       }
#       if (p < 0.05){
#         acq_rf_statistics[["wilcox"]] = c(acq_rf_statistics[["wilcox"]], '*')
#       } else {
#         acq_rf_statistics[["wilcox"]] = c(acq_rf_statistics[["wilcox"]], '')
#       }
#     } else {
#       acq_rf_statistics[["wilcox"]] = c(acq_rf_statistics[["wilcox"]], '')
#     }
#   }
# }
# acq_rf_statistics = data.frame(acq_rf_statistics)
# rm(acq_rf_last)
# 
# 
# dfs5 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64) %>%
#   group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
#   summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
#                      'test_balanced_accuracy', 'enrichment'), 
#                    list(mean = mean, sd = sd, se = se))) %>% ungroup()
# 
# sfig5a = ggplot(subset(dfs5, dataset == 'PKM2' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
#   geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(size=0.6) + 
#   coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
#   scale_y_continuous(breaks = seq(0,7, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
# 
# sfig5b = ggplot(subset(dfs5, dataset == 'PKM2' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
#   geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Mean ROC AUC\non the test set', x='', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(size=0.6) +
#   coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
#   scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
#   scale_x_continuous(breaks = seq(1,15, by=2)) +
#   default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))
# 
# sfig5c = ggplot(subset(dfs5, dataset == 'ALDH1' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
#   geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(size=0.6) + 
#   coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
#   scale_y_continuous(breaks = seq(0,7, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))
# 
# sfig5d = ggplot(subset(dfs5, dataset == 'ALDH1' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
#   geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Mean ROC AUC\non the test set', x='', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(size=0.6) +
#   coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
#   scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
#   scale_x_continuous(breaks = seq(1,15, by=2)) +
#   default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))
# 
# sfig5e = ggplot(subset(dfs5, dataset == 'VDR' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
#   geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Mean enrichment in\nacquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(size=0.6) + 
#   coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
#   scale_y_continuous(breaks = seq(0,7, by=1)) +
#   default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))
# 
# sfig5f = ggplot(subset(dfs5, dataset == 'VDR' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
#   geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Mean ROC AUC\non the test set', x='Active learning cycle', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(size=0.6) +
#   coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
#   scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
#   scale_x_continuous(breaks = seq(1,15, by=2)) +
#   default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))
# 
# 
# sfig5 = plot_grid(sfig5a, sfig5b, 
#                  sfig5c, sfig5d,
#                  sfig5e, sfig5f,
#                  labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
#                  rel_widths = c(1, 1.125),
#                  ncol=2, label_size=8)
# sfig5
# 
# dev.print(pdf, 'al_v4_sfig5.pdf', width = 88/25.4, height = 110/25.4)
# rm(sfig5a, sfig5b, sfig5c, sfig5d, sfig5e, sfig5f, dfs5)



##### S6: Other datasets #####

df6 = df %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

# data.frame(datasets = c('FEN1', 'GBA', 'KAT2A', 'IDH1', 'OPRK1', 'ADRB2'),
#            hits = c(100, 55 ,54, 11, 9 ,5))

# first row
sfig6a = ggplot(subset(dfs6, dataset == 'FEN1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

sfig6b = ggplot(subset(dfs6, dataset == 'FEN1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

sfig6c = ggplot(subset(dfs6, dataset == 'GBA' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))

sfig6d = ggplot(subset(dfs6, dataset == 'GBA' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

# second row
sfig6e = ggplot(subset(dfs6, dataset == 'KAT2A' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfig6f = ggplot(subset(dfs6, dataset == 'KAT2A' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfig6g = ggplot(subset(dfs6, dataset == 'IDH1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig6h = ggplot(subset(dfs6, dataset == 'IDH1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

# third row
sfig6i = ggplot(subset(dfs6, dataset == 'OPRK1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfig6j = ggplot(subset(dfs6, dataset == 'OPRK1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfig6k = ggplot(subset(dfs6, dataset == 'ADRB2' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 6), expand=F) + 
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))

sfig6l = ggplot(subset(dfs6, dataset == 'ADRB2' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

sfig6 = plot_grid(sfig6a, sfig6b, sfig6c, sfig6d,
                  sfig6e, sfig6f, sfig6g, sfig6h,
                  sfig6i, sfig6j, sfig6k, sfig6l,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'), 
                  ncol=4, label_size=8)
sfig6

# 180 mm/ 88 mm
dev.print(pdf, 'al_v4_sfig6.pdf', width = 180/25.4, height = 123/25.4)
rm(fig6a, fig6b, fig6c, fig6d, fig6e, fig6f)


##### S7. ROC AUC #####

###### Statistics Acq functions #######

acq_roc_last = subset(df, train_cycle >= 0 & total_mols_screened == 1000 & batch_size == 64 & bias == 'Innate' & n_start == 64)
acq_roc_statistics = list()
for (dataset_ in c('PKM2', 'ALDH1', 'VDR')){
  for (algo in c('mlp', 'gcn', 'rf')){
    for (acq1 in acq_funcs){
      acq_roc_statistics[['datasets']] = c(acq_roc_statistics[['datasets']] , dataset_)
      acq_roc_statistics[['algos']] = c(acq_roc_statistics[['algos']] , algo)
      acq_roc_statistics[['acqs']] = c(acq_roc_statistics[['acqs']] , acq1)
      
      a_ = subset(acq_roc_last, dataset == dataset_ & architecture == algo & acquisition_method == acq1)$test_roc_auc
      b_ = subset(acq_roc_last, dataset == dataset_ & architecture == algo & acquisition_method == "Random")$test_roc_auc
      
      if (length(a_) + length(b_) > 0 & length(a_) == length(b_)){
        
        p = wilcox.test(a_, b_, paired=T)$p.value
        
        if (is.na(p)){
          p=1
        }
        if (p < 0.05){
          acq_roc_statistics[["wilcox"]] = c(acq_roc_statistics[["wilcox"]], '*')
        } else {
          acq_roc_statistics[["wilcox"]] = c(acq_roc_statistics[["wilcox"]], '')
        }
      } else {
        acq_roc_statistics[["wilcox"]] = c(acq_roc_statistics[["wilcox"]], '')
      }
    }
  }
}
acq_roc_statistics = data.frame(acq_roc_statistics)
rm(acq_roc_last)


dfs7 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64) %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

sfig7a = ggplot(subset(dfs7, dataset == 'PKM2' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

sfig7b = ggplot(subset(dfs7, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

sfig7c = ggplot(subset(dfs7, dataset == 'PKM2' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfig7d = ggplot(subset(dfs7, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

sfig7e = ggplot(subset(dfs7, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

sfig7f = ggplot(subset(dfs7, dataset == 'ALDH1' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig7g = ggplot(subset(dfs7, dataset == 'VDR' & architecture == 'mlp'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

sfig7h = ggplot(subset(dfs7, dataset == 'VDR' & architecture == 'gcn'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
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

sfig7i = ggplot(subset(dfs7, dataset == 'VDR' & architecture == 'rf'), aes(x = train_cycle, y=test_roc_auc_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = test_roc_auc_mean - test_roc_auc_se, ymax = test_roc_auc_mean + test_roc_auc_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.6) + 
  coord_cartesian(xlim = c(1, 15.5), ylim = c(0.4, 0.7), expand=F) +
  scale_y_continuous(breaks = seq(0.4, 0.7, by=0.05)) +
  scale_x_continuous(breaks = seq(1,15, by=2)) +
  scale_linetype_manual(values=c("dotted"))+
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))
sfig7 = plot_grid(sfig7a, sfig7b, sfig7c,
                  sfig7d, sfig7e, sfig7f,
                  sfig7g, sfig7h, sfig7i,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'), 
                  ncol=3, label_size=8)
sfig7

# 180 mm/ 88 mm
dev.print(pdf, 'al_v4_sfig7.pdf', width = 110/25.4, height = 90/25.4)


##### S8: substructures in other datasets #####

dfs8 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64)
dfs8_statistics = list()

# statistics for s8
for (dataset_ in c('PKM2', "ALDH1", "VDR")){
  for (algo in c('mlp', 'gcn', 'rf')){
    for (acq in acq_funcs){
      dfs8_statistics[['dataset']] = c(dfs8_statistics[['dataset']], dataset_)
      dfs8_statistics[['algo']] = c(dfs8_statistics[['algo']], algo)
      dfs8_statistics[['acq']] = c(dfs8_statistics[['acq']], acq)
      
      p_rand = wilcox.test(subset(dfs8, architecture == algo & dataset == dataset_ & total_mols_screened==1000 & acquisition_method == acq)$unique_patterns,
                           subset(dfs8, architecture == algo & dataset == dataset_ & total_mols_screened==1000 & acquisition_method == 'Random')$unique_patterns,
                           paired=T)$p.value
      if (is.na(p_rand)){
        p_rand = 1
      }
      if (p_rand < 0.05){
        dfs8_statistics[['wilx_rand']] = c(dfs8_statistics[['wilx_rand']], '*')
      } else {
        dfs8_statistics[['wilx_rand']] = c(dfs8_statistics[['wilx_rand']], '')
      }
    }
  }
}
dfs8_statistics = data.frame(dfs8_statistics)

df8 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64) %>%
  group_by(architecture, dataset, acquisition_method, total_mols_screened) %>%
  summarise(across(.cols = where(is.numeric), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

sfig8a = ggplot(subset(df8, architecture == 'mlp' & dataset == 'PKM2'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig8b = ggplot(subset(df8, architecture == 'gcn' & dataset == 'PKM2'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig8c = ggplot(subset(df8, architecture == 'rf' & dataset == 'PKM2'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig8d = ggplot(subset(df8, architecture == 'mlp' & dataset == 'ALDH1'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig8e = ggplot(subset(df8, architecture == 'gcn' & dataset == 'ALDH1'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig8f = ggplot(subset(df8, architecture == 'rf' & dataset == 'ALDH1'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='\n', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig8g = ggplot(subset(df8, architecture == 'mlp' & dataset == 'VDR'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = 'Unique patterns in\nall acquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig8h = ggplot(subset(df8, architecture == 'gcn' & dataset == 'VDR'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

sfig8i = ggplot(subset(df8, architecture == 'rf' & dataset == 'VDR'), aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
  geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
              color=NA, alpha=0.1) +
  labs(y = '\n', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(linewidth=0.6) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
  scale_y_continuous(breaks = seq(38,48, by=2)) +
  default_plot_theme + theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

# 
sfig8 = plot_grid(sfig8a, sfig8b, sfig8c,
                  sfig8d, sfig8e, sfig8f,
                  sfig8g, sfig8h, sfig8i,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'),
                  ncol=3, label_size=8)
sfig8

dev.print(pdf, 'al_v4_sfig8.pdf', width = 110/25.4, height = 100/25.4)


##### S9: start size #####

dfs9 = subset(df, total_mols_screened == 1000 & batch_size == 64 & bias == 'Innate')
dfs9 = subset(dfs9, acquisition_method %in% c("Random", "Similarity", "BALD (least mutual information)"))

dfs9 = dfs9 %>%
  group_by(acquisition_method, n_start, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()


sfig9a = ggplot(subset(dfs9, dataset == 'PKM2' & architecture == 'mlp'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"))
# u r b l
sfig9b = ggplot(subset(dfs9, dataset == 'PKM2' & architecture == 'gcn'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"))

sfig9c = ggplot(subset(dfs9, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"))

sfig9d = ggplot(subset(dfs9, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"))

sfig9e = ggplot(subset(dfs9, dataset == 'VDR' & architecture == 'mlp'),  aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = 'Mean enrichment in\n1000 acquired molecules', x='molecules in the start set', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"))

sfig9f = ggplot(subset(dfs9, dataset == 'VDR' & architecture == 'gcn'),  aes(x = n_start, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture, group=acquisition_method))+
  geom_pointrange(aes(ymin=enrichment_mean-enrichment_se, ymax=enrichment_mean+enrichment_se), size = 0, linewidth=0.75, alpha=0.5) + 
  geom_point(size=1) +
  labs(y = '', x='molecules in the start set', color = 'Method', fill = 'Method')+
  scale_color_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  scale_fill_manual(values=c("#005f73", "#bbbbbb", "#f8cd48")) +
  coord_cartesian(xlim = c(0.5, 6), ylim = c(0, 6), expand=F) +
  scale_y_continuous(breaks = seq(0,6, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"))



sfig9 = plot_grid(sfig9a, sfig9b, 
                  sfig9c, sfig9d,
                  sfig9e, sfig9f,
                  labels = c('a', 'b', 'c', 'd', 'e', 'f'), 
                  ncol=2, label_size=8)
sfig9

# 180 mm/ 88 mm
dev.print(pdf, 'al_v4_sfig9.pdf', width = 88/25.4, height = 123/25.4)



##### S10: low-yield datasets #######

dfs10 = subset(df, batch_size == 64 & bias == 'Innate' & total_mols_screened==1000 & n_start==64) %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

# df of yields
dataset_yield =  data.frame(dataset = c('ALDH1', 'VDR', 'PKM2', 'FEN1', 'GBA', 'KAT2A', 'IDH1', 'OPRK1', 'ADRB2'),
                            yield = c(4986, 239, 223, 100, 55 ,54, 11, 9 ,5))
dataset_yield$label = paste0(dataset_yield$dataset, '\n(', dataset_yield$yield, ' hits)')
dfs10$yield = dataset_yield$yield[match(df10$dataset, dataset_yield$dataset)]
dfs10$yield_label = dataset_yield$label[match(df10$dataset, dataset_yield$dataset)]
x_axis_values = round(sort(unique(log(dfs10$yield))), 2)

sfig10a = ggplot(subset(dfs10, architecture == 'mlp'), aes(x = log(yield), y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture)) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='log(dataset yield)', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_smooth(method=lm, se=T, fullrange=TRUE, alpha=0.2, level=0.95, size=0.6)+
  coord_cartesian(xlim = c(x_axis_values[1], 8.6), ylim = c(0, 8), expand=F) +
  scale_y_continuous(breaks = seq(0,8, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0, 0.25, 0.25), "cm"))

sfig10b = ggplot(subset(dfs10, architecture == 'gcn'), aes(x = log(yield), y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture)) +
  labs(y = '', x='log(dataset yield)', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_smooth(method=lm, se=T, fullrange=TRUE, alpha=0.2, level=0.95, size=0.6)+
  coord_cartesian(xlim = c(x_axis_values[1], 8.6), ylim = c(0, 8), expand=F) +
  scale_y_continuous(breaks = seq(0,8, by=1)) +
  default_plot_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0), "cm"))

sfig10 = plot_grid(sfig10a, sfig10b, 
                 labels = c('a', 'b'), 
                 ncol=2, label_size=8)

sfig10

dev.print(pdf, 'al_v4_sfig10.pdf', width = 88/25.4, height = 44/25.4)


# # statistics
# for (algo in c('mlp', 'gcn')){
#   for (acq in acq_funcs){
# 
#     p <- cor.test(subset(dfs10, architecture == algo & acquisition_method == acq)$enrichment_mean,
#                   subset(dfs10, architecture == algo & acquisition_method == acq)$yield,
#                   method="kendall")$p.value
#     if (p < 0.05){print(paste(algo, acq, p))}
#   }
# }



#### Junkjard ####


# dataset_ = 'VDR'
# 
# acq_funcs
# 
# wilcox.test(subset(dfs8, architecture == 'rf' & dataset == dataset_ & total_mols_screened==1000 & acquisition_method == 'Explore (most uncertain)')$unique_patterns,
#             subset(dfs8, architecture == 'gcn' & dataset == dataset_ & total_mols_screened==1000 & acquisition_method == 'Explore (most uncertain)')$unique_patterns,
#             paired=T)$p.value
# 

# dfx = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & architecture == 'mlp' & dataset == 'ALDH1' & n_start == 64)

# Old s fig 8 (the not normalized)
# dfxx = dfx %>%
#   group_by(acquisition_method, total_mols_screened) %>%
#   summarise(across(.cols = where(is.numeric), 
#                    list(mean = mean, sd = sd, se = se))) %>% ungroup()
# 
# figsx = ggplot(dfxx, aes(x = total_mols_screened, y=hit_unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
#   geom_ribbon(aes(ymin = hit_unique_patterns_mean - hit_unique_patterns_se, ymax = hit_unique_patterns_mean + hit_unique_patterns_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Unique substructures in\nacquired hits', x='n molecules screened', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(linewidth=0.6) + 
#   coord_cartesian(xlim = c(64, 1010), expand=F) +
#   # scale_y_continuous(breaks = seq(15,45, by=5)) +
#   default_plot_theme + theme(plot.margin = margins_d)
# 
# dev.print(pdf, 'al_Sup_fig5.pdf', width = 44/25.4, height = 44/25.4)


# # Line plots #
# 
# pattern_occurence <- read_csv("pattern_occurence_ALDH1.csv", 
#                               col_types = cols(...1 = col_skip(), index = col_skip()))
# pattern_occurence$total_patterns = pattern_occurence$hit_patterns + pattern_occurence$nonhit_patterns
# 
# dfx = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & architecture == 'rf' & dataset == 'ALDH1' & n_start == 64)
# 
# margins_c = unit(c(0.25, 0.25, 0, 0.25), "cm")
# margins_d = unit(c(0.25, 0.25, 0, 0), "cm")
# margins_e = unit(c(0.25, 0.25, 0, 0.25), "cm")
# margins_f = unit(c(0.25, 0.25, 0, 0.25), "cm")
# margins_g = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
# margins_h = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
# 
# 
# 
# # Statistics
# dfxc_stats = subset(dfx, total_mols_screened==1000)
# dfx_stats = list()
# for (acq in acq_funcs){
#   dfx_stats[['subplot']] = c(dfx_stats[['subplot']], 'f')
#   dfx_stats[['acq']] = c(dfx_stats[['acq']], acq)
#   
#   p_rand = wilcox.test(subset(dfxc_stats, acquisition_method == acq)$unique_patterns,
#                        subset(dfxc_stats, acquisition_method == 'Random')$unique_patterns,
#                        paired=T)$p.value
#   if (is.na(p_rand)){
#     p_rand = 1
#   }
#   if (p_rand < 0.05){
#     dfx_stats[['wilx_rand']] = c(dfx_stats[['wilx_rand']], '*')
#   } else {
#     dfx_stats[['wilx_rand']] = c(dfx_stats[['wilx_rand']], '')
#   }
#   
#   p_sim = wilcox.test(subset(dfxc_stats, acquisition_method == acq)$unique_patterns,
#                       subset(dfxc_stats, acquisition_method == 'Similarity')$unique_patterns,
#                       paired=T)$p.value
#   if (is.na(p_sim)){
#     p_sim = 1
#   }
#   if (p_sim < 0.05){
#     dfx_stats[['wilx_sim']] = c(dfx_stats[['wilx_sim']], '*')
#   } else {
#     dfx_stats[['wilx_sim']] = c(dfx_stats[['wilx_sim']], '')
#   }
# }
# 
# 
# ## b ##
# dfxc = dfx %>%
#   group_by(acquisition_method, total_mols_screened) %>%
#   summarise(across(.cols = where(is.numeric), 
#                    list(mean = mean, sd = sd, se = se))) %>% ungroup()
# 
# fig3c = ggplot(dfxc, aes(x = total_mols_screened, y=unique_patterns_mean, color=acquisition_method, fill=acquisition_method))+
#   geom_ribbon(aes(ymin = unique_patterns_mean - unique_patterns_se, ymax = unique_patterns_mean + unique_patterns_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Unique patterns in\nall acquired molecules', x='\n', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(linewidth=0.6) + 
#   coord_cartesian(xlim = c(64, 1010), ylim = c(38, 48), expand=F) + 
#   scale_y_continuous(breaks = seq(38,48, by=2)) +
#   default_plot_theme + theme(plot.margin = margins_c)
# 
# ## d ##
# tot_hit_patterns_found_per_method = dfx %>%
#   group_by(acquisition_method, seed) %>%
#   summarise(across(.cols = hit_unique_patterns, list(max = max))) %>%
#   ungroup()
# 
# perc_of_tot_patterns = c()
# for (i_row in 1:nrow(dfx)){
#   row = dfx[i_row, ]
#   tot_hit_patterns = subset(tot_hit_patterns_found_per_method, acquisition_method == row$acquisition_method & seed == row$seed)$hit_unique_patterns_max
#   perc = row$hit_unique_patterns/tot_hit_patterns * 100
#   perc_of_tot_patterns = c(perc_of_tot_patterns, perc)
# }
# dfx$perc_of_tot_patterns = perc_of_tot_patterns
# 
# dfxd = dfx %>%
#   group_by(acquisition_method, total_mols_screened) %>%
#   summarise(across(.cols = where(is.numeric), 
#                    list(mean = mean, sd = sd, se = se))) %>% ungroup()
# 
# fig3d = ggplot(dfxd, aes(x = total_mols_screened, y=perc_of_tot_patterns_mean, color=acquisition_method, fill=acquisition_method))+
#   geom_ribbon(aes(ymin = perc_of_tot_patterns_mean - perc_of_tot_patterns_se, ymax = perc_of_tot_patterns_mean + perc_of_tot_patterns_se),
#               color=NA, alpha=0.1) +
#   labs(y = 'Unique patterns in\nacquired hits (% of total)', x='\n', color = 'Method', fill = 'Method')+
#   scale_color_manual(values=custom_colours) +
#   scale_fill_manual(values=custom_colours) +
#   geom_line(linewidth=0.6) + 
#   coord_cartesian(xlim = c(64, 1010), expand=F) +
#   default_plot_theme + theme(plot.margin = margins_d)
# 

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

## Supp. Fig. 1: UMAP ##

# df_UMAP_PKM2 <- read_csv("UMAP_PKM2.csv")
# df_UMAP_PKM2 = df_UMAP_PKM2[order(df_UMAP_PKM2$y, decreasing=F),]
# df_UMAP_PKM2$y = factor(df_UMAP_PKM2$y)
# df_UMAP_PKM2$y_alpha = as.numeric(df_UMAP_PKM2$y )
# df_UMAP_PKM2$y_size = df_UMAP_PKM2$y_alpha
# df_UMAP_PKM2$y_size[df_UMAP_PKM2$y_size == 2] = 0.25
# df_UMAP_ALDH1$y_alpha[df_UMAP_PKM2$y_alpha == 1] = 0.05
# 
# sfig1a = ggplot(data=df_UMAP_PKM2, aes(x=UMAP1, y=UMAP2, color=y))+
#   geom_point(size=df_UMAP_PKM2$y_size, shape=16, alpha=df_UMAP_PKM2$y_alpha) +
#   labs(title='PKM2')+
#   scale_color_manual(values=c("#c3d6d9", "#005F73")) +
#   default_plot_theme
# 
# df_UMAP_ALDH1 <- read_csv("UMAP_ALDH1.csv")
# df_UMAP_ALDH1 = df_UMAP_ALDH1[order(df_UMAP_ALDH1$y, decreasing=F),]
# df_UMAP_ALDH1$y = factor(df_UMAP_ALDH1$y)
# df_UMAP_ALDH1$y_alpha = as.numeric(df_UMAP_ALDH1$y )
# df_UMAP_ALDH1$y_size = df_UMAP_ALDH1$y_alpha
# df_UMAP_ALDH1$y_size[df_UMAP_ALDH1$y_size == 2] = 0.25
# df_UMAP_ALDH1$y_alpha[df_UMAP_ALDH1$y_alpha == 1] = 0.05
# 
# sfig1b = ggplot(data=df_UMAP_ALDH1, aes(x=UMAP1, y=UMAP2, color=y))+
#   geom_point(size=df_UMAP_ALDH1$y_size, shape=16, alpha=df_UMAP_ALDH1$y_alpha) +
#   labs(title='ALDH1')+
#   scale_color_manual(values=c("#c3d6d9", "#005F73")) +
#   default_plot_theme
# 
# df_UMAP_VDR <- read_csv("UMAP_VDR.csv")
# df_UMAP_VDR = df_UMAP_VDR[order(df_UMAP_VDR$y, decreasing=F),]
# df_UMAP_VDR$y = factor(df_UMAP_VDR$y)
# df_UMAP_VDR$y_alpha = as.numeric(df_UMAP_VDR$y )
# df_UMAP_VDR$y_size = df_UMAP_VDR$y_alpha
# df_UMAP_VDR$y_size[df_UMAP_VDR$y_size == 2] = 0.25
# df_UMAP_VDR$y_alpha[df_UMAP_VDR$y_alpha == 1] = 0.05
# 
# sfig1c = ggplot(data=df_UMAP_VDR, aes(x=UMAP1, y=UMAP2, color=y))+
#   geom_point(size=df_UMAP_VDR$y_size, shape=16, alpha=df_UMAP_VDR$y_alpha) +
#   labs(title='VDR')+
#   scale_color_manual(values=c("#c3d6d9", "#005F73")) +
#   default_plot_theme
# 
# plot_grid(sfig1a, sfig1b, sfig1c, 
#           labels = c('a', 'b', 'c'), 
#           ncol=3, label_size=8)
# 
# dev.print(pdf, 'al_sfig1_umap.pdf', width = 180/25.4, height = 60/25.4)
