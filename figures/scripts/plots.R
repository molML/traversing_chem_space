
library(readr)
library(dplyr)
library(data.table)
library(readr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(data.table)
library(patchwork)
library(scico)



se <- function(x, na.rm = FALSE) {sd(x, na.rm=na.rm) / sqrt(sum(1*(!is.na(x))))}
zip <- function(...) {mapply(list, ..., SIMPLIFY = FALSE)}


default_plot_theme = theme(
  panel.background = element_blank(),
  plot.title = element_text(size=7, hjust = 0.5, face = "plain"),
  axis.text.y = element_text(size=6, face="plain", colour = "#1e3648"),
  axis.text.x = element_text(size=6, face="plain", colour = "#1e3648"),
  axis.title.x = element_text(size=7, face="plain", colour = "#1e3648"),
  axis.title.y = element_text(size=7, face="plain", colour = "#1e3648"),
  axis.ticks = element_line(color="#1e3648", size=0.45),
  axis.line.x.bottom=element_line(color="#1e3648", size=0.45),
  axis.line.y.left=element_line(color="#1e3648", size=0.45),
  legend.key = element_blank(),
  legend.position = 'None',
  legend.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_text(size=7),
  legend.text = element_text(size=7),
  legend.spacing.y = unit(0., 'cm'),
  legend.key.size = unit(0.25, 'cm'),
  legend.key.width = unit(0.5, 'cm'))


custom_colours = c('#653f88', '#926db4', '#c3a3e0',
                   '#3e5a6a', '#5c8095', '#85b8d6',
                   '#2e5950', '#3d7368', '#60ae9e',
                   '#7b5435', '#c28a5e', '#efc57b',
                   '#766d40', '#a89b5a', '#dccd7e')


custom_colours = c("#5c8095", "#60ae9e", "#85b8d6", "#c28a5e", "#bbbbbb", "#efc57b")

umap_cols = c("#bbbbbb", "#5c8095", "#dccd7e")
bias_colours = c("#5c8095", "#efc57b", "#c28a5e")
# col_gradient = c("#3d7368", "#5c8095", "#85b8d6", "#bbbbbb", "#efc57b", "#c28a5e")

col_gradient = c("#7b5435", "#c28a5e", "#efc57b", "#85b8d6", "#5c8095", "#3e5a6a")
# col_gradient = c("#3e5a6a", "#5c8095", "#85b8d6", "#60ae9e", "#3d7368", "#2e5950")


acq_funcs = c("BALD (least mutual information)", 
              "Exploit (best predictions)", 
              "Exploit, no retraining", 
              "Explore (most uncertain)", 
              "Random",
              "Similarity")

acq_cols = list(custom_colours, acq_funcs)

setwd("~/Dropbox/PycharmProjects/traversing_chem_space")

source('figures/scripts/data_processing.R')
source('figures/scripts/fig2.R')
source('figures/scripts/fig3.R')
source('figures/scripts/fig4.R')
source('figures/scripts/exfig1.R')
source('figures/scripts/supplementary.R')


