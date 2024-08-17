
# compute scaffold stuff in Python
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

se <- function(x, na.rm = FALSE) {sd(x, na.rm=na.rm) / sqrt(sum(1*(!is.na(x))))}
zip <- function(...) {mapply(list, ..., SIMPLIFY = FALSE)}


#### Data for Fig 2: Structural diversity ####

df2 = subset(df,  batch_size == 64 & dataset == 'ALDH1' & n_start == 64) # train_cycle >= 0 &
df2 <- subset(df2, architecture %in% c('gcn', 'mlp'))
df2 <- subset(df2, !(acquisition_method == 'Similarity' & architecture == 'gcn'))


write_csv(df2, 'figures/data/fig2.csv')


#### Data for Fig 3: acq function  ####

df3 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & n_start == 64)

write_csv(df3, 'figures/data/fig3.csv')


##### Data for Fig 4: molecular properties ####

df4_ridge = subset(read_csv("figures/data/properties_ridge.csv"), dataset == 'ALDH1' & hit == 1)
df4_ridge$train_cycle = factor(df4_ridge$train_cycle, levels=unique(df4_ridge$train_cycle))
df4_ridge = reshape2::melt(df4_ridge, id.vars = c("train_cycle", "seed", "acquisition_method", 'architecture', 'bias'), variable.name = "metric", measure.vars = 9:13)
df4_ridge$acquisition_method = gsub('random', 'Random', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('similarity', 'Similarity', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('exploitation no retrain', 'Exploit, no retraining', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('exploitation', 'Exploit (best predictions)', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('bald', 'BALD (least mutual information)', df4_ridge$acquisition_method)
df4_ridge$acquisition_method = gsub('exploration', 'Explore (most uncertain)', df4_ridge$acquisition_method)
df4_ridge = subset(df4_ridge, acquisition_method %in% c("Random", 'BALD (least mutual information)', "Exploit (best predictions)"))
df4_ridge = subset(df4_ridge, bias == 'random' & architecture == 'mlp' & seed == 2 & metric %in% c('TPSA', 'MolWt'))

write_csv(df4_ridge, 'figures/data/fig4_de.csv')

###

pattern_occurence <- read_csv("figures/data/pattern_occurence_ALDH1.csv", col_types = cols(...1 = col_skip(), index = col_skip()))
pattern_occurence$total_patterns = pattern_occurence$hit_patterns + pattern_occurence$nonhit_patterns
df4 = subset(df, train_cycle >= 0 & batch_size == 64 & bias == 'Innate' & architecture == 'mlp' & dataset == 'ALDH1' & n_start == 64)

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


## h ##
pattern_enrichment_aldehyde = c()
for (x in zip(df4$total_mols_screened, df4['Aldehyde'][[1]])){
  ntot_mols = x[[1]]
  found = x[[2]]
  expected = subset(pattern_occurence, pattern == 'Aldehyde')$total_patterns*ntot_mols/100000
  enrichment = found/expected
  pattern_enrichment_aldehyde = c(pattern_enrichment_aldehyde, enrichment)
}
df4$pattern_enrichment_aldehyde = pattern_enrichment_aldehyde


## f ##
expected_nr_of_aldehydes_in_bacth = subset(pattern_occurence, pattern == 'Aldehyde')$total_patterns*64/100000
expected_aldehydes_per_hit = subset(pattern_occurence, pattern == 'Aldehyde')$hit_patterns/4986
print(paste0(expected_nr_of_aldehydes_in_bacth, "% of hits, ", expected_aldehydes_per_hit, "% of non-hits"))

aldehyde_hit_enrichment = c()
for (x in zip(df4$hits_discovered, df4['hit_Aldehyde'][[1]])){
  ntot_hits = x[[1]]
  found = x[[2]]
  
  expected = ntot_hits*subset(pattern_occurence, pattern == 'Aldehyde')$hit_patterns/4986
  enrichment = found/expected
  aldehyde_hit_enrichment = c(aldehyde_hit_enrichment, enrichment)
}
df4$aldehyde_hit_enrichment = aldehyde_hit_enrichment


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

write_csv(df4, 'figures/data/fig4_fghijk.csv')


##### Data for Fig 5: low data ####


df5 = subset(df, train_cycle <= 15 & batch_size == 64 & bias == 'Innate')
df5 = subset(df5, acquisition_method %in% c("Similarity", "BALD (least mutual information)"))
df5 = subset(df5, dataset %in% c("PKM2" ,"ALDH1", "VDR"))
df5 = df5[c(2, 32, 34, 35, 36, 38, 41, 153)]

df5$id = paste(df5$architecture, df5$n_start, df5$seed, df5$dataset)

df5_sim = subset(df5, acquisition_method == "Similarity")
df5_bald = subset(df5, acquisition_method == "BALD (least mutual information)")

df5_bald = df5_bald[df5_bald$id %in% df5_sim$id, ]
df5_sim = df5_sim[df5_sim$id %in% df5_bald$id, ]

df5_sim = df5_sim[order(df5_sim$id),]
df5_bald = df5_bald[order(df5_bald$id),]

# difference in enrichment
df5_sim$enrichment = df5_bald$enrichment - df5_sim$enrichment 

write_csv(df5_sim, 'figures/data/exfig1.csv')


####


