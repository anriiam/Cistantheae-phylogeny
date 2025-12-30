### Running sNMF from VCF files ###

# Anri Chomentowska, 2025
# adapted slightly from https://connor-french.github.io/intro-pop-structure-r/
rm(list = ls())

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")

# load libraries
library(adegenet)
library(vcfR)
library(tidyverse)
library(LEA)


###### Prepping the data ######
setwd("~/Dropbox/Yale/Research/Chapter1/VCF_files")

# read vcf
vcf_file <- "calyptridium_cdsref_min11_0.75.recode.vcf"
vcf <- read.vcfR(vcf_file)

# recode the genotypes according to the .geno system
vcf_geno <- vcf %>% 
  extract_gt_tidy() %>% 
  select(-gt_DP, -gt_CATG, -gt_GT_alleles) %>% 
  mutate(gt1 = str_split_fixed(gt_GT, "/", n = 2)[,1],
         gt2 = str_split_fixed(gt_GT, "/", n = 2)[,2],
         geno_code = case_when(
           # homozygous for reference allele = 0
           gt1 == 0 & gt2 == 0 ~ 0,
           # heterozygous = 1
           gt1 == 0 & gt2 == 1 ~ 1,
           gt1 == 1 & gt2 == 0 ~ 1,
           # homozygous for alternate allele = 2
           gt1 == 1 & gt2 == 1 ~ 2,
           # missing data = 9
           gt1 == "" | gt2 == "" ~ 9
         )) %>% 
  select(-gt_GT, -gt1, -gt2)

# rotate the table so individuals are columns and genotypes are rows
geno_snmf <- vcf_geno %>% 
  pivot_wider(names_from = Indiv, values_from = geno_code) %>% 
  select(-Key)

# take a look at the data
geno_snmf

# look up sample order:
sample_ids <- colnames(geno_snmf)
sample_df <- data.frame(ID = sample_ids)

# write the order of samples to csv
setwd("~/Dropbox/Yale/Research/Chapter1/sNMF")
write.csv(sample_df, "calyp_order_cds.csv", row.names = FALSE) # outside of R, add a 'plot_order' column to match phylogeny order (use numbers) and add a 'pop' column and add description of the subclade

# write the new genotype file
write.table(geno_snmf, 
            "~/Dropbox/Yale/Research/Chapter1/sNMF/calyp_geno_cdsref.geno",
            col.names = FALSE,
            row.names = FALSE,
            sep = "")


###### Running sNMF ######
setwd("~/Dropbox/Yale/Research/Chapter1/sNMF")

# re-read your .geno file and run snmf
data_snmf <- snmf(input.file = "calyp_geno_cdsref.geno",
                  K = 1:20,
                  entropy = TRUE,
                  repetitions = 100,
                  project = "new",
                  alpha = 100
)

# quick plot of cross-entropy estimates and save this
plot(data_snmf, cex = 1.2, col = "lightblue", pch = 19)

# choose the best K
ce <-  cross.entropy(data_snmf, K = 5)
ce

# select the run with the lowest cross-entropy for K = 5
best_run <- which.min(ce)
best_run

# extract the Q-matrix (admixture/ancestry matrix)
q_mat <- LEA::Q(data_snmf, K = 5, run = best_run) 

# add column ("population") names
colnames(q_mat) <- paste0("P", 1:5) #change this to 1:K as appropriate
head(q_mat)


###### Plotting sNMF results######
setwd("~/Dropbox/Yale/Research/Chapter1/sNMF")

# to plot in the right order, read the updated .csv order file
pops <- read.csv("calyp_order_cds.csv")

# convert the Q matrix to a data frame
q_df <- q_mat %>% 
  as_tibble() %>% 
  # optional for plotting
  mutate(individual = pops$ID,
         identity = pops$pop,
         order = pops$plot_order)

# take a look
q_df

# transform long-form to plot proportions of ancestry assignment
q_df_long <- q_df %>% 
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 
q_df_long

# quick plot
q_df_plot <- q_df_long %>% 
  # arrange the data set by the plot order
  arrange(order) %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))
q_df_plot

# nicer plot
q_df_plot %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_viridis_d() +
  labs(fill = "Group") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6, margin = margin(t = -5)),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  )


###### Reload previously run analyses ######
setwd("~/Dropbox/Yale/Research/Chapter1/sNMF")

# load project and find best K
project <- load.snmfProject("srosu_genomeref.snmfProject")

plot(project, cex = 1.2, col = "lightblue", pch = 19)

ce <-  cross.entropy(project, K = 5)
ce

best_run <- which.min(ce)
best_run

# Q matrix
q_mat <- LEA::Q(project, K = 5, run = best_run) 

colnames(q_mat) <- paste0("P", 1:5) #change this to 1:K as appropriate

head(q_mat)

# load order spreadsheet
setwd("~/Dropbox/Yale/Research/Chapter1/sNMF")
pops <- read.csv("srosu_order_genome.csv")

# convert the Q matrix to a data frame
q_df <- q_mat %>% 
  as_tibble() %>% 
  mutate(individual = pops$ID,
         species = pops$pop,
         order = pops$plot_order)
q_df

# transform long-form
q_df_long <- q_df %>% 
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 
q_df_long

# plots
q_df_plot <- q_df_long %>% 
  arrange(order) %>% 
  mutate(individual = forcats::fct_inorder(factor(individual)))
q_df_plot

q_df_plot %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_viridis_d() +
  labs(fill = "Group") +
  theme_minimal() +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6, margin = margin(t = -5)),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  )



