### Cistantheae phylogenetic analyses ###
# Anri Chomentowska, 2025 #

#Set up
rm(list = ls())

set.seed(123)
setwd("~/Dropbox/Yale/Research/Chapter1") #change

#load libraries
library(ape)
library(phytools)
library(geiger)
library(tidyverse)
library(ggtree)
library(tidytree)
library(viridis)
library(mapdata)
library(ggrepel)
library(treeio)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(ggmap)
library(maps)
library(geosphere) 
library(caper)
library(tibble)
library(purrr)
library(broom)
library(coda)
library(phylolm)
library(FactoMineR)
library(factoextra)


###### Time calibrating the tree ####### 
ml.tree <- read.tree("FINAL_v14_trimmed.tree")
plotTree(ml.tree,fsize=0.4)

#find node numbers
ggtree(ml.tree) + 
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 2) +
  geom_tiplab(size=2) +
  geom_treescale(x=0, y=45) +
  xlim(NA, 0.15)

#Nodes: (phemeranthus, rest), (calandrinia+lewisiopsis, rest), (paniculata, rest)
nodes=c(53,55,99)

##Option 1: range, upper Arakaki 2011 / lower timetree.org
age.min=c(31,24.4,2.4)
age.max=c(39.9,30.9,9.9)

##Option 2: no range, just Arakaki 2011
#USE THIS
age.min=c(39.9,30.9,9.9)
age.max=c(39.9,30.9,9.9)

#make sure we have right nodes
plotTree(ml.tree,fsize=0.4)
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(obj$xx[nodes],obj$yy[nodes],pch=21,
       bg=palette()[1:5],cex=2)

#create data frame containing calibration points
calibration<-makeChronosCalib(ml.tree,node=nodes,
                              age.min=age.min,age.max=age.max)
calibration

#fit model using penalized-likelihood
pl.tree<-chronos(ml.tree,calibration=calibration)
pl.tree

par(mfrow = c(1, 1))

plotTree(pl.tree,direction="leftwards",
         xlim=c(40,-10),ftype="i",mar=c(4.1,1.1,0.1,1.1),
         fsize=0.4)
axis(1)
title(xlab="millions of years before present")
abline(v=seq(0,40,by=10),lty="dotted",col="grey")

#adjust lambda
#very little smoothing
#using discrete model rather than model = "correlated"

lambda.tree<-chronos(ml.tree,calibration=calibration, model = "discrete",
                     lambda=0.1)

plotTree(lambda.tree,direction="leftwards",
         xlim=c(40,-10),ftype="i",mar=c(4.1,1.1,0.1,1.1),
         fsize=0.4)
axis(1)
title(xlab="millions of years before present")
abline(v=seq(0,40,by=10),lty="dotted",col="grey")

dev.new(width=5, height=4, unit="in")
plot(pl.tree$edge.length,lambda.tree$edge.length,
     pch=21,bg="grey",cex=1.2,bty="n",
     xlab=expression(paste(lambda,"= 1")),
     ylab=expression(paste(lambda,"= 0.1")))
lines(c(0,40),c(0,40))
legend("topleft","1:1 line",lty="solid",bty="n")
grid()
title(main=expression(paste(
  "Comparison of edge lengths with two different ",
  lambda," values")))

write.tree(lambda.tree, file = "~/Dropbox/Yale/Research/Chapter1/FINAL_v14_time.tree")



####### trait evolution - discrete #######
###load data
trim_tree <- read.tree("FINAL_v14_time_cleanname.tree") #labels cleaned
df_tree <- read.csv("FINAL_v14_2025_LH_bioclim.csv")
View(df_tree)

#let's quickly plot the tree
plotTree(trim_tree,direction="leftwards",
         xlim=c(40,-10),ftype="i",mar=c(4.1,1.1,0.1,1.1),
         fsize=0.4)
axis(1)
title(xlab="millions of years before present")
abline(v=seq(0,40,by=10),lty="dotted",col="grey")

#just take out life history
df_LH <- df_tree %>% select(species, life_history)
df_LH <- df_LH %>%
  column_to_rownames(var = "species")

#reorder the data frame to match the tree tip labels
df_LH <- df_LH[trim_tree$tip.label, , drop = FALSE]
View(df_LH)

#create the named life history vector
LH <- setNames(df_LH$life_history, rownames(df_LH))
head(LH)

#phew, the data is now in a condition to be worked with, so let's plot
plotTree(trim_tree,fsize=0.4,ftype="i",lwd=1,offset=0.5)

cols<-setNames(c("#E41A1C", "#377EB8"),sort(unique(LH)))
tiplabels(pie=to.matrix(LH,sort(unique(LH))),piecol=cols,cex=0.2)
add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8)


###Models
#define and fit models of character evolution
model_er<-fitMk(trim_tree,LH,model="ER",pi="fitzjohn")
#model_sym<-fitMk(trim_tree,LH,model="IR",pi="fitzjohn")
model_ard<-fitMk(trim_tree,LH,model="ARD",pi="fitzjohn")

#anova to average models
model_aov<-anova(model_er,model_ard)

#save the model averages
write.csv(model_aov, "LH_model_anova.csv")

#run ancr
tree_ancr<-ancr(model_aov)
#tree_ancr<-ancr(model_er)
tree_ancr
result <- tree_ancr$ace

#write results
write.csv(result, "LH_result_ER.csv")

#plot ancr
cols<-setNames(c("#E41A1C", "#377EB8"),sort(unique(LH)))
plot(tree_ancr,
     args.plotTree=list(fsize=0.4,offset=0.5),
     args.nodelabels=list(cex=0.35, piecol=cols),
     args.tiplabels=list(cex=0.2, piecol=cols),
     legend=FALSE)
legend(2,50,
       sort(unique(LH)),pch=16,col=cols,
       horiz=FALSE,cex=0.8,bty="n",pt.cex=2,
       y.intersp=1.2)


####### trait evolution - continuous ####### 
#prep
df_tree <- df_tree %>%
  column_to_rownames(var = "species") %>%
  .[trim_tree$tip.label, , drop = FALSE]           # Reorder by tree tip
View(df_tree)

#choose continuous trait to model
df_select <- df_tree %>% dplyr::select(bio14_mean)
cont_trait <- as.matrix(df_select)[,1]

#ML estimation of ancestral states under BM
fit.BM <- phytools::fastAnc(trim_tree, cont_trait, vars=TRUE, CI=TRUE)

obj<-contMap(trim_tree, cont_trait, plot=FALSE)
obj

#choose appropriate color pallete for trait 
custom_colors <- colorRampPalette(c("#377EB8", "gray80", "#E41A1C"))(100)
custom_colors <- colorRampPalette(c("#4D2600", "#D2B48C", "#F5F5DC"))(100)
custom_colors <- rev(viridis::viridis(100))
custom_colors <- viridis::viridis(100)

obj <- setMap(obj, custom_colors)

plot(obj,legend=15,sig=2,fsize=c(0.7,0.9))



####### evolutionary correlation between discrete and continuous traits ####### 
###PhyloANOVA
#create a named vector for continuous trait of choice
trait_vec <- setNames(df_tree$bio11_mean, rownames(df_tree))

#create a named factor for life history (binary trait)
lh_vec <- setNames(df_tree$life_history, rownames(df_tree))

#run phylogenetic ANOVA
phylANOVA(trim_tree, x = lh_vec, y = trait_vec, nsim = 1000)

#define color palette
life_history_colors <- c("Annual" = "#E41A1C", "Perennial" = "#377EB8")

#plot
ggplot(df_tree, aes(x = LH, y = bio11_mean, color = LH)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) + 
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_color_manual(values = life_history_colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Life History",
    y = "Mean Annual Precipitation (bio11)",
    color = "Life History"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )


###PGLS
#move rownames (species names) into a column called species
df_tree_rownames <- df_tree %>%
  rownames_to_column(var = "species")

#make trim_tree nodeless
trim_tree_nodeless <- trim_tree
trim_tree_nodeless$node.label <- NULL

#prep comparative data
comp_data <- comparative.data(
  phy = trim_tree_nodeless,
  data = df_tree_rownames,
  names.col = "species",
  vcv = TRUE,
  na.omit = FALSE
)

#convert life history to numeric (Annual = 0, Perennial = 1)
comp_data$data$LH_numeric <- ifelse(comp_data$data$life_history == "Perennial", 1, 0)

#fit PGLS model (choose one cont variable to test)
pgls_model <- pgls(bio1_mean ~ LH_numeric, data = comp_data)

#check model summary
summary(pgls_model)

#nowfor the loop!
df_tree_rownames$LH_numeric <- ifelse(df_tree_rownames$life_history == "Perennial", 1, 0)

#select just continuous traits
continuous_traits <- df_tree_rownames %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-LH_numeric) %>%  # don't include predictor!!!
  names()

#create empty list to store model outputs
results_list <- list()

#loop through each trait
for (trait in continuous_traits) {
  #prep data
  temp_df <- df_tree_rownames %>%
    dplyr::select(species, LH_numeric, all_of(trait))
  
  #drop NA rows
  temp_df <- na.omit(temp_df)
  
  #match tree to data
  temp_comp <- comparative.data(
    phy = trim_tree_nodeless,
    data = temp_df,
    names.col = "species",
    vcv = TRUE,
    na.omit = FALSE
  )
  
  #build model formula dynamically
  model_formula <- as.formula(paste(trait, "~ LH_numeric"))
  
  #run PGLS with lambda estimated
  model <- tryCatch(
    {
      pgls(model_formula, data = temp_comp, lambda = "ML")
    },
    error = function(e) return(NULL)
  )
  
  #extract results for successful model run
  if (!is.null(model)) {
    model_summary <- summary(model)
    
    results_list[[trait]] <- tibble(
      trait = trait,
      estimate = coef(model)[2],
      std_error = model_summary$coefficients[2, "Std. Error"],
      t_value = model_summary$coefficients[2, "t value"],
      p_value = model_summary$coefficients[2, "Pr(>|t|)"],
      r_squared = model_summary$r.squared,
      lambda = model$param["lambda"]
    )
  }
}

#combine all results into one data frame
results_df <- bind_rows(results_list)

#view results
print(results_df)

#write results
write.csv(results_df, "PGLS_results.csv")

#apply FDR correction to raw p-values
results_df <- results_df %>%
  mutate(
    p_value_fdr = p.adjust(p_value, method = "fdr"),
    significant = p_value_fdr < 0.05
  )

#view top results
head(results_df[order(results_df$p_value_fdr), ])

#save new version
write.csv(results_df, "PGLS_results_with_FDR.csv", row.names = FALSE)


###Threshbayes
#convert to numeric (0 = Annual, 1 = Perennial)
lh_trait <- ifelse(df_tree$life_history == "Perennial", 1, 0)
names(lh_trait) <- rownames(df_tree)

#choose continuous trait
bio11_trait <- df_tree$bio11_mean
names(bio11_trait) <- rownames(df_tree)

#combine discrete and continuous traits into 2-column matrix
trait_matrix <- cbind(lh_trait, bio11_trait)
View(trait_matrix)

#run MCMC
set.seed(111)
tb_result_bio11 <- threshBayes(
  tree = trim_tree_nodeless,
  X = trait_matrix,
  types = c("disc", "cont"),
  ngen = 2e6
)

#plot posterior distribution of correlation
plot(tb_result_bio11)

mcmc_result <- tb_result_bio11$par  # Posterior samples
ngen <- nrow(mcmc_result)  # should be 10,000 if sample every 100 over 1e6
burnin <- round(0.1 * ngen) # 10%

post_burnin <- mcmc_result[(burnin + 1):ngen, "r"]
mean_corr <- mean(post_burnin)
cred_int <- quantile(post_burnin, probs = c(0.025, 0.975))

cat("Mean posterior correlation:", round(mean_corr, 3), "\n")
cat("95% credible interval:", round(cred_int[1], 3), "to", round(cred_int[2], 3), "\n")

#estimate posterior correlation
par(mfrow = c(1, 1))
plot(density(post_burnin, bw = 0.05), 
     xlab = "Evolutionary Correlation (r)", 
     main = "Posterior Density for Correlation (bio11 vs Life History)",
     ylim = c(0, 4))  # <- adjust this to make plot taller

#add lines for mean and zero
abline(v = 0, lty = "dashed", col = "red")
abline(v = mean_corr, col = "blue", lwd = 2)

#extract the correlation values post-burnin
mcmc_r <- tb_result_bio11$par[(burnin + 1):nrow(tb_result_bio11$par), "r"]

#convert to mcmc object
r_chain <- as.mcmc(mcmc_r)

#effective sample size
ess <- effectiveSize(r_chain)
print(ess)

#95% highest posterior density interval
hpd <- HPDinterval(r_chain)
print(hpd)

#a function to automate above:
summarize_threshBayes <- function(tb_result_bio11, trait_name, burnin_prop = 0.1) {
  mcmc_result <- tb_result_bio11$par
  ngen <- nrow(mcmc_result)
  burnin <- round(burnin_prop * ngen)
  post_burnin <- mcmc_result[(burnin + 1):ngen, "r"]
  r_chain <- as.mcmc(post_burnin)
  
  tibble::tibble(
    trait = trait_name,
    mean_corr = mean(post_burnin),
    ci_lower = quantile(post_burnin, 0.025),
    ci_upper = quantile(post_burnin, 0.975),
    ess = effectiveSize(r_chain)[[1]],
    hpd_lower = HPDinterval(r_chain)[1],
    hpd_upper = HPDinterval(r_chain)[2]
  )
}

#write a summary
summary_tb_bio11 <- summarize_threshBayes(tb_result_bio11, trait_name = "bio11_mean")
print(summary_tb_bio11)

#save the summary
write.csv(summary_tb_bio11, "threshBayes_summary_bio11_2mil.csv")


###Phylolm
#phylolm::phylolm() with model = "OU" / both regression and ancestral states

#will need
lh_trait <- ifelse(df_tree$life_history == "Perennial", 1, 0)
names(lh_trait) <- rownames(df_tree)

#choose continuous trait
trait <- "bio_11"
trait_formula_intercept <- as.formula(paste0(trait, " ~ 1"))
trait_formula_regression <- as.formula(paste0(trait, " ~ LH_numeric"))

#make sure species are rownames
df_tree_rownames_true <- df_tree_rownames %>%
  column_to_rownames("species")

##intercept-only models
fit_bm <- phylolm(trait_formula_intercept, data = df_tree_rownames_true, phy = trim_tree_nodeless, model = "BM")
fit_ou <- phylolm(trait_formula_intercept, data = df_tree_rownames_true, phy = trim_tree_nodeless, model = "OUfixedRoot")
#for elevation only,
fit_ou <- phylolm(trait_formula_intercept, data = df_tree_rownames_true, phy = trim_tree_nodeless, model = "OUfixedRoot", upper.bound = 100)

##regression models
fit_bm_reg <- phylolm(trait_formula_regression, data = df_tree_rownames_true, phy = trim_tree_nodeless, model = "BM")
fit_ou_reg <- phylolm(trait_formula_regression, data = df_tree_rownames_true, phy = trim_tree_nodeless, model = "OUfixedRoot")
#for elevation only,
fit_ou_reg <- phylolm(trait_formula_regression, data = df_tree_rownames_true, phy = trim_tree_nodeless, model = "OUfixedRoot", upper.bound = 100)

#calculate pseudo-R²: 1 - (RSS_reg / RSS_null)
rss_null <- sum(residuals(fit_ou)^2)
rss_reg <- sum(residuals(fit_ou_reg)^2)
pseudo_r2 <- 1 - (rss_reg / rss_null)

#ancestral states (OU)
ou_ancestral_states <- fit_ou$anc.reconstruction
trait_combined <- c(
  setNames(df_tree_rownames_true[[trait]], rownames(df_tree_rownames_true)), # tips
  ou_ancestral_states # nodes
)

#plot trait evolution
trait_map <- contMap(trim_tree_nodeless, trait_combined, plot = FALSE)

#choose color palette
#temp
trait_map <- setMap(trait_map, colors = colorRampPalette(c("#377EB8", "gray80", "#E41A1C"))(100))
#moisture
trait_map <- setMap(trait_map, colors <- rev(viridis::viridis(100)))
#elevation
trait_map <- setMap(trait_map, colors = colorRampPalette(c("#4D2600", "#D2B48C", "#F5F5DC"))(100))

plot(trait_map, legend = FALSE, fsize = 0.7, lwd = 4, outline = FALSE)
add.color.bar(
  10, trait_map$cols,
  title = paste(trait),
  lims = trait_map$lims,
  digits = 1,
  prompt = FALSE,
  x = 2, y = 50,
  lwd = 5, fsize = 0.8,
  outline = FALSE
)

tip_states <- lh_trait[trim_tree_nodeless$tip.label]
tip_colors <- ifelse(tip_states == 0, "#E41A1C", "#377EB8")
tiplabels(pch = 21, bg = tip_colors, cex = 0.9)

#save AICs and log-likelihoods results
model_comp_df <- tibble::tibble(
  model = c("BM (intercept)", "OU (intercept)", "BM (regression)", "OU (regression)"),
  AIC = c(AIC(fit_bm), AIC(fit_ou), AIC(fit_bm_reg), AIC(fit_ou_reg)),
  logLik = c(fit_bm$logLik, fit_ou$logLik, fit_bm_reg$logLik, fit_ou_reg$logLik)
)

#save OU regression summary results
reg_ou_summary <- summary(fit_ou_reg)
reg_results_df <- tibble::tibble(
  term = rownames(reg_ou_summary$coefficients),
  estimate = reg_ou_summary$coefficients[, "Estimate"],
  std_error = reg_ou_summary$coefficients[, "StdErr"],
  t_value = reg_ou_summary$coefficients[, "t.value"],
  p_value = reg_ou_summary$coefficients[, "p.value"],
  pseudo_R2 = pseudo_r2
)

#write to .csv files
write.csv(model_comp_df, paste0("phylolm_", trait, "_model_comparison.csv"), row.names = FALSE)
write.csv(reg_results_df, paste0("phylolm_", trait, "_ou_regression_summary.csv"), row.names = FALSE)

#read in compiled phylolm results
phylolm_df <- read.csv("phylolm_summary.csv")

#apply FDR correction
phylolm_df <- phylolm_df %>%
  mutate(p_value_fdr = p.adjust(p_value, method = "fdr"),
         significant = p_value_fdr < 0.05)

#write updated version
write.csv(phylolm_df, "phylolm_results_with_FDR.csv", row.names = FALSE)



####### climate PCA ####### 
#extract and scale climate traits
climate_data <- df_tree_rownames %>%
  dplyr::select(species, where(is.numeric)) %>%
  dplyr::select(-LH_numeric) %>%  # remove numeric LH if it’s there
  column_to_rownames("species")  # set species as rownames
View(climate_data)

#store grouping variable i.e. life history
life_history_vec <- df_tree_rownames$life_history
names(life_history_vec) <- df_tree_rownames$species

#make sure above vector is a factor and matches rownames
group_vec <- factor(life_history_vec[rownames(climate_data)])

#run PCA using FactoMineR
climate_pca <- PCA(climate_data, scale.unit = TRUE, graph = FALSE)

#plot
fviz_pca_ind(climate_pca,
             label = "none",
             geom.ind = "point",
             col.ind = group_vec,
             palette = c("Annual" = "#E41A1C", "Perennial" = "#377EB8"),
             addEllipses = TRUE,
             ellipse.type = "confidence") +
  labs(title = "PCA of Climatic Traits Colored by Life History")

fviz_pca_var(climate_pca,
             col.var = "contrib",
             gradient.cols = c("lightblue", "blue", "darkblue"),
             repel = TRUE) +
  labs(title = "Trait Contributions to PCA Axes")

#tidy version of loadings/correlations with PC axes
loadings_df <- as.data.frame(climate_pca$var$cor) %>%
  rownames_to_column("trait") %>%
  arrange(desc(abs(Dim.1))) %>%
  mutate(top_PC1 = row_number() <= 5) %>%
  arrange(desc(abs(Dim.2))) %>%
  mutate(top_PC2 = row_number() <= 5)

#view top contributing traits to PC1 and/or PC2
loadings_save <- loadings_df %>% filter(top_PC1 | top_PC2)

#write loadings
write.csv(loadings_save, "PCA_results_topPC1PC2.csv")

#contributions of variables to PC1 and PC2
contribution_PC1 <- climate_pca$var$contrib[, 1:2] %>%
  round(2) %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  arrange(desc(Dim.1)) %>%
  head(10)

contribution_PC2 <- climate_pca$var$contrib[, 1:2] %>%
  round(2) %>%
  as.data.frame() %>%
  rownames_to_column("trait") %>%
  arrange(desc(Dim.2)) %>%
  head(10)

#write contributions
write.csv(contribution_PC1, "PCA_results_contribution_PC1.csv")
write.csv(contribution_PC2, "PCA_results_contribution_PC2.csv")


