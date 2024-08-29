#--------------------------------------------------------------#
# Prepare workspace
#--------------------------------------------------------------#

# Load libraries
library(ggExtra)
library(ggridges)
library(ggrepel)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
library(hrbrthemes)
library(ggpubfigs)
library(seqinr)
library(ggseqlogo)
library(cowplot)
library(stringr)
library(treemapify)
library(mgcv)
library(shadowtext)
library(ggtree)
library(treeio)
library(ape)
library(phytools)
library(ggtree)


# Set wd
setwd("~/Dropbox/Publications/starfleet/analyses/Figure3/tyr")

# define order of classes, based on JGI mycocosm tree
classOrder <- c("Pucciniomycotina", "Pucciniomycetes", "Ustilaginomycotina", "Ustilaginomycetes", "Agaricomycetes", "Dacrymycetes", "Tremellomycetes", "Wallemiomycetes", "Pezizomycetes", "Orbiliomycetes", "Eurotiomycetes", "Dothideomycetes", "Lecanoromycetes", "Leotiomycetes", "Sordariomycetes", "Xylonomycetes", "Saccharomycotina", "Taphrinomycotina", "Mucoromycota", "Glomeromycotina", "Mortierellomycotina", "Mucormycotina", "Zoopagomycota", "Zoopagomycotina", "Entomophthoromycotina", "Kickxellomycotina", "Chytridiomycota", "Blastocladiomycota", "Chytridiomycetes", "Monoblepharidomycetes", "Neocallimastigomycetes", "Microsporidia", "Cryptomycota")

#--------------------------------------------------------------#
# Figure 2
#--------------------------------------------------------------#

#--------------------------------------------------------------#
# Starfish annotate, insert and flank benchmarking
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure2/benchmarkingStarfish//")

# precision = TP / TP + FP
# recall = TP / TP + FN
# TP: starfish prediction matches manual annotation
# FP: starfish prediction does not match manual annotation (>1kb difference for boundaries, >1bp difference for direct repeat)
# FN: no feature predicted by starfish

starfishBenchmark <- read.table("starfishBenchmarkStats.txt", sep = "\t", header = T)
starfishBenchmark$module <- factor(starfishBenchmark$module, levels = c("captain", "boundary", "DR"))
ggplot(starfishBenchmark) +
  geom_col(aes(x = benchmark, y = F1_score), fill = "darkgray", width = 0.75) + 
  geom_text(aes(x = benchmark, y = as.double(F1_score) + 0.025, label = F1_score), angle = 0, size = 4) +
  facet_grid(.~ module, scales = "free",  space = "free") +
  ggpubfigs::theme_grey() + 
  theme(text  = element_text(size = 15, family = "Helvetica"),
        axis.text.x=element_text(angle=45)) +
  labs(title = "Starfish performance on 42 manually annotated elements", y = "F1 score", x = "Starfish module")


#--------------------------------------------------------------#
# FIGURE 3
#--------------------------------------------------------------#

#--------------------------------------------------------------#
# Tyr family pairwise identities
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure3/tyr/")

famPID <- read.table("funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.pid_byfam.txt", sep = "\t", header = T)

ggplot(famPID %>%
         mutate(percentID = pid * 100) %>%
         select(comparisonType, percentID)) +
  geom_violin(aes(x = fct_rev(comparisonType), y = percentID), fill = "darkgray") + 
  geom_boxplot(aes(x =  fct_rev(comparisonType), y = percentID), fill = "white", width = 0.2, outlier.size = 0, outlier.fill = NA, outlier.color = NA) + 
  ggpubfigs::theme_grey() + 
  theme(text  = element_text(size = 18, family = "Helvetica")) +
  scale_x_discrete(labels = c('within family', 'between family')) +
  labs(title = "Pairwise seqID of YR refs, from p1-512 alignment", y = "Pairwise sequence identity (%)", x = "Comparison type")

famPID %>%
  group_by(comparisonType) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(min = min,
                          median = median,
                          mean = mean,
                          stdev = sd,
                          q25 = ~quantile(., 0.25),
                          q75 = ~quantile(., 0.75),
                          max = max))) 

#--------------------------------------------------------------#
# Tyr taxonomic distribution
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure3/tyr/")

# read in data
tyrCounts <- read.table("tyrsPerGenome.txt", header = T, sep = "\t")
tyrCounts$taxonID <- factor(tyrCounts$taxonID, levels = classOrder)

# traditional boxplot of YR counts per genome
ggplot(tyrCounts %>%
         select(taxonID, tyrCount) %>%
         filter(taxonID %in% classOrder)) +
  geom_jitter(aes(x = fct_rev(taxonID), y = tyrCount), width = 0.25, height = 0.25, shape = 21, fill = "gray", color = "#3B3B3B", alpha = 0.4, size = 2) +
  geom_boxplot(aes(x =  fct_rev(taxonID), y = tyrCount), fill = "white", width = 0.75, size = 0.6, outlier.size = 0, outlier.fill = NA, outlier.color = NA) + 
  ggpubfigs::theme_grey() + 
  labs(title = "YR taxonomic distribution", x = "Class ID", y = "YR count per genome") +
  ylim(0,30) +
  coord_flip() +
  theme(text  = element_text(size = 6, family = "Helvetica"),
        axis.text.y = element_text(hjust=0))
  
print(n = 30, tyrCounts %>%
  select(taxonID, tyrCount) %>%
  filter(taxonID %in% classOrder) %>%
  group_by(taxonID) %>%
  summarize(tyrTotal = sum(tyrCount)))

# number of starships/genome by number of genomes per species
# in all unamibugous species with at least 1 annotated element
# looks like element annotation also depends on other things, like quality

tyrSpecies <- read.table("tyrsPerSpecies.txt", sep = "\t", header = T)
elementSpecies <- tyrSpecies %>%
  group_by(taxonID) %>%
  tally(elementCount) %>%
  dplyr::filter(n > 0) %>%
  dplyr::filter(!grepl('sp\\.', taxonID)) %>%
  select(taxonID)
ggplot(tyrSpecies %>%
         filter(taxonID %in% elementSpecies$taxonID)) +
  geom_jitter(aes(x = taxonCount, y = elementCount), shape = 21, fill = "gray", color = "#3B3B3B", alpha = 0.4, size = 2) +
  ggpubfigs::theme_grey() + 
  labs(title = "elements per genome as a funciton of genomes per species", x = "Genomes per genus", y = "Elements per genome") +
  theme(text  = element_text(size = 6, family = "Helvetica"),
        axis.text.y = element_text(hjust=0))

#--------------------------------------------------------------#
# Tyr family taxonomic class enrichment
#--------------------------------------------------------------#

setwd("../../Figure3/tyrTaxonomicDistribution/")

# contingency table of the form:
#	                focal taxon	not focal taxon
# focal tyr family	p11				p10
# not focal tyr   	p01				p00

tyrContingency <- read.table("tyr8235_seq2fam_contingencyTable.txt", sep = "\t", header = T)

tyrContingencyData <- tyrContingency %>%
  select(p11,p10,p01,p00)

# fisher test on each row function
row_fisher <- function(row, alt = 'greater', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           p_val = f$p.value,
           or = f$estimate[[1]],
           or_ll = f$conf.int[1],
           or_ul = f$conf.int[2]))
}

# run row_fisher and add to original df
tyrFisherP <- data.frame(t(apply(tyrContingencyData, 1, row_fisher)))
tyrContingency['p_val']= tyrFisherP['p_val']

# correct for multiple testing
p_val_BH <- as.data.frame(p.adjust(tyrContingency$p_val, method = "BH")) %>%
  rename(p_val_BH = 1)
tyrContingency['p_val_BH']= p_val_BH['p_val_BH']

# add dummy rows to represent taxa that are missing YRs (can't think of a nicer way)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam01", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam02", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam03", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam04", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam05", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam06", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam07", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam08", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam09", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam10", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam11", "Pucciniomycetes", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam01", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam02", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam03", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam04", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam05", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam06", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam07", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam08", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam09", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam10", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam11", "Saccharomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam01", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam02", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam03", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam04", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam05", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam06", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam07", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam08", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam09", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam10", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam11", "Taphrinomycotina", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam01", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam02", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam03", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam04", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam05", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam06", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam07", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam08", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam09", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam10", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam11", "Mucoromycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam01", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam02", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam03", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam04", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam05", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam06", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam07", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam08", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam09", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam10", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam11", "Zoopagomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam01", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam02", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam03", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam04", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam05", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam06", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam07", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam08", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam09", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam10", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)
tyrContingency[nrow(tyrContingency) + 1,] <- list("fam11", "Chytridiomycota", 0, 0, 0, 0, 0, 0, 0, 1, 1)

# reorder taxonomic class
tyrContingency$taxonID <- factor(tyrContingency$taxonID, levels = rev(classOrder))

# calculate heatmap of prevalence
# draw the tiles twice with opposite alphas to ensure the red outlines are drawn over adjacent outlines
ggplot(tyrContingency, aes(tyrFam, taxonID)) +
  geom_tile(aes(fill = taxonPrevalence, color=ifelse(p_val_BH<=0.05, '#ff264e', '#F0F0F0'), alpha=ifelse(p_val_BH<=0.05, 0, 1))) +
  geom_tile(aes(fill = taxonPrevalence, color=ifelse(p_val_BH<=0.05, '#ff264e', 'NA'), alpha=ifelse(p_val_BH<=0.05, 1, 0)), linewidth = 1) +
  geom_shadowtext(aes(label=ifelse(p11 > 0, p11, element_blank()), color = "white"), size = 3) +
  scale_color_identity() +
  coord_equal() + 
  guides(alpha="none") +
  scale_fill_gradient(low = "white", high = "black", na.value = "white") +
  theme_minimal() +
  geom_text(aes(x = 12, y = taxonID, label = p11 + p01), angle = 0, size = 4) +
  theme(axis.text.x=element_text(angle=45,vjust = 0.75),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='#eeeeee')) +
  labs(title = "YR taxonomic enrichment (8235 set)", y = "Class ID", x = "YR family")

write_tsv(tyrContingency, file = "tyr8235_seq2fam_contingencyTable_pvalues.txt", quote = "none")


#--------------------------------------------------------------#
# Tyr cladeily sequence logos
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure3/famLogos/")

clade01 <- read.fasta("clade01.p1-512.clipkit", as.string = TRUE, seqonly = TRUE)
clade02 <- read.fasta("clade02.p1-512.clipkit", as.string = TRUE, seqonly = TRUE)
clade03 <- read.fasta("clade03.p1-512.clipkit", as.string = TRUE, seqonly = TRUE)

# extract Core Binding (CB) domain 
clade01.CB <- tibble(CB = character())
for (i in 1:length(clade01)) {
  seq <- clade01[i] %>%
    str_sub(76,166)
  clade01.CB <- add_row(clade01.CB, CB = seq)
}

# extract Catalytic (CAT) domain 
clade01.CAT <- tibble(CAT = character())
for (i in 1:length(clade01)) {
  seq <- clade01[i] %>%
    str_sub(198,469)
  clade01.CAT <- add_row(clade01.CAT, CAT = seq)
}
clade01.CAT.logo <- ggseqlogo(clade01.CAT, seq_type = "aa", col_scheme = "chemistry") +
  theme(axis.text.x = element_blank(), legend.position = "none")

clade02.CAT <- tibble(CAT = character())
for (i in 1:length(clade02)) {
  seq <- clade02[i] %>%
    str_sub(198,469)
  clade02.CAT <- add_row(clade02.CAT, CAT = seq)
}
clade02.CAT.logo <- ggseqlogo(clade02.CAT, seq_type = "aa", col_scheme = "chemistry") +
  theme(axis.text.x = element_blank(), legend.position = "none")

clade03.CAT <- tibble(CAT = character())
for (i in 1:length(clade03)) {
  seq <- clade03[i] %>%
    str_sub(198,469)
  clade03.CAT <- add_row(clade03.CAT, CAT = seq)
}
clade03.CAT.logo <- ggseqlogo(clade03.CAT, seq_type = "aa", col_scheme = "chemistry") +
  theme(legend.position = "none")

plot_grid(clade01.CAT.logo, clade02.CAT.logo, clade03.CAT.logo, ncol = 1, align = 'v')

# extract all 6 active sites
clade01.active <- tibble(active = character())
for (i in 1:length(clade01)) {
  activeSite1 <- clade01[i] %>%
    str_sub(247,247)
  activeSite2 <- clade01[i] %>%
    str_sub(303,303)
  activeSite3 <- clade01[i] %>%
    str_sub(428,428)
  activeSite4 <- clade01[i] %>%
    str_sub(431,431)
  activeSite5 <- clade01[i] %>%
    str_sub(456,456)
  activeSite6 <- clade01[i] %>%
    str_sub(466,466)
  allSites <- paste0(activeSite1, activeSite2, activeSite3, activeSite4, activeSite5, activeSite6)
  clade01.active <- add_row(clade01.active, active = allSites)
}
clade01.active.logo <- ggseqlogo(clade01.active, seq_type = "aa", col_scheme = "chemistry") +
  theme(axis.text.x = element_blank(), legend.position = "none")

clade02.active <- tibble(active = character())
for (i in 1:length(clade02)) {
  activeSite1 <- clade02[i] %>%
    str_sub(247,247)
  activeSite2 <- clade02[i] %>%
    str_sub(303,303)
  activeSite3 <- clade02[i] %>%
    str_sub(428,428)
  activeSite4 <- clade02[i] %>%
    str_sub(431,431)
  activeSite5 <- clade02[i] %>%
    str_sub(456,456)
  activeSite6 <- clade02[i] %>%
    str_sub(466,466)
  allSites <- paste0(activeSite1, activeSite2, activeSite3, activeSite4, activeSite5, activeSite6)
  clade02.active <- add_row(clade02.active, active = allSites)
}
clade02.active.logo <- ggseqlogo(clade02.active, seq_type = "aa", col_scheme = "chemistry") +
  theme(axis.text.x = element_blank(), legend.position = "none")

clade03.active <- tibble(active = character())
for (i in 1:length(clade03)) {
  activeSite1 <- clade03[i] %>%
    str_sub(247,247)
  activeSite2 <- clade03[i] %>%
    str_sub(303,303)
  activeSite3 <- clade03[i] %>%
    str_sub(428,428)
  activeSite4 <- clade03[i] %>%
    str_sub(431,431)
  activeSite5 <- clade03[i] %>%
    str_sub(456,456)
  activeSite6 <- clade03[i] %>%
    str_sub(466,466)
  allSites <- paste0(activeSite1, activeSite2, activeSite3, activeSite4, activeSite5, activeSite6)
  clade03.active <- add_row(clade03.active, active = allSites)
}
clade03.active.logo <- ggseqlogo(clade03.active, seq_type = "aa", col_scheme = "chemistry") +
  theme(axis.text.x = element_blank(), legend.position = "none")

plot_grid(clade01.active.logo, clade02.active.logo, clade03.active.logo, ncol = 1, align = 'v')

# extract all 3 patches and 2 boxes
patchBox <- tibble(famID = character(), regionID = character(), seq = character())

for (i in 1:length(clade01)) {
  patch1 <- clade01[i] %>%
    str_sub(205,214)
  box1 <- clade01[i] %>%
    str_sub(237, 280)
  patch2 <- clade01[i] %>%
    str_sub(302,318)
  patch3 <- clade01[i] %>%
    str_sub(381,386)
  box2 <- clade01[i] %>%
    str_sub(426,468)
  patchBox <- add_row(patchBox, seq = patch1, famID = "clade01", regionID = "patch1")
  patchBox <- add_row(patchBox, seq = patch2, famID = "clade01", regionID = "patch2")
  patchBox <- add_row(patchBox, seq = patch3, famID = "clade01", regionID = "patch3")
  patchBox <- add_row(patchBox, seq = box1, famID = "clade01", regionID = "box1")
  patchBox <- add_row(patchBox, seq = box2, famID = "clade01", regionID = "box2")
}

for (i in 1:length(clade02)) {
  patch1 <- clade02[i] %>%
    str_sub(205,214)
  box1 <- clade02[i] %>%
    str_sub(237, 280)
  patch2 <- clade02[i] %>%
    str_sub(302,318)
  patch3 <- clade02[i] %>%
    str_sub(381,386)
  box2 <- clade02[i] %>%
    str_sub(426,468)
  patchBox <- add_row(patchBox, seq = patch1, famID = "clade02", regionID = "patch1")
  patchBox <- add_row(patchBox, seq = patch2, famID = "clade02", regionID = "patch2")
  patchBox <- add_row(patchBox, seq = patch3, famID = "clade02", regionID = "patch3")
  patchBox <- add_row(patchBox, seq = box1, famID = "clade02", regionID = "box1")
  patchBox <- add_row(patchBox, seq = box2, famID = "clade02", regionID = "box2")
}

for (i in 1:length(clade03)) {
  patch1 <- clade03[i] %>%
    str_sub(205,214)
  box1 <- clade03[i] %>%
    str_sub(237, 280)
  patch2 <- clade03[i] %>%
    str_sub(302,318)
  patch3 <- clade03[i] %>%
    str_sub(381,386)
  box2 <- clade03[i] %>%
    str_sub(426,468)
  patchBox <- add_row(patchBox, seq = patch1, famID = "clade03", regionID = "patch1")
  patchBox <- add_row(patchBox, seq = patch2, famID = "clade03", regionID = "patch2")
  patchBox <- add_row(patchBox, seq = patch3, famID = "clade03", regionID = "patch3")
  patchBox <- add_row(patchBox, seq = box1, famID = "clade03", regionID = "box1")
  patchBox <- add_row(patchBox, seq = box2, famID = "clade03", regionID = "box2")
}

# patch 1 logos
clade01.patch1 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade01", regionID == "patch1") %>%
              select(seq)) +
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade02.patch1 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade02", regionID == "patch1") %>%
              select(seq)) +
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade03.patch1 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade03", regionID == "patch1") %>%
              select(seq)) +
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

# box 1 logos
clade01.box1 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade01", regionID == "box1") %>%
              select(seq)) +
  annotate('text', x=11, y = 3.5, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade02.box1 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade02", regionID == "box1") %>%
              select(seq)) +
  annotate('text', x=11, y = 3.5, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade03.box1 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade03", regionID == "box1") %>%
              select(seq)) +
  annotate('text', x=11, y = 3.5, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

# patch 2 logos
clade01.patch2 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade01", regionID == "patch2") %>%
              select(seq)) +
  annotate('text', x=2, y = 4, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade02.patch2 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade02", regionID == "patch2") %>%
              select(seq)) +
  annotate('text', x=2, y = 4, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade03.patch2 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade03", regionID == "patch2") %>%
              select(seq)) +
  annotate('text', x=2, y = 4, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

# patch 3 logos
clade01.patch3 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade01", regionID == "patch3") %>%
              select(seq)) +
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade02.patch3 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade02", regionID == "patch3") %>%
              select(seq)) +
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade03.patch3 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade03", regionID == "patch3") %>%
              select(seq)) +
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

# box 2 logos
clade01.box2 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade01", regionID == "box2") %>%
              select(seq)) +
  annotate('text', x=3, y = 4, label="*", size = 18) + 
  annotate('text', x=6, y = 4, label="*", size = 18) + 
  annotate('text', x=31, y = 4, label="*", size = 18) + 
  annotate('text', x=41, y = 4, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade02.box2 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade02", regionID == "box2") %>%
              select(seq)) +
  annotate('text', x=3, y = 4, label="*", size = 18) + 
  annotate('text', x=6, y = 4, label="*", size = 18) + 
  annotate('text', x=31, y = 4, label="*", size = 18) + 
  annotate('text', x=41, y = 4, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

clade03.box2 <- ggplot() + 
  geom_logo(patchBox %>%
              filter(famID == "clade03", regionID == "box2") %>%
              select(seq)) +
  annotate('text', x=3, y = 4, label="*", size = 18) + 
  annotate('text', x=6, y = 4, label="*", size = 18) + 
  annotate('text', x=31, y = 4, label="*", size = 18) + 
  annotate('text', x=41, y = 4, label="*", size = 18) + 
  theme_logo() +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none", axis.title.y = element_blank())

plot_grid(clade01.patch1, clade01.box1, clade01.patch2, clade01.patch3, clade01.box2,
          clade02.patch1, clade02.box1, clade02.patch2, clade02.patch3, clade02.box2,
          clade03.patch1, clade03.box1, clade03.patch2, clade03.patch3, clade03.box2,
          rel_widths = c(15,49,23, 11, 48), nrow = 3, ncol = 5, align = 'v')

plot_grid(clade01.box1, clade01.patch2, clade01.box2,
          clade02.box1, clade02.patch2, clade02.box2,
          clade03.box1, clade03.patch2, clade03.box2,
          rel_widths = c(49,23,48), nrow = 3, ncol = 3, align = 'v')

#--------------------------------------------------------------#
# FIGURE 4
#--------------------------------------------------------------#

#--------------------------------------------------------------#
# Cargo OG annotation frequency
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure4/cargoAnn/")

ogCOG <- read.table("ref349_genes_cluster.ann", sep = "\t", header = T)

# define cog2info vectors

cog2cat <- c("d"="Cellular processes and signaling",
"m"="Cellular processes and signaling",
"n"="Cellular processes and signaling",
"o"="Cellular processes and signaling",
"t"="Cellular processes and signaling",
"u"="Cellular processes and signaling",
"v"="Cellular processes and signaling",
"w"="Cellular processes and signaling",
"y"="Cellular processes and signaling",
"z"="Cellular processes and signaling",
"a"="Information storage and processing",
"b"="Information storage and processing",
"j"="Information storage and processing",
"k"="Information storage and processing",
"l"="Information storage and processing",
"c"="Metabolism",
"e"="Metabolism",
"f"="Metabolism",
"g"="Metabolism",
"h"="Metabolism",
"i"="Metabolism",
"p"="Metabolism",
"q"="Metabolism",
"r"="Poorly characterized",
"s"="Poorly characterized")

ggplot(ogCOG %>%
         dplyr::count(cogID) %>%
         dplyr::mutate(category = recode(cogID, !!!cog2cat)),
       aes(area = n, label = cogID, subgroup = category)) +
  geom_treemap(fill = "lightgray", size = 3, alpha = 0.25, colour = "#585858") +
  geom_treemap_text( colour = "#585858",
                     place = "centre",
                     size = 15,
                     grow = TRUE) +
  geom_treemap_subgroup_border(colour = "#585858", size = 5) +
  geom_treemap_subgroup_text(place = "centre", grow = TRUE,
                             alpha = 0.75, colour = "black",
                             fontface = "italic")


#--------------------------------------------------------------#
# Captain patristic vs cargo jaccard distances 
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure4/starshipDist/")

# calculate avg within fam, group and hap patristic distance
pairwisePatristic <- read.table("all_superfam_captains.kpicsg.lengthsOnly.clipkit.patristic_byfam.txt", sep = "\t", header = T)
withinFamMean <- pairwisePatristic %>%
  filter(famType == "withinFam") %>%
  dplyr::summarize(withinFamMean = mean(patristicDist))

withinGroupMean <- pairwisePatristic %>%
  filter(groupType == "withinGroup") %>%
  dplyr::summarize(withinGroupMean = mean(patristicDist))

# plot patristic distance vs jaccard distance in cargo eggnogOG
pairwiseDists <- read.table("ref349_merged_patristic_jaccard_distances.txt", sep = "\t", header = T)

# convert jaccard dist to jaccard similarity
distPlot <- ggplot(pairwiseDists %>%
                     mutate(cargoSim = 1 - cargoDist),
                   aes(x = patristicDist, y = cargoSim)) +
  geom_point(shape = 21, fill = "darkgray", color = "#3B3B3B", alpha = 0.6, size = 2) +
  geom_smooth(method = "gam", color = "red", linewidth = 1.5) +
  ggpubfigs::theme_grey() + 
  geom_vline(xintercept = withinFamMean[1,1], linetype = 2, color = "black", linewidth = 1 ) +
  geom_vline(xintercept = withinGroupMean[1,1], linetype = 2, color = "black", linewidth = 1 ) +
  theme(text  = element_text(size = 18, family = "Helvetica"),
        panel.background = element_rect(fill = 'white', color = 'white')) +
  labs(title = "Captain relatedness vs. cargo similarity", x = "Pairwise patristic distance between captains", y = "Pairwise cargo Jaccard similarity")
ggMarginal(distPlot, type="density", fill = "gray")


# calculate a generalized additive model
# https://m-clark.github.io/generalized-additive-models/application.html
pairwiseDists.mod_gam1 <- gam(cargoDist ~ s(patristicDist, bs = "cs"), data = pairwiseDists)
summary(pairwiseDists.mod_gam1)

# calculate summary stats of cargo similarity within groups and families
# we have a lot of rows with missing values because while all captains are in pairwisePatristic,
# only the captains belonging to the 349 reference elements are in pairwiseDists
pairwiseCombined <- pairwisePatristic %>%
  select(queryID, targetID, groupType, famType) %>% 
  left_join(pairwiseDists, by = c("queryID" = "queryID", "targetID" = "subjectID")) %>%
  mutate(cargoSim = 1 - cargoDist) %>%
  remove_missing() 

# by YR group
pairwiseCombined %>%
  group_by(groupType) %>%
  summarize(median=median(cargoSim, na.rm=TRUE),
            mean=mean(cargoSim, na.rm=TRUE),
            standard_deviation=sd(cargoSim, na.rm=TRUE))

# by YR family
pairwiseCombined %>%
  group_by(famType) %>%
  summarize(median=median(cargoSim, na.rm=TRUE),
            mean=mean(cargoSim, na.rm=TRUE),
            standard_deviation=sd(cargoSim, na.rm=TRUE))

# plot patristic distance vs percent alignable sequence (nucmer -mum -l 1000 -i 50)
pairwiseAlignable <- read.table("ref349_merged_patristic_nucmer_alignable.txt", sep = "\t", header = T)

distPlot <- ggplot(pairwiseAlignable,
                   aes(x = patristicDist, y = percAlignable)) +
  geom_point(shape = 21, fill = "darkgray", color = "#3B3B3B", alpha = 0.6, size = 2) +
  geom_smooth(method = "gam", color = "red", linewidth = 1.5) +
  ggpubfigs::theme_grey() + 
  geom_vline(xintercept = withinFamMean[1,1], linetype = 2, color = "black", linewidth = 1 ) +
  geom_vline(xintercept = withinGroupMean[1,1], linetype = 2, color = "black", linewidth = 1 ) +
  theme(text  = element_text(size = 18, family = "Helvetica"),
        panel.background = element_rect(fill = 'white', color = 'white')) +
  labs(title = "Captain relatedness vs. percent shared sequence (-mum -l 1000 -i 50)", x = "Pairwise patristic distance between captains", y = "Percent shared sequence")
ggMarginal(distPlot, type="density", fill = "gray")

# https://m-clark.github.io/generalized-additive-models/application.html
pairwiseAlignable.mod_gam1 <- gam(percAlignable ~ s(patristicDist, bs = "cs"), data = pairwiseAlignable)
summary(pairwiseAlignable.mod_gam1)

# calculate summary stats of perc alignable within groups and families
# we have a lot of rows with missing values because while all captains are in pairwisePatristic,
# only the captains belonging to the 349 reference elements are in pairwiseDists
pairwiseCombined <- pairwisePatristic %>%
  select(queryID, targetID, groupType, famType) %>% 
  left_join(pairwiseAlignable, by = c("queryID" = "queryID", "targetID" = "subjectID")) %>%
  remove_missing() 

# by YR group
pairwiseCombined %>%
  group_by(groupType) %>%
  summarize(median=median(percAlignable, na.rm=TRUE),
            mean=mean(percAlignable, na.rm=TRUE),
            standard_deviation=sd(percAlignable, na.rm=TRUE))

# by YR family
pairwiseCombined %>%
  group_by(famType) %>%
  summarize(median=median(percAlignable, na.rm=TRUE),
            mean=mean(percAlignable, na.rm=TRUE),
            standard_deviation=sd(percAlignable, na.rm=TRUE))


#--------------------------------------------------------------#
# Starship lengths, by fam
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure4/starshipLengths/")

# don't show starships classified at "clade" level only
starshipLengths <- read.table("mycodb.349refShips.lengths.fam.txt", sep = "\t", header = T) %>%
  dplyr::filter(!grepl('clade', famID))
  
# careful: you can't get a sense of the absolute differences in density because it appears to be standardized across groups
ggplot(starshipLengths) + 
  stat_density_ridges(aes(x = length, y = famID), alpha = 0.9, from = 15000, to = 600000) +
  theme_ridges() + 
  theme(text  = element_text(size = 16, family = "Helvetica"), axis.text.x = element_text(angle=45,hjust=0.95,vjust=1)) +
  labs(title = "Lengths of 349 reference Starships", y = "family ID", x = "length (kilobases)") +
  scale_x_continuous(labels=c("0" = "0", "200000" = "200","400000" = "400", "600000" = "600"))


# go for boxplot and violin plots instead
ggplot(starshipLengths) + 
  geom_violin(aes(x = length, y = famID), fill = "darkgray") +
  geom_boxplot(aes(x = length, y = famID), fill = "white", width = 0.55) +
  geom_text(data = starshipLengths %>%
              group_by(famID) %>%
              tally(),
              aes(x = -15000, y = famID, label = n), angle = 0, size = 4) +
  ggpubfigs::theme_grey() + 
  theme(text  = element_text(size = 16, family = "Helvetica"), axis.text.x = element_text(angle=45,hjust=0.95,vjust=1)) +
  labs(title = "Lengths of 349 reference Starships", y = "family ID", x = "length (kilobases)") +
  scale_x_continuous(labels=c("0" = "0", "100000" = "100", "200000" = "200","300000" = "300","400000" = "400", "500000" = "500"))

# some summary stats
starshipLengths %>%
  group_by(famID) %>%
  summarize(median = median(length, na.rm = TRUE))


#--------------------------------------------------------------#
# FIGURE 5
#--------------------------------------------------------------#

#--------------------------------------------------------------#
# Starship insertion site features
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure5/insertionSites/")

# all categories of interest
siteCats <- read.table("222emptySites_lt30_all_categories.txt", sep = "\t", header = T) %>%
  filter(familyID != "clade01" & familyID != "clade03" & familyID != "fam11")

siteCats$category <- factor(siteCats$category, levels = c("other", "predicted gene", "AT rich", "5s rDNA"))
ggplot(siteCats) +
  geom_bar(aes(x = familyID, fill = category)) + 
  ggpubfigs::theme_grey() + 
  theme(text  = element_text(size = 15, family = "Helvetica"),
        axis.text.x=element_text(angle=45)) +
  scale_fill_manual(values = c("darkgray", "#8da0cb", "#66c2a5", "#fc8d62")) +
  labs(title = "Features of 222 high quality insertion sites <=30bp from the 349 reference elements", y = "Number of insertion sites", x = "Starship family")
  
#--------------------------------------------------------------#
# Direct repeats mapped on the tree of manually verified caps
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure5/directRepeats//")

# read in captain tree (116 manually verified elements spanning captain diversity)
captainTree <- unroot(as.phylo(read.iqtree("116refCaptains.kpicsg.clipkit.treefile")))

# read in 5s rDNA and AT-rich associations
captain5s <- read.table("5s_captainIDs.ids", sep = "\t", header = T)
ATrich  <- read.table("lowGC_captainIDs.ids", sep = "\t", header = T)

# midpoint root the tree
outgroup <- c("laspus2_tyr7798", "fusvir6_tyr6937", "talcel1_GAM38547.1", "asppar10_tyr2350", "fuscul4_tyr5444", "aspluc13_tyr1861", "aspfum13_tyr1277", "colfru6_tyr4318", "orboli11_KAF3111326.1", "fussol1_232585", "botdew1_KAF7931899.1", "macpha1_749149")
captainTree.rooted <- midpoint.root(captainTree)

# collapse all nodes with SH-aLRT < 80% and UFboot < 95%
# captainTree.rooted <- midpoint.root(as.polytomy(captainTree, feature='node.label', fun=function(x) (str_match(x, "^\\d+") < 80 & str_match(x, "\\d+$") < 95)))

# color captain tree by clade
fam01.node <- getMRCA(captainTree.rooted, c("hiscap4_QSS52209.1", "parbra1_6028", "cocimm1_898", "aspasp11_tyr520", "aspasp8_tyr630", "aspawa1_GCB27642.1", "metmaj1_tyr8194", "fuscom1_510842", "fusvir6_tyr6937", "kredeu2_tyr7716", "spobra1_10170_10169", "metbru1_2543", "fusoxy46_KAH7490251.1"))
fam02.node <-  getMRCA(captainTree.rooted, c("asppro1_tyr2437", "paevar1_195182", "aspfis1_EAW15512.1", "asplen2_tyr1741", "aspnom1_158252", "talcel1_GAM38547.1", "aspcal2_761163", "aspcal2_tyr804", "penchr2_71351", "aspche2_tyr880"))
fam03.node <-  getMRCA(captainTree.rooted, c("aspsyd1_145891", "aspmel1_tyr1928", "aspfis1_EAW21814.1", "aspfum3_9656", "asppar10_tyr2350", "aspind2_339748", "aspwes1_tyr3208", "blader2_tyr3388", "cocimm1_4925", "orboli2_TGJ70084.1"))
fam04.node <-  getMRCA(captainTree.rooted, c("altalt9_tyr85", "fuscul4_tyr5444"))
fam05.node <-  getMRCA(captainTree.rooted, c("aspfla7_tyr1211", "aspnig3_284515", "aspluc13_tyr1861", "aspfum4_KMK60551.1", "aspfla13_tyr1195", "aspnig13_tyr2056", "cocpos5_KMM69757.1", "cocimm1_185", "aspter6_tyr2909", "altlon1_tyr173", "aspawa1_GCB26061.1", "aspfum6_EDP49496.1", "aspjen1_tyr1568", "cocimm4_KMP10201.1", "aspche2_tyr886", "aspfum6_EDP50805.1", "colfru1_KAE9577683.1", "fussuc1_tyr6724", "aspfum6_EDP47816.1", "bipsor1_EMD65867.1"))
fam06.node <-  getMRCA(captainTree.rooted, c("aspnig15_tyr2064", "aspsoj2_tyr2670", "aspfum15_tyr1309", "aspory10_tyr2258", "aspfla10_tyr1186", "aspfum4_KMK54996.1", "aspnid3_tyr2021", "aspste3_tyr2739", "aspoch5_tyr2231", "asppar4_tyr2392", "fonnub1_8242", "fusspo3_tyr6697", "colaen1_KAF5522205.1", "fusfuj10_tyr5560", "altgai1_tyr156", "trirub8_EZG02207.1", "aspfla12_tyr1191", "aspfum13_tyr1277", "aspnig3_283767"))
fam07.node <-  getMRCA(captainTree.rooted, c("aspnig4_7123", "aspara1_tyr457", "aspfla12_tyr1193", "aspfum15_tyr1308", "aspche1_tyr850", "aspoch6_tyr2235", "penchr1_145210_145211_145212", "pensp.9_KAF4763348.1", "pensp.9_KAF4765341.1", "colfru1_KAE9574038.1", "colfru6_tyr4318", "fusver1_19403", "fusoxy38_RKK18157.1", "fushos1_tyr5945", "fusoxy23_2378", "ilydes1_KAH7008698.1"))
fam09.node <- getMRCA(captainTree.rooted, c("pyrory1_EHA56134.1", "pyrory2_QBZ54310.1", "pyrory2_QBZ61703.1", "verdah3_443143", "fusequ2_tyr5478", "fusfus52_tyr5838", "fusxyl3_KAG5773805.1", "fusoxy45_KAG7424752.1", "fussol1_232585", "fusoxy43_KAG7003254.1", "kredeu2_504220"))
fam10.node <- getMRCA(captainTree.rooted, c("altalt8_tyr84", "altalt9_tyr86", "altsol2_tyr196", "macpha1_749149", "altalt11_tyr67"))

# annotate all strongly supported branches in bold
captainTree.gg <- ggtree(captainTree.rooted, aes(size = str_match(label, "^\\d+") >= 80 | str_match(label, "\\d+$") >= 95 | isTip == T)) +
  geom_tiplab(size=0) +
  scale_size_manual(values=c(0.75, .25)) +
  layout_rectangular() +
  theme_tree(legend.position=c(.05, -.85)) + 
  geom_hilight(mapping=aes(subset = node %in% fam01.node), fill = "#d7301f", type = "rect") +
  geom_hilight(mapping=aes(subset = node %in% fam02.node), fill = "#f03b20", type = "rect") +
  geom_hilight(mapping=aes(subset = node %in% fam03.node), fill = "#fd8d3c", type = "rect") +
  geom_hilight(mapping=aes(subset = node %in% fam04.node), fill = "#fecc5c", type = "rect") +
  geom_hilight(mapping=aes(subset = node %in% fam05.node), fill = "#fee391", type = "rect") +
  geom_hilight(mapping=aes(subset = node %in% fam06.node), fill = "#8279B8", type = "rect") +
  geom_hilight(mapping=aes(subset = node %in% fam07.node), fill = "#9e9ac8", type = "rect") +
  geom_hilight(mapping=aes(subset = node %in% fam09.node), fill = "#197b41", type = "rect") +
  geom_hilight(mapping=aes(subset = node %in% fam10.node), fill = "#31a354", type = "rect") +
  geom_tippoint(mapping=aes(subset = label %in% captain5s$captainID), shape = 21, size = 2.5, fill = "white", color = "black") +
  scale_x_ggtree() 
msaplot(captainTree.gg, "direct_repeat_cores_captains.fasta", offset=0.1, width=1.5)

#--------------------------------------------------------------#
# Direct repeat alignments onto 5s rDNA
#--------------------------------------------------------------#

# Set wd
setwd("../../Figure5/directRepeats//")

# read in data, and reorder captainIDs according to clade
DRcoords <- read.table("cap26_targetSiteCoords.txt", sep = "\t", header = T)
captainOrder <- c("spobra1_10170_10169", "aspcal2_tyr804", "laspus2_tyr7798", "aspche2_tyr886", "aspnig3_284515", "aspfla10_tyr1186", "aspfum4_KMK54996.1", "aspfum15_tyr1309", "aspnid3_tyr2021", "aspnig15_tyr2064", "aspoch5_tyr2231", "aspory10_tyr2258", "asppar4_tyr2392", "aspsoj2_tyr2670", "aspste3_tyr2739", "aspfla12_tyr1191", "altalt8_tyr84", "macpha1_749149", "annsty1_389805", "fusequ2_tyr5478", "fusfus52_tyr5838", "fusoxy43_KAG7003254.1", "fusoxy45_KAG7424752.1", "fussol1_232585", "fusxyl3_KAG5773805.1", "pyrory1_EHA56134.1", "pyrory2_QBZ54310.1")
DRcoords$captainID <- factor(DRcoords$captainID, levels = rev(captainOrder))

ggplot(data = DRcoords %>%
         remove_missing()) +
  geom_segment(aes(y = captainID, yend = captainID, x = begin_5s - 0.3, xend = end_5s + 0.3), color = "black", linewidth = 5) + 
  geom_segment(aes(y = captainID, yend = captainID, x = begin_5s, xend = end_5s, color = familyID), linewidth = 4) + 
  scale_color_manual(values = c("clade01" = "darkgray", "fam01-5" = "#fee391", "fam02-1" = "#8279B8", "fam03-1" = "#197b41", "fam03-4" = "#addd8e")) +
  ggpubfigs::theme_grey() + 
  labs(x = "5s rDNA coordinates") +
  theme(text  = element_text(size = 8, family = "Helvetica"),
        panel.background = element_rect(fill = '#eeeeee', color = 'white')) + 
  xlim(0,123)
  
  





