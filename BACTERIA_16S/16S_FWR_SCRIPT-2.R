###################################################################################
############################6: DATA EXPLORATION & FIGURES 1-2######################
###################################################################################

### Clear workspace ###
rm(list=ls())

library("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")

#Re-read the data to create the Phyloseq object
count_tab <- read.table("16S_ASVs_countsF_USE.tsv", sep="\t", header=T, row.names=1) 
tax_tab <- read.table("16S_ASVs_Taxonomy_F_USE.tsv", sep="\t", header=T, row.names=1)
sample_info_tab <- read.table("16S_METADATA_use.txt", sep="\t", header=TRUE, row.names=1, fileEncoding = "UTF-8")
asvs.t = t(count_tab) # t(x) = Transpose 'otu_table-class' or 'phyloseq-class'

# Create a phyloseq object
OTU <- otu_table(asvs.t, taxa_are_rows=F)
SAM <- sample_data(sample_info_tab,errorIfNULL=TRUE)
TAX <- tax_table(as.matrix(tax_tab), errorIfNULL=TRUE)
data_phylo <- phyloseq(OTU, TAX, SAM)

#Rarefaction curve
library(vegan)
rarecurve(asvs.t, step = 100, cex = 0.75, las = 1) #the total number of reads per sample vary between samples

#Library size
sum_seq <- rowSums(asvs.t)
plot(sum_seq, ylim=c(0,100000), main=c("Number of counts per sample"), xlab=c("Samples"))

min(sum_seq)
max(sum_seq) # there is differences in the library size for samples

########################################################################################
###############################Figure 1: Alpha diversity################################
########################################################################################

# Load necessary libraries
library(ggsci)# collection of 'ggplot2' color palettes inspired by plots in scientific journals
library(Rmisc)#Contains many functions useful for data analysis and utility operations
library(tidyverse)
library(lme4)
library(multcomp)
library(emmeans)
library(car)
library(janitor)
library(broom)
library (stats)

#To make data frames with the previously made phyloseq object
sample.data.16S = data.frame(sample_data(data_phylo))

alpha.16S = estimate_richness(data_phylo)

##########################################
###########A. SHANNON DIVERSITY###########
##########################################

##Alpha diversity (Shannon diversity index) across different sections 
#(e.g., Leaf, Thatch, Rhizosphere, Bulk) in your microbiome dataset, considering the factor "Management" as a predictor. 

sample.data.16S$Shannon = alpha.16S$Shannon

#handling missing values using na.omit()
sample.data.16S.Shannon <- na.omit(sample.data.16S)

sample.data.16S$Section = factor(sample.data.16S$Section, levels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere"))
sample.data.16S$Management = factor(sample.data.16S$Management, levels = c("Greens", "Home Lawn", "Unmanaged", "Fairway"))

#install.packages("multcompView")
library("multcompView")

sample.data.16S.Shannon = sample.data.16S %>% 
  nest(data = -Section) %>% 
  mutate(Section = factor(Section, levels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere")),
         glm = map(data, ~ glm(Shannon ~ Management, data = ., family = quasipoisson)),
         anova_result = map(glm, ~ Anova(.x, type = "II", test = "F") %>% broom::tidy()),  # Store ANOVA results
         anova_summary = map(anova_result, summary),  # Optionally, summarize ANOVA results
         emmeans = map(data, ~ emmeans(glm(Shannon ~ Management, data = ., family = quasipoisson),  ~  Management, type='response')),
         stat.grp = map(data, ~ multcomp::cld(emmeans(glm(Shannon ~ Management, data = . , family = quasipoisson),  ~  Management, type='response'), Letter="ABCDEFGHIJKLMNOPQRSTUVWZYZ")), 
         disease.CI = map(data, ~as.data.frame(summary(emmeans(glm(Shannon ~ Management, data = . , family = quasipoisson),  ~  Management, type='response'))))) %>%
  unnest(c(stat.grp))

# Print the first few rows to verify if the 'anova' component contains the ANOVA results
print(head(sample.data.16S.Shannon))

sample.data.16S.Shannon$.group = gsub(" ", "",sample.data.16S.Shannon$.group)

##########################################
##############Plotting Shannon############
##########################################

BShannon <- ggplot() +
  geom_errorbar(data = sample.data.16S.Shannon, aes(x=Section, ymin = asymp.LCL, ymax = asymp.UCL, color=Management), 
                position = position_dodge(width = 0.8), width = 0, alpha=0.7, linewidth=4, show.legend = F) +
  geom_point(data = sample.data.16S.Shannon,aes(x = Section, y = rate, color = Management), 
             position = position_dodge(width = 0.8), size = 6, show.legend = FALSE, alpha = 0.99) +
  geom_text(data = sample.data.16S.Shannon, 
            aes(label = .group, x = Section, group = Management, y = asymp.UCL+(0.004*max(asymp.UCL))), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 6, hjust = 0.5) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 'twodash', alpha = 0.4)+
  xlab('') + ylab("Shannon Diversity Index") +
  ylim(NA, max(sample.data.16S.Shannon$asymp.UCL+
                 (0.04*max(sample.data.16S.Shannon$asymp.UCL))))+
  theme_light(base_size = 20, base_family = "Times") + 
  theme(panel.background = element_rect(fill = '#FAFAFA')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("Greens" = "#d95f0e", "Home Lawn" = "#8856a7", "Unmanaged" = "#c51b8a", "Fairway" = "#117864"),
                     name = "Management")

# Increase the plot margins
BShannon <- BShannon +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_BShannon <- ggdraw() +
  draw_plot(BShannon) +  # Add the plot
  draw_label("Bacteria", fontface = "bold", size = 18, x = 0.5, y = 0.98) +  # Add the title
  draw_label("D", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

ggsave("Shannon_16SF.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

# Print the structure of the ANOVA results
print(sample.data.16S.Shannon$anova_result[[1]])

# Extract p-values from ANOVA results
p_values_anova <- sample.data.16S.Shannon %>%
  mutate(p_value_anova = map_dbl(anova_result, ~ filter(.x, term == "Management")$p.value))

##########################################
############B. OBSERVED DIVERSITY#########
##########################################

#To make data frames with the previously made phyloseq object
sample.data.16S$Observed = alpha.16S$Observed

#handling missing values using na.omit() #Run after line 128
sample.data.16S.Observed <- na.omit(sample.data.16S)

sample.data.16S.Observed = sample.data.16S %>% 
  nest(data = -Section) %>% mutate(Section = factor(Section, levels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere"),
                                                    labels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere")),
                                   glm = map(data, ~ glm(Observed ~ Management, data = . , family = quasipoisson)),
                                   anova = map(data, ~ Anova(glm(Observed ~ Management, data = . , family = quasipoisson))),
                                   emmeans = map(data, ~ emmeans(glm(Observed ~ Management, data = . , family = quasipoisson),  ~  Management, type='response')),
                                   stat.grp = map(data, ~ multcomp::cld(emmeans(glm(Observed ~ Management, data = . , family = quasipoisson),  ~  Management, type='response'), Letter="ABCDEFGHIJKLMNOPQRSTUVWZYZ")), 
                                   disease.CI = map(data, ~as.data.frame(summary(emmeans(glm(Observed ~ Management, data = . , family = quasipoisson),  ~ Management, type='response'))))) %>%
  unnest(c(stat.grp))

sample.data.16S.Observed$.group = gsub(" ", "",sample.data.16S.Observed$.group)

##########################################
############## Plotting Observed##########
##########################################

BObserved <- ggplot() +
  geom_errorbar(data = sample.data.16S.Observed, aes(x=Section, ymin = asymp.LCL, ymax = asymp.UCL, color=Management), 
                position = position_dodge(width = 0.8), width = 0, alpha=0.7, linewidth=4, show.legend = F) +
  geom_point(data = sample.data.16S.Observed, aes(x = Section, y = rate, color = Management), 
             position = position_dodge(width = 0.8), size = 6, show.legend = FALSE, alpha = 0.99) +
  geom_text(data = sample.data.16S.Observed, aes(label = .group, x = Section, group = Management, y = asymp.UCL+(0.004*max(asymp.UCL))), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 6, hjust = 0.5) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 'twodash', alpha = 0.4)+
  xlab('') + ylab("Observed Richness") +
  ylim(NA, max(sample.data.16S.Observed$asymp.UCL+
                 (0.08*max(sample.data.16S.Observed$asymp.UCL))))+
  theme_light(base_size = 20, base_family = "Times") + 
  theme(panel.background = element_rect(fill = '#FAFAFA')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("Greens" = "#d95f0e", "Home Lawn" = "#8856a7", "Unmanaged" = "#c51b8a", "Fairway" = "#117864"),
                     name = "Management")

# Increase the plot margins
BObserved <- BObserved +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_BObserved <- ggdraw() +
  draw_plot(BObserved) +  # Add the plot
  draw_label("Bacteria", fontface = "bold", size = 18, x = 0.5, y = 0.98) +  # Add the title
  draw_label("C", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

ggsave("Observed_16SF.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

#######################################################################################
###############################Figure 2: Beta diversity################################
########################################################################################

# Rarefy to an even depth
set.seed(111)  # setting seed for reproducibility
data_phylo.rarefied = rarefy_even_depth(data_phylo)

# Normalize read counts (this gives relative abundance)
data_phylo.16S.norm = transform_sample_counts(data_phylo.rarefied, function(x) x / sum(x))
data_phylo.16S.norm.nmds = ordinate(data_phylo.16S.norm , method = 'NMDS', distance = 'bray', k = 3, maxtry = 1000)

data_phylo.16S.norm
                    
##########################################
#################Plot#####################
##########################################
                                              
library(RColorBrewer)

beta16S <- plot_ordination(data_phylo.16S.norm, data_phylo.16S.norm.nmds, color = "Management", shape = "Section") + 
  geom_point(aes(color = Management, shape = Section), size = 7) +
  scale_color_brewer(palette = "Dark2") +  # Use colorblind-friendly palette from RColorBrewer
  theme_bw(base_size = 20, base_family = "Times")
beta16S

# Access stress value from NMDS ordination object
stress_value <- data_phylo.16S.norm.nmds$stress
print(paste("Stress value:", stress_value))

data_otu_rarefied = data.frame(otu_table(data_phylo.rarefied))

# Load required package
library(vegan)

# Permanova test using the vegan package
permanova_result <- adonis2(data_otu_rarefied~Management, data=sample_info_tab, permutations=999, method="bray")
permanova_result
permanova_result2 <- adonis2(data_otu_rarefied~Management + Section, data=sample_info_tab, permutations=999, method="bray")
permanova_result2

# Plot the ordination with ggplot
beta16S <- plot_ordination(data_phylo.16S.norm , data_phylo.16S.norm.nmds, color = "Management", shape = 'Section') + 
  geom_point(size = 7) + 
  scale_color_brewer(palette = "Dark2")  +
  theme_bw(base_size = 20, base_family = "Times") +
  # Add the R2, p-value, and stress annotations within the plot
  annotate("text", x = -1.4, y = -1.2, 
           label = paste("R2 =", 0.289), 
           vjust = 0, hjust = 0, size = 3, color = "black") +
  annotate("text", x = -1.4, y = -1.3, 
           label = paste("p-value =", 0.001), 
           vjust = 0, hjust = 0, size = 3, color = "black") +
  annotate("text", x = -1.4, y = -1.4, 
           label = paste("Stress =", 0.0575), 
           vjust = 0, hjust = 0, size = 3, color = "black")
beta16S

# Increase the plot margins
beta16S <- beta16S +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_beta16S <- ggdraw() +
  draw_plot(beta16S) +  # Add the plot
  draw_label("Bacteria", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("B", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_beta16S

ggsave("BetaDiv_16SF_rare.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

#END 
#################################################################################################################################
#################################################################################################################################
