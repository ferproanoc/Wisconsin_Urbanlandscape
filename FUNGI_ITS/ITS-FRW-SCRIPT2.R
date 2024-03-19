###################################################################################
############################6: DATA EXPLORATION & FIGURES 1-2######################
###################################################################################

### Clear workspace ###
rm(list=ls())

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")

#Re-read the data to create the Phyloseq object
count_tab <- read.table("ASVs_ITS_countsFORWARD_USE.tsv", sep="\t", header=T, row.names=1) 
tax_tab <- read.table("ASVs_taxa_ITSFRW_USE.tsv", sep="\t", header=T, row.names=1)
sample_info_tab <- read.table("ITS_METADATA_use.txt", sep="\t", header=TRUE, row.names=1, fileEncoding = "UTF-8")

asvs.t = t(count_tab) # t(x) = Transpose 'otu_table-class' or 'phyloseq-class'

# Create a phyloseq object
OTU <- otu_table(asvs.t, taxa_are_rows=F)
SAM <- sample_data(sample_info_tab,errorIfNULL=TRUE)
TAX <- tax_table(as.matrix(tax_tab), errorIfNULL=TRUE)
data_phylo <- phyloseq(OTU, TAX, SAM)
data_phylo 

#Rarefaction curve
library(vegan)
rarecurve(asvs.t, step = 100, cex = 0.75, las = 1) #the total number of reads per sample vary between samples
#rarecurve(count_tab, step = 100, cex = 0.75, las = 1)

#Library size
sum_seq <- rowSums(asvs.t)
plot(sum_seq, ylim=c(0,100000), main=c("Number of counts per sample"), xlab=c("Samples"))
sum_seq

min(sum_seq)
max(sum_seq) # there is differences in the library size for the different samples

#######################################################################################
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

#To make data frames with the previously made phyloseq object
sample.data.ITS = data.frame(sample_data(data_phylo))

alpha.ITS = estimate_richness(data_phylo)

######################################
###########A. SHANNON DIVERSITY#######
######################################

##Alpha diversity (Shannon diversity index) across different sections 
#(e.g., Leaf, Thatch, Rhizosphere, Bulk) in your microbiome dataset, considering the factor "Management" as a predictor. 

sample.data.ITS$Shannon = alpha.ITS$Shannon

#handling missing values using na.omit()
sample.data.ITS.Shannon <- na.omit(sample.data.ITS)

sample.data.ITS$Section = factor(sample.data.ITS$Section, levels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere"))
sample.data.ITS$Management = factor(sample.data.ITS$Management, levels = c("Greens", "Home Lawn", "Unmanaged", "Fairway"))

#install.packages("multcompView")
library(multcompView)

sample.data.ITS.Shannon = sample.data.ITS %>% 
  nest(data = -Section) %>% 
  mutate(Section = factor(Section, levels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere"),
                                            labels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere")),
                    glm = map(data, ~ glm(Shannon ~ Management, data = ., family = quasipoisson)),
                    anova = map(data, ~ Anova(glm(Shannon ~ Management, data = . , family = quasipoisson))),
                    emmeans = map(data, ~ emmeans(glm(Shannon ~ Management, data = ., family = quasipoisson),  ~  Management, type='response')),
                    stat.grp = map(data, ~ multcomp::cld(emmeans(glm(Shannon ~ Management, data = . , family = quasipoisson),  ~  Management, type='response'), Letter="ABCDEFGHIJKLMNOPQRSTUVWZYZ")), 
                    disease.CI = map(data, ~as.data.frame(summary(emmeans(glm(Shannon ~ Management, data = . , family = quasipoisson),  ~  Management, type='response'))))) %>%
  unnest(c(stat.grp))

sample.data.ITS.Shannon$.group = gsub(" ", "",sample.data.ITS.Shannon$.group)

######################################
###### Plotting Shannon
######################################

#install.packages("cowplot")
library(cowplot)
packageVersion("cowplot")
citation("cowplot")

FShannon <- ggplot() +
  geom_errorbar(data = sample.data.ITS.Shannon, aes(x=Section, ymin = asymp.LCL, ymax = asymp.UCL, color=Management), 
                position = position_dodge(width = 0.8), width = 0, alpha=0.7, linewidth=4, show.legend = FALSE) +
  geom_point(data = sample.data.ITS.Shannon,aes(x = Section, y = rate, color = Management), 
             position = position_dodge(width = 0.8), size = 6, show.legend = FALSE, alpha = 0.99) +
  geom_text(data = sample.data.ITS.Shannon, aes(label = .group, x = Section, group = Management, y = asymp.UCL+(0.004*max(asymp.UCL))), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 6, hjust = 0.5) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 'twodash', alpha = 0.4)+
  xlab('') + ylab("Shannon Diversity Index") +
  ylim(NA, max(sample.data.ITS.Shannon$asymp.UCL+(0.04*max(sample.data.ITS.Shannon$asymp.UCL)))) +
  theme_light(base_size = 20, base_family = "Times") + 
  theme(panel.background = element_rect(fill = '#FAFAFA')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("Greens" = "#d95f0e", "Home Lawn" = "#8856a7", "Unmanaged" = "#c51b8a", "Fairway" = "#117864"),
                     name = "Management")

FShannon
# Increase the plot margins
FShannon <- FShannon +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_FShannon <- ggdraw() +
  draw_plot(FShannon) +  # Add the plot
  draw_label("Fungi", fontface = "bold", size = 18, x = 0.5, y = 0.98) +  # Add the title
  draw_label("B", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_FShannon

ggsave("Shannon_ITSF.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

######################################
########B. OBSERVED DIVERSITY#########
######################################

#To make data frames with the previously made phyloseq object
sample.data.ITS$Observed = alpha.ITS$Observed

#handling missing values using na.omit() #Run after line 128
sample.data.ITS.Observed <- na.omit(sample.data.ITS)

sample.data.ITS.Observed = sample.data.ITS %>% 
  nest(data = -Section) %>% mutate(Section = factor(Section, levels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere"),
                                             labels = c("Leaf", "Thatch", "Rhizoplane", "Rhizosphere")),
                            glm = map(data, ~ glm(Observed ~ Management, data = . , family = quasipoisson)),
                            anova = map(data, ~ Anova(glm(Observed ~ Management, data = . , family = quasipoisson))),
                            emmeans = map(data, ~ emmeans(glm(Observed ~ Management, data = . , family = quasipoisson),  ~  Management, type='response')),
                            stat.grp = map(data, ~ multcomp::cld(emmeans(glm(Observed ~ Management, data = . , family = quasipoisson),  ~  Management, type='response'), Letter="ABCDEFGHIJKLMNOPQRSTUVWZYZ")), 
                            disease.CI = map(data, ~as.data.frame(summary(emmeans(glm(Observed ~ Management, data = . , family = quasipoisson),  ~ Management, type='response'))))) %>%
  unnest(c(stat.grp))

sample.data.ITS.Observed$.group = gsub(" ", "",sample.data.ITS.Observed$.group)

######################################
# Plotting Observed
######################################

FObserved <-  ggplot() +
  geom_errorbar(data = sample.data.ITS.Observed, aes(x=Section, ymin = asymp.LCL, ymax = asymp.UCL, color=Management), 
                position = position_dodge(width = 0.8), width = 0, alpha=0.7, linewidth=4, show.legend = F) +
  geom_point(data = sample.data.ITS.Observed, aes(x = Section, y = rate, color = Management), 
             position = position_dodge(width = 0.8), size = 6, show.legend = FALSE, alpha = 0.99) +
  geom_text(data = sample.data.ITS.Observed, aes(label = .group, x = Section, group = Management, y = asymp.UCL+(0.004*max(asymp.UCL))), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 6, hjust=0.5) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), lty = 'twodash', alpha = 0.4)+
  xlab('') + ylab("Observed Richness") +
  ylim(NA, max(sample.data.ITS.Observed$asymp.UCL+
                 (0.08*max(sample.data.ITS.Observed$asymp.UCL))))+
  theme_light(base_size = 20, base_family = "Times") + 
  theme(panel.background = element_rect(fill = '#FAFAFA')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("Greens" = "#d95f0e", "Home Lawn" = "#8856a7", "Unmanaged" = "#c51b8a", "Fairway" = "#117864"),
                     name = "Management")
FObserved

# Increase the plot margins
FObserved <- FObserved +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_FObserved <- ggdraw() +
  draw_plot(FObserved) +  # Add the plot
  draw_label("Fungi", fontface = "bold", size = 18, x = 0.5, y = 0.98) +  # Add the title
  draw_label("A", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_FObserved

ggsave("Observed_ITSF.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

#######################################################################################
###############################Figure 2: Beta diversity################################
########################################################################################

# Rarefy to an even depth
set.seed(111)  # setting seed for reproducibility
data_phylo.rarefied = rarefy_even_depth(data_phylo)

# Normalize read counts (this gives relative abundance)
data_phylo.ITS.norm = transform_sample_counts(data_phylo.rarefied, function(x) x / sum(x))
data_phylo.ITS.norm.nmds = ordinate(data_phylo.ITS.norm , method = 'NMDS', distance = 'bray', k = 3, maxtry = 1000)

######################################
################Plot
######################################
                                              
library(RColorBrewer)

plot_ordination(data_phylo.ITS.norm, data_phylo.ITS.norm.nmds, color = "Management", shape = "Section") + 
  geom_point(aes(color = Management, shape = Section), size = 7) +
  scale_color_brewer(palette = "Dark2") +  # Use colorblind-friendly palette from RColorBrewer
  theme_bw(base_size = 20, base_family = "Times")

data_otu_rarefied = data.frame(otu_table(data_phylo.rarefied))

# Access stress value from NMDS ordination object
stress_value <- data_phylo.ITS.norm.nmds$stress
print(paste("Stress value:", stress_value))

# Load required package
library(vegan)

# Permanova test using the vegan package
permanova_result <- adonis2(data_otu_rarefied~Management, data=sample_info_tab, permutations=999, method="bray")
permanova_result

# Permanova test using the vegan package
permanova_result2 <- adonis2(data_otu_rarefied~Management+Section, data=sample_info_tab, permutations=999, method="bray")
permanova_result2


# Plot the ordination with ggplot
betaITS <- plot_ordination(data_phylo.ITS.norm , data_phylo.ITS.norm.nmds, color = "Management", shape = 'Section') + 
  geom_point(size = 7) + 
  scale_color_brewer(palette = "Dark2")  +
  theme_bw(base_size = 20, base_family = "Times") +
  # Add the R2, p-value, and stress annotations within the plot
  annotate("text", x = -1.4, y = -1.2, 
           label = paste("R2 =", 0.353), 
           vjust = 0, hjust = 0, size = 3, color = "black") +
  annotate("text", x = -1.4, y = -1.3, 
           label = paste("p-value =", 0.001), 
           vjust = 0, hjust = 0, size = 3, color = "black") +
  annotate("text", x = -1.4, y = -1.4, 
           label = paste("Stress =", 0.0846), 
           vjust = 0, hjust = 0, size = 3, color = "black")
betaITS

# Increase the plot margins
betaITS <- betaITS +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_betaITS <- ggdraw() +
  draw_plot(betaITS) +  # Add the plot
  draw_label("Fungi", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("A", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_betaITS

ggsave("BetaDiv_ITSF.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

#END 
#################################################################################################################################
#################################################################################################################################
