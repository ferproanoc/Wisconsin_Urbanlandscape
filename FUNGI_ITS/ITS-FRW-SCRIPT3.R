###############################################################################################################################################
#########################################################ITS_DATASET R ANALYSIS-3###################################################################
##############################################################################################################################################
######Fer Proano Cuenca #################################################################FEB 2024#############################################

###################################################################################
######################7:COMMUNITY COMPOSITION - Figure 3 ##########################
###################################################################################

### Clear workspace ###
rm(list=ls())

#Re-read the data to create the Phyloseq object
count_tab <- read.table("ASVs_ITS_countsFORWARD_USE.tsv", sep="\t", header=T, row.names=1) 
tax_tab <- read.table("ASVs_taxa_ITSFRW_USE.tsv", sep="\t", header=T, row.names=1)
sample_info_tab <- read.table("ITS_METADATA_use.txt", sep="\t", header=TRUE, row.names=1, fileEncoding = "UTF-8")

asvs.t = t(count_tab) # t(x) = Transpose 'otu_table-class' or 'phyloseq-class'

# Print unique phyla before filtering
print(unique(tax_tab$Phylum))

# Exclude specified phyla from the taxonomy table
phyla_to_remove <- c("Cercozoa", "Metazoa_phy_Incertae_sedis") #"Rozellomycota"
tax_tab_filtered <- tax_tab[!tax_tab$Phylum %in% phyla_to_remove, ]

# Create a phyloseq object
OTU <- otu_table(asvs.t, taxa_are_rows=F)
SAM <- sample_data(sample_info_tab,errorIfNULL=TRUE)
#TAX <- tax_table(as.matrix(tax_tab), errorIfNULL=TRUE)
TAX <- tax_table(as.matrix(tax_tab_filtered), errorIfNULL=TRUE)
data_phylo <- phyloseq(OTU, TAX, SAM)

set.seed(111) # keep result reproductive
data_phylo.rarefied = rarefy_even_depth(data_phylo)
#data_phylo.rarefied = prune_samples(sample_sums(data_phylo.rarefied) > 5000,data_phylo.rarefied)

#normalize phyloseq object 
ps.norm = transform_sample_counts(data_phylo.rarefied, function(x) x / sum(x))

#Creating remainder category 
# agglomerate taxa
glom <- tax_glom(ps.norm, taxrank = 'Phylum')

# create dataframe from phyloseq object
dat <- psmelt(glom)

# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)

library(ggplot2)
library(dplyr)
library(plyr)

# Group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(dat, ~Phylum, function(x) c(median = median(x$Abundance)))

# Find Phyla whose rel. abund. is less than 1%
remainder <- medians[medians$median <= 0.001, ]$Phylum
remainder 
# Change their name to "Others"
dat$Phylum[dat$Phylum %in% remainder] <- 'Others'

# Group and summarize by Section, Management, Phylum, and Sample
dat_grouped <- aggregate(Abundance ~ Section + Management + Phylum + Sample, data = dat, FUN = sum)

unique_combinations <- unique(dat_grouped[, c("Section", "Management", "Phylum")])
print(unique_combinations)

# Summarize the data using aggregate()
dat_summary <- aggregate(Abundance ~ Section + Management + Phylum, data = dat_grouped, FUN = mean)

# Print the resulting dataframe
print(head(dat_summary))

# Color palette for stacked bar graph:
cbbPalette <- c("#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#999999", "#D55E00", "#9933CC")

#stacked bar plot 
colnames(dat_summary)[colnames(dat_summary) == "Abundance"] <- "mean_abundance"

# Manually specify the order of Phylum levels
phylum_order <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Glomeromycota", "Mortierellomycota", "Rozellomycota", "Others")

# Convert Phylum to a factor with the custom order
dat_summary$Phylum <- factor(dat_summary$Phylum, levels = phylum_order)

# Stacked bar plot
ITScom <- ggplot(dat_summary, aes(x = Section, y = mean_abundance, fill = Phylum)) + 
  geom_col() +
  facet_grid(. ~ Management, scales = "free") + 
  ylab("Relative Abundance") + xlab(NULL)+
  theme_grey(base_size = 20, base_family = "Times") + 
  scale_fill_manual(values = cbbPalette) + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
ITScom 

# Increase the plot margins
ITScom <- ITScom +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

library(cowplot)

# Create a ggdraw object
plot_ITScom <- ggdraw() +
  draw_plot(ITScom) +  # Add the plot
  draw_label("Fungi", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("A", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_ITScom

ggsave("ITSF_Phyla_community.tiff", width = 12, height = 8.2, units = "in", dpi = 600)

#############################End#############################################################################################


#####################################################################Not run#######################################################
#################################################Another way tp create a stacked plot###############################################
library("phyloseq")
library("ggplot2")
library("tidyverse")
library("BiocManager")
BiocManager::install("microbiome")
library("devtools") 
install_github("microbiome/microbiome")
library("microbiome")
library("scales")

#Re-read the data to create the Phyloseq object
count_tab <- read.table("ASVs_ITS_countsFORWARD_USE.tsv", sep="\t", header=T, row.names=1) 
tax_tab <- read.table("ASVs_taxa_ITSFRW_USE.tsv", sep="\t", header=T, row.names=1)
sample_info_tab <- read.table("ITS_METADATA.txt", sep="\t", header=TRUE, row.names=1, fileEncoding = "UTF-8")

asvs.t = t(count_tab) # t(x) = Transpose 'otu_table-class' or 'phyloseq-class'

# Create a phyloseq object
OTU <- otu_table(asvs.t, taxa_are_rows=F)
SAM <- sample_data(sample_info_tab,errorIfNULL=TRUE)
TAX <- tax_table(as.matrix(tax_tab), errorIfNULL=TRUE)
data_phylo <- phyloseq(OTU, TAX, SAM)
data_phylo 

# Calculate the Phylum quantity 
set.seed(111) # keep result reproductive
data_phylo.rarefied = rarefy_even_depth(data_phylo)

Fung.phylum <- data_phylo.rarefied %>%
  tax_glom(taxrank = "Phylum") %>%  # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%  # Melt to long format
  filter(Abundance > 0.01) %>% # Filter out low abundance taxa
  arrange(Phylum) # Sort data frame alphabetically by phylum

# Set colors for plotting
Class_phylum <- c("turquoise","sienna","tomato","peru","blue3","coral2","firebrick","green3","yellow3","seagreen4","orange","red","deeppink","cyan","darkorange3","darkviolet","red3","red4","grey81","seagreen1","darkorange4","yellow","yellow4","green","green4","hotpink","blue4","purple","purple3","purple4","tan","tan3","maroon","tan4","black","grey50","grey91","pink","navy","pink3","lawngreen","pink4","lightskyblue")

#Organize x axis
Fung.phylum$Section_new <- factor(Fung.phylum$Section,levels = c("Thatch", "Bulk", "Leaf", "Rhizosphere"))
Fung.phylum$Management_new <- factor(Fung.phylum$Management,levels = c("Greens", "Fairway", "Home Lawn", "Unmanaged"))

#Plot #phylum labels based on abundance

ggplot(Fung.phylum, aes(x = Section_new, y = Abundance, fill = reorder(Phylum, Abundance))) +
  facet_wrap(~Management_new, nrow = 1) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = Class_phylum) +
  theme_bw() +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = NULL, y = "Relative Abundance") +  # Change y-axis label
  theme(axis.title.x = element_blank()) +  # Remove x-axis title completely
  theme(legend.title = element_blank())  # Remove legend title

ggsave("ITSPhyla_community_rare.tiff", dpi = 900, width = 12, height = 7.5, units = 'in')

unique(tax_tab[, "Phylum"])
table(tax_tab[, "Phylum"])
#####################################################################################################################################

