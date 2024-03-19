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

##End
#############################################################################################
#############################################################################################
