###################################################################################
######################7:COMMUNITY COMPOSITION - Figure 3 ##########################
###################################################################################

### Clear workspace ###
rm(list=ls())

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
data_phylo 

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
packageVersion("dplyr")

# Group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(dat, ~Phylum, function(x) c(median = median(x$Abundance)))

# Find Phyla whose rel. abund. is less than 1%
remainder <- medians[medians$median <= 0.02, ]$Phylum
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
cbbPalette <- c( "#E69F00", "#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#9933CC", "#999999", "sienna","tomato","peru","blue3","coral2","firebrick")

#data analysis 
#stacked bar plot 
colnames(dat_summary)[colnames(dat_summary) == "Abundance"] <- "mean_abundance"

# Manually specify the order of Phylum levels
phylum_order <- c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Chloroflexi", "Myxococcota", "Planctomycetota", "Proteobacteria", "Verrucomicrobiota", "Others")

# Convert Phylum to a factor with the custom order
dat_summary$Phylum <- factor(dat_summary$Phylum, levels = phylum_order)

# Stacked bar plot
b16Scom <- ggplot(dat_summary, aes(x = Section, y = mean_abundance, fill = Phylum)) + 
  geom_col() +
  facet_grid(. ~ Management, scales = "free") + 
  ylab("Relative Abundance") + xlab(NULL)+
  theme_gray(base_size = 20, base_family = "Times") + 
  scale_fill_manual(values = cbbPalette, na.value = "grey") + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
b16Scom 

# Increase the plot margins
b16Scom  <- b16Scom  +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_b16Scom  <- ggdraw() +
  draw_plot(b16Scom ) +  # Add the plot
  draw_label("Bacteria", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("B", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_b16Scom 

ggsave("16SF_Phyla_community.tiff", width = 12, height = 8.2, units = "in", dpi = 600)

#END
#############################################################################################
#############################################################################################


