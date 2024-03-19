#############################################################################################################################################
#########################################################ITS_DATASET R ANALYSIS-5###################################################################
##############################################################################################################################################
#######Fer Proano Cuenca ################################################################FEB 2024#############################################

###############################################################################
######################9: R2 BY MANAGEMENT INTENSITY -Figure 6##################
###############################################################################

##################ABOVE GROUND HABITAT: LEAVES AND THATCH######################

### Clear workspace ###
rm(list=ls())

# Packages needed for analysis
library(dplyr)
library(tibble)
library(phyloseq)
library(ape)
library(vegan)
library(FSA)
library(eulerr)

# Packages needed for plotting
library(ggplot2)
library(grid)
library(gridExtra)

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

# Section is the column in your sample_info_tab that denotes the habitat
sample_info_tab$Habitat <- factor(ifelse(sample_info_tab$Section %in% c("Thatch", "Leaf"), "Above_ground", "Below_ground"))

# Update the phyloseq object with the new metadata
SAM <- sample_data(sample_info_tab, errorIfNULL=TRUE)
data_phylo <- phyloseq(OTU, TAX, SAM)
data_phylo
# Assuming "Management" is the column in your sample_info_tab that denotes the management type
sample_info_tab$Management <- factor(sample_info_tab$Management)

# Update the phyloseq object with the new metadata
SAM <- sample_data(sample_info_tab, errorIfNULL=TRUE)
data_phylo <- phyloseq(OTU, TAX, SAM)
data_phylo 

# Create phyloseq object for above ground habitat
data_phylo_above<- subset_samples(data_phylo, Habitat == "Above_ground")

# Create phyloseq objects for specific management types
management_types <- c("Unmanaged", "Home Lawn", "Greens", "Fairway")

# Subset for each management type
data_phylo_unmanaged <- subset_samples(data_phylo_above, Management == "Unmanaged")
data_phylo_home_lawn <- subset_samples(data_phylo_above, Management == "Home Lawn")
data_phylo_greens <- subset_samples(data_phylo_above, Management == "Greens")
data_phylo_fairway <- subset_samples(data_phylo_above, Management == "Fairway")

###########ABOVE GROUND - UNMANAGED:##########
# Rarefy to an even depth
set.seed(111)
data_phylo_unmanaged.rare <- rarefy_even_depth(data_phylo_unmanaged)

# Normalize read counts
data_phylo_unmanaged.norm <- transform_sample_counts(data_phylo_unmanaged.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_unmanaged.rare)

# Calculate the number of individuals in the meta community (Average read depth)
N <- mean(apply(OTU.table, 1, sum))

# Calculate the average relative abundance of each taxa across communities
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p.df = data.frame(p) %>%
  rownames_to_column(var="OTU")

# Calculate the occurrence frequency of each taxa
OTU.table.bi <- 1*(OTU.table>0)
freq.table <- apply(OTU.table.bi, 2, mean)
freq.table <- freq.table[freq.table != 0]
freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

#Combine
C <- inner_join(p.df,freq.df, by="OTU") %>%
  arrange(p)
# Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
C.no0 <- C %>%
  filter(freq != 0, p != 0)

#Calculate the limit of detection
d <- 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
p.list <- C.no0$p
freq.list <- C.no0$freq
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.01))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats.unmanaged <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                                 Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                                 Samples=nrow(OTU.table), Richness=length(p.list), 
                                 Detect=d)

##############ABOVE - HOME LAWN:###############
# Rarefy to an even depth
set.seed(111)
data_phylo_home_lawn.rare <- rarefy_even_depth(data_phylo_home_lawn)

# Normalize read counts
data_phylo_home_lawn.norm <- transform_sample_counts(data_phylo_home_lawn.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_home_lawn.rare)

# Calculate the number of individuals in the meta community (Average read depth)
N <- mean(apply(OTU.table, 1, sum))

# Calculate the average relative abundance of each taxa across communities
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p.df = data.frame(p) %>%
  rownames_to_column(var="OTU")

# Calculate the occurrence frequency of each taxa
OTU.table.bi <- 1*(OTU.table>0)
freq.table <- apply(OTU.table.bi, 2, mean)
freq.table <- freq.table[freq.table != 0]
freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

#Combine
C <- inner_join(p.df,freq.df, by="OTU") %>%
  arrange(p)
# Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
C.no0 <- C %>%
  filter(freq != 0, p != 0)

#Calculate the limit of detection
d <- 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
p.list <- C.no0$p
freq.list <- C.no0$freq
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.01))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats.home_lawn <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                                 Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                                 Samples=nrow(OTU.table), Richness=length(p.list), 
                                 Detect=d)

#################ABOVE GROUND - GREENS:################
# Rarefy to an even depth
set.seed(111)
data_phylo_greens.rare <- rarefy_even_depth(data_phylo_greens)

# Normalize read counts
data_phylo_greens.norm <- transform_sample_counts(data_phylo_greens.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_greens.rare)

# Calculate the number of individuals in the meta community (Average read depth)
N <- mean(apply(OTU.table, 1, sum))

# Calculate the average relative abundance of each taxa across communities
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p.df = data.frame(p) %>%
  rownames_to_column(var="OTU")

# Calculate the occurrence frequency of each taxa
OTU.table.bi <- 1*(OTU.table>0)
freq.table <- apply(OTU.table.bi, 2, mean)
freq.table <- freq.table[freq.table != 0]
freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

#Combine
C <- inner_join(p.df,freq.df, by="OTU") %>%
  arrange(p)
# Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
C.no0 <- C %>%
  filter(freq != 0, p != 0)

#Calculate the limit of detection
d <- 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
p.list <- C.no0$p
freq.list <- C.no0$freq
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.01))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats.greens <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                              Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                              Samples=nrow(OTU.table), Richness=length(p.list), 
                              Detect=d)

#################ABOVE GROUND - FAIRWAY:################
# Rarefy to an even depth
set.seed(111)
data_phylo_fairway.rare <- rarefy_even_depth(data_phylo_fairway)

# Normalize read counts
data_phylo_fairway.norm <- transform_sample_counts(data_phylo_fairway.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_fairway.rare)

# Calculate the number of individuals in the meta community (Average read depth)
N <- mean(apply(OTU.table, 1, sum))

# Calculate the average relative abundance of each taxa across communities
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p.df = data.frame(p) %>%
  rownames_to_column(var="OTU")

# Calculate the occurrence frequency of each taxa
OTU.table.bi <- 1*(OTU.table>0)
freq.table <- apply(OTU.table.bi, 2, mean)
freq.table <- freq.table[freq.table != 0]
freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

#Combine
C <- inner_join(p.df,freq.df, by="OTU") %>%
  arrange(p)
# Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
C.no0 <- C %>%
  filter(freq != 0, p != 0)

#Calculate the limit of detection
d <- 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
p.list <- C.no0$p
freq.list <- C.no0$freq
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.01))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats.fairway <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                               Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                               Samples=nrow(OTU.table), Richness=length(p.list), 
                               Detect=d)

################################BELOW GROUND HABITAT############################

# Create phyloseq object for terrestrial habitat
data_phylo_below <- subset_samples(data_phylo, Habitat == "Below_ground")

# Create phyloseq objects for specific management types
management_types <- c("Unmanaged", "Home Lawn", "Greens", "Fairway")

# Subset for each management type #.b = below ground
data_phylo_unmanaged.b <- subset_samples(data_phylo_below, Management == "Unmanaged")
data_phylo_home_lawn.b <- subset_samples(data_phylo_below, Management == "Home Lawn")
data_phylo_greens.b <- subset_samples(data_phylo_below, Management == "Greens")
data_phylo_fairway.b <- subset_samples(data_phylo_below, Management == "Fairway")

######################BELOW GROUND - UNMANAGED:################

# Rarefy to an even depth
set.seed(111)
data_phylo_unmanaged.b.rare <- rarefy_even_depth(data_phylo_unmanaged.b)

# Normalize read counts
data_phylo_unmanaged.b.norm <- transform_sample_counts(data_phylo_unmanaged.b.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_unmanaged.b.rare)

# Calculate the number of individuals in the meta community (Average read depth)
N <- mean(apply(OTU.table, 1, sum))

# Calculate the average relative abundance of each taxa across communities
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p.df = data.frame(p) %>%
  rownames_to_column(var="OTU")

# Calculate the occurrence frequency of each taxa
OTU.table.bi <- 1*(OTU.table>0)
freq.table <- apply(OTU.table.bi, 2, mean)
freq.table <- freq.table[freq.table != 0]
freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

#Combine
C <- inner_join(p.df,freq.df, by="OTU") %>%
  arrange(p)
# Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
C.no0 <- C %>%
  filter(freq != 0, p != 0)

#Calculate the limit of detection
d <- 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
p.list <- C.no0$p
freq.list <- C.no0$freq
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.01))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats.unmanaged.b <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                                   Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                                   Samples=nrow(OTU.table), Richness=length(p.list), 
                                   Detect=d)

#################BELOW GROUND - HOME LAWN:#####################
# Rarefy to an even depth
set.seed(111)
data_phylo_home_lawn.b.rare <- rarefy_even_depth(data_phylo_home_lawn.b)

# Normalize read counts
data_phylo_home_lawn.b.norm <- transform_sample_counts(data_phylo_home_lawn.b.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_home_lawn.b.rare)

# Calculate the number of individuals in the meta community (Average read depth)
N <- mean(apply(OTU.table, 1, sum))

# Calculate the average relative abundance of each taxa across communities
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p.df = data.frame(p) %>%
  rownames_to_column(var="OTU")

# Calculate the occurrence frequency of each taxa
OTU.table.bi <- 1*(OTU.table>0)
freq.table <- apply(OTU.table.bi, 2, mean)
freq.table <- freq.table[freq.table != 0]
freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

#Combine
C <- inner_join(p.df,freq.df, by="OTU") %>%
  arrange(p)
# Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
C.no0 <- C %>%
  filter(freq != 0, p != 0)

#Calculate the limit of detection
d <- 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
p.list <- C.no0$p
freq.list <- C.no0$freq
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.01))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats.home_lawn.b <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                                   Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                                   Samples=nrow(OTU.table), Richness=length(p.list), 
                                   Detect=d)

#######################BELOW GROUND - GREENS:#######################
# Rarefy to an even depth
set.seed(111)
data_phylo_greens.b.rare <- rarefy_even_depth(data_phylo_greens.b)

# Normalize read counts
data_phylo_greens.b.norm <- transform_sample_counts(data_phylo_greens.b.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_greens.b.rare)

# Calculate the number of individuals in the meta community (Average read depth)
N <- mean(apply(OTU.table, 1, sum))

# Calculate the average relative abundance of each taxa across communities
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p.df = data.frame(p) %>%
  rownames_to_column(var="OTU")

# Calculate the occurrence frequency of each taxa
OTU.table.bi <- 1*(OTU.table>0)
freq.table <- apply(OTU.table.bi, 2, mean)
freq.table <- freq.table[freq.table != 0]
freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

#Combine
C <- inner_join(p.df,freq.df, by="OTU") %>%
  arrange(p)
# Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
C.no0 <- C %>%
  filter(freq != 0, p != 0)

#Calculate the limit of detection
d <- 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
p.list <- C.no0$p
freq.list <- C.no0$freq
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.01))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats.greens.b <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                                Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                                Samples=nrow(OTU.table), Richness=length(p.list), 
                                Detect=d)

###########################BELOW GROUND  - FAIRWAY:###################
# Rarefy to an even depth
set.seed(111)
data_phylo_fairway.b.rare <- rarefy_even_depth(data_phylo_fairway.b)

# Normalize read counts
data_phylo_fairway.b.norm <- transform_sample_counts(data_phylo_fairway.b.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_fairway.b.rare)

# Calculate the number of individuals in the meta community (Average read depth)
N <- mean(apply(OTU.table, 1, sum))

# Calculate the average relative abundance of each taxa across communities
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p.df = data.frame(p) %>%
  rownames_to_column(var="OTU")

# Calculate the occurrence frequency of each taxa
OTU.table.bi <- 1*(OTU.table>0)
freq.table <- apply(OTU.table.bi, 2, mean)
freq.table <- freq.table[freq.table != 0]
freq.df = data.frame(OTU=names(freq.table), freq=freq.table)

#Combine
C <- inner_join(p.df,freq.df, by="OTU") %>%
  arrange(p)
# Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
C.no0 <- C %>%
  filter(freq != 0, p != 0)

#Calculate the limit of detection
d <- 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
p.list <- C.no0$p
freq.list <- C.no0$freq
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.01))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats.fairway.b <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                                 Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                                 Samples=nrow(OTU.table), Richness=length(p.list), 
                                 Detect=d)

################################PLOT####################

#VALUES FROM THE FITSTATS DATA - Assuming you have Rsqr values for above and below habitats
Above_ground <- c(0.31472, 0.60403, 0.49735, 0.61221)  # Replace with your actual data
Below_ground <- c(0.46382, 0.60000, 0.21318, 0.42656)  # Replace with your actual data

# Assuming you have another variable
Management_intensity<- c("Unmanaged", "Home Lawn", "Greens", "Fairway")  # Replace with your actual data

# Create a data frame
data_df <- data.frame(
  Management_intensity = Management_intensity,
  Above_ground = Above_ground,
  Below_ground = Below_ground)

# Print the data frame
print(data_df)

# Reshape the data to long format for ggplot
library(tidyr)
data_long <- gather(data_df, key = "Type", value = "R2_value", -Management_intensity)

# Create a scatter plot with different colors for aerial and terrestrial
ITS_R2 <- ggplot(data_long, aes(x = Management_intensity, y = R2_value, color = Type)) +
  geom_point(size = 5) +
  labs(y = expression("Model fit ("~R^2~")"),  
       x = NULL,
       color = "") +
  scale_color_manual(values = c("orange", "purple")) +  # Set custom colors
  theme_bw(base_size = 20, base_family = "Times") +
  ylim(0, 1)   # Set y-axis limits from 0 to 1
ITS_R2

ITS_R2   <- ITS_R2 +
  theme(plot.margin = margin(40, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_ITS_R2   <- ggdraw() +
  draw_plot(ITS_R2) +  # Add the plot
  draw_label("Fungi", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("A", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_ITS_R2 

ggsave(file = "R2_Neutral_ITSF.tiff", dpi = 900, width = 12, height = 10, units = 'in')

#####END######################################################################################