##########################################################################################################
##################8:SLOAN NEUTRAL MODEL - Figure 4-5 Whole community, below and above ground ##############
##########################################################################################################

### Clear workspace ###
citation("RColorBrewer")
packageVersion("RColorBrewer")

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

# Rarefy to an even depth
set.seed(111)  # setting seed for reproducibility
data_phylo.rare = rarefy_even_depth(data_phylo)

# Normalize read counts (this gives relative abundance)
data_phylo.norm = transform_sample_counts(data_phylo.rare, function(x) x/sum(x))

library("Hmisc")
library("minpack.lm")
library("stats4")

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo.rare)

####################################Method 1
                                          
# calculate the average number of individuals per community (number of sequences per sumple?) #
N <- mean(apply(OTU.table, 1, sum))

# calculate the relative abundance of each OTU across communities (across sample places) #
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N

# calculate the occurence frequency of each OTU across communities (across sample places) #
spp.bi <- 1*(OTU.table>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq !=0]

# combine the relative abundance and occurence frequency data #
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C <- C[!(apply(C, 1, function(y) any(y==0))),]
p <- C[,2]
freq <- C[,3]

# assign each OTU names to p and freq #
names(p) <- C[,1]
names(freq) <- C[,1]

# calcuation of the limit of detection #
d = 1/N

# estimation of migration parameter using Non-linear least squares #
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))

# get prediction values #
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)

# get 95% confidence interval using Wilson score #
pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)

# calculate goodness-of-fit (R-squared) #
Rsqr <- 1 - sum((freq - freq.pred)^2)/sum((freq - mean(freq))^2)

# arrange results #
B <- cbind(p, freq, freq.pred, pred.ci[,2:3])

# merge prediction results and keystone list #
#B <- merge(A, KSDF, by=0)
colnames(B) <- c("p", "freq", "freq.pred", "ci.lower", "ci.upper")

# write prediction data #
write.csv(B, "Result_Prediction_16S.tsv")

m.fit 
Rsqr

# set point colour #
Col <- rep(1,nrow(B))
B <- cbind(B, Col)

B[B$freq > B$ci.upper,]$Col <- 2
B[B$freq < B$ci.lower,]$Col <- 3

B[B$Col==1,]$Col <- "Neutral"
B[B$Col==2,]$Col <- "Above prediction"
B[B$Col==3,]$Col <- "Below prediction"

# create prediction line data #
C <- B[!duplicated(B$p),]

# Set up ggplot with aesthetics
p <- ggplot(B, aes(x=log10(p), y=freq, colour=Col))
p <- p + geom_point()

# Set manual scale for shapes and sizes
p <- p + scale_shape_manual(values=c(2,16))
p <- p + scale_size_manual(values=c(3,1))

# Set manual scale for colors
p <- p + scale_colour_manual(values=c("darkgreen", "orange", "purple"))

# Add lines to the plot
p <- p + geom_line(data=C, aes(x=log10(p), y=freq.pred), col="black", size=1.0)
p <- p + geom_line(data=C, aes(x=log10(p), y=ci.upper), col="black", lty=2, size=1.0)
p <- p + geom_line(data=C, aes(x=log10(p), y=ci.lower), col="black", lty=2, size=1.0)

# Add text annotation
p <- p + annotate("text", label="R^2*' = 0.614'",x=-4.8, y=0.95, size=5, parse=TRUE)

# Set theme options
p <- p + theme_bw(base_size = 16, base_family = "Times")
p <- p + theme(
  axis.text.x = element_text(size=rel(1.5)),
  axis.title.x = element_text(size=rel(1.5)),
  axis.text.y = element_text(size=rel(1.5)),
  axis.title.y = element_text(size=rel(1.5)),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  legend.position=c(0.05, 0.75),  # Position the legend at the bottom-right corner
  legend.justification=c(0.1,0.8),
  legend.text=element_text(size=10)
)

# Add axis labels
p <- p + xlab(expression(paste(log[10], "(Mean relative abundance)")))
p <- p + ylab("Occurence frequency")

# Hide guide legends for scale, size, and color
p <- p + guides(scale="none",
                size=FALSE,
                pch=guide_legend(title=NULL),
                colour=guide_legend(title=NULL))
p

# Increase the plot margins
bacterial.model  <- p  +
  theme(plot.margin = margin(25, 5, 5, 5, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_bacterial.model  <- ggdraw() +
  draw_plot(bacterial.model) +  # Add the plot
  draw_label("Bacteria", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("D", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_bacterial.model

ggsave(file = "Neutral_metacommunity_16SF.tiff", dpi = 900, width = 16, height = 12, units = 'in')

#########################################################
###########Figure 5: above and below ground##############
#########################################################

# Section is the column in your sample_info_tab that denotes the habitat
sample_info_tab$Habitat <- factor(ifelse(sample_info_tab$Section %in% c("Thatch", "Leaf"), "Above_ground", "Below_ground"))

# Update the phyloseq object with the new metadata
SAM <- sample_data(sample_info_tab, errorIfNULL=TRUE)
data_phylo <- phyloseq(OTU, TAX, SAM)

# Create phyloseq object for above habitat
data_phylo_above <- subset_samples(data_phylo, Habitat == "Above_ground")

################################ABOVE GROUND################################

# Rarefy to an even depth
set.seed(111)
data_phylo_above.rare <- rarefy_even_depth(data_phylo_above)

# Normalize read counts
data_phylo_above.norm <- transform_sample_counts(data_phylo_above.rare, function(x) x/sum(x))

library("Hmisc")
library("minpack.lm")
library("stats4")

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_above.rare)

# calculate the average number of individuals per community (number of sequences per sumple?) #
N <- mean(apply(OTU.table, 1, sum))

# calculate the relative abundance of each OTU across communities (across sample places) #
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N

# calculate the occurence frequency of each OTU across communities (across sample places) #
spp.bi <- 1*(OTU.table>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq !=0]

# combine the relative abundance and occurence frequency data #
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C <- C[!(apply(C, 1, function(y) any(y==0))),]
p <- C[,2]
freq <- C[,3]

# assign each OTU names to p and freq #
names(p) <- C[,1]
names(freq) <- C[,1]

# calcuation of the limit of detection #
d = 1/N

# estimation of migration parameter using Non-linear least squares #
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))

# get prediction values #
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)

# get 95% confidence interval using Wilson score #
pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)

# calculate goodness-of-fit (R-squared) #
Rsqr <- 1 - sum((freq - freq.pred)^2)/sum((freq - mean(freq))^2)

# arrange results #
B <- cbind(p, freq, freq.pred, pred.ci[,2:3])

# merge prediction results and keystone list #
#B <- merge(A, KSDF, by=0)
colnames(B) <- c("p", "freq", "freq.pred", "ci.lower", "ci.upper")

m.fit 
Rsqr

# set point colour #
Col <- rep(1,nrow(B))
B <- cbind(B, Col)

B[B$freq > B$ci.upper,]$Col <- 2
B[B$freq < B$ci.lower,]$Col <- 3

B[B$Col==1,]$Col <- "Neutral"
B[B$Col==2,]$Col <- "Above prediction"
B[B$Col==3,]$Col <- "Below prediction"

# create prediction line data #
C <- B[!duplicated(B$p),]

# Set up ggplot with aesthetics
p <- ggplot(B, aes(x=log10(p), y=freq, colour=Col))
p <- p + geom_point()

# Set manual scale for shapes and sizes
p <- p + scale_shape_manual(values=c(2,16))
p <- p + scale_size_manual(values=c(3,1))

# Set manual scale for colors
p <- p + scale_colour_manual(values=c("darkgreen", "orange", "purple"))

# Add lines to the plot
p <- p + geom_line(data=C, aes(x=log10(p), y=freq.pred), col="black", size=1.0)
p <- p + geom_line(data=C, aes(x=log10(p), y=ci.upper), col="black", lty=2, size=1.0)
p <- p + geom_line(data=C, aes(x=log10(p), y=ci.lower), col="black", lty=2, size=1.0)

# Add text annotation
p <- p + annotate("text", label="R^2*' = 0.668'",x=-4.8, y=0.95, size=5, parse=TRUE)

# Set theme options
p <- p + theme_bw(base_size = 16, base_family = "Times")
p <- p + theme(
  axis.text.x = element_text(size=rel(1.5)),
  axis.title.x = element_text(size=rel(1.5)),
  axis.text.y = element_text(size=rel(1.5)),
  axis.title.y = element_text(size=rel(1.5)),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  legend.position=c(0.05,0.75),  # Position the legend at the bottom-right corner
  legend.justification=c(0.1,0.8),
  legend.text=element_text(size=10)
)

# Add axis labels
p <- p + xlab(expression(paste(log[10], "(Mean relative abundance)")))
p <- p + ylab("Occurence frequency")

# Hide guide legends for scale, size, and color
p <- p + guides(scale="none",
                size=FALSE,
                pch=guide_legend(title=NULL),
                colour=guide_legend(title=NULL))
p

# Increase the plot margins
above16S.model  <- p  +
  theme(plot.margin = margin(25, 5, 5, 5, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_above16S.model  <- ggdraw() +
  draw_plot(above16S.model) +  # Add the plot
  draw_label("Above ground", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("E", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_above16S.model

ggsave(file = "Neutral_above_ground_16SF.tiff", dpi = 900, width = 16, height = 12, units = 'in')

##########################################BELOW GROUND#########################

# Create phyloseq object for below ground habitat
data_phylo_below <- subset_samples(data_phylo, Habitat == "Below_ground")
# Rarefy to an even depth
set.seed(111)
data_phylo_below.rare <- rarefy_even_depth(data_phylo_below)

# Normalize read counts
data_phylo_below.norm <- transform_sample_counts(data_phylo_below.rare, function(x) x/sum(x))

# Extract the OTU table from the phyloseq object
OTU.table = otu_table(data_phylo_below.rare)

# calculate the average number of individuals per community (number of sequences per sumple?) #
N <- mean(apply(OTU.table, 1, sum))

# calculate the relative abundance of each OTU across communities (across sample places) #
p.m <- apply(OTU.table, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N

# calculate the occurence frequency of each OTU across communities (across sample places) #
spp.bi <- 1*(OTU.table>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq !=0]

# combine the relative abundance and occurence frequency data #
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C <- C[!(apply(C, 1, function(y) any(y==0))),]
p <- C[,2]
freq <- C[,3]

# assign each OTU names to p and freq #
names(p) <- C[,1]
names(freq) <- C[,1]

# calcuation of the limit of detection #
d = 1/N

# estimation of migration parameter using Non-linear least squares #
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))

# get prediction values #
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)

# get 95% confidence interval using Wilson score #
pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)

# calculate goodness-of-fit (R-squared) #
Rsqr <- 1 - sum((freq - freq.pred)^2)/sum((freq - mean(freq))^2)

# arrange results #
B <- cbind(p, freq, freq.pred, pred.ci[,2:3])

# merge prediction results and keystone list #
#B <- merge(A, KSDF, by=0)
colnames(B) <- c("p", "freq", "freq.pred", "ci.lower", "ci.upper")

m.fit 
Rsqr

# set point colour #
Col <- rep(1,nrow(B))
B <- cbind(B, Col)

B[B$freq > B$ci.upper,]$Col <- 2
B[B$freq < B$ci.lower,]$Col <- 3

B[B$Col==1,]$Col <- "Neutral"
B[B$Col==2,]$Col <- "Above prediction"
B[B$Col==3,]$Col <- "Below prediction"

# create prediction line data #
C <- B[!duplicated(B$p),]

# Set up ggplot with aesthetics
p <- ggplot(B, aes(x=log10(p), y=freq, colour=Col))
p <- p + geom_point()

# Set manual scale for shapes and sizes
p <- p + scale_shape_manual(values=c(2,16))
p <- p + scale_size_manual(values=c(3,1))

# Set manual scale for colors
p <- p + scale_colour_manual(values=c("darkgreen", "orange", "purple"))

# Add lines to the plot
p <- p + geom_line(data=C, aes(x=log10(p), y=freq.pred), col="black", size=1.0)
p <- p + geom_line(data=C, aes(x=log10(p), y=ci.upper), col="black", lty=2, size=1.0)
p <- p + geom_line(data=C, aes(x=log10(p), y=ci.lower), col="black", lty=2, size=1.0)

# Add text annotation
p <- p + annotate("text", label="R^2*' = 0.595'",x=-4.8, y=0.95, size=5, parse=TRUE)

# Set theme options
p <- p + theme_bw(base_size = 16, base_family = "Times")
p <- p + theme(
  axis.text.x = element_text(size=rel(1.5)),
  axis.title.x = element_text(size=rel(1.5)),
  axis.text.y = element_text(size=rel(1.5)),
  axis.title.y = element_text(size=rel(1.5)),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  legend.position=c(0.05,0.75),  # Position the legend at the bottom-right corner
  legend.justification=c(0.1,0.8),
  legend.text=element_text(size=10)
)

# Add axis labels
p <- p + xlab(expression(paste(log[10], "(Mean relative abundance)")))
p <- p + ylab("Occurence frequency")

# Hide guide legends for scale, size, and color
p <- p + guides(scale="none",
                size=FALSE,
                pch=guide_legend(title=NULL),
                colour=guide_legend(title=NULL))
p

# Increase the plot margins
below16S.model  <- p  +
  theme(plot.margin = margin(25, 5, 5, 5, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_below16S.model  <- ggdraw() +
  draw_plot(below16S.model) +  # Add the plot
  draw_label("Below ground", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("F", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_below16S.model

ggsave(file = "Neutral_below_ground_16SF.tiff", dpi = 900, width = 16, height = 12, units = 'in')

#END
#################################################################################################
#################################################################################################
