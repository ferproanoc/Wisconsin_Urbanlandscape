###############################################################################################################################################
#########################################################ITS_DATASET R ANALYSIS-4###################################################################
##############################################################################################################################################
######Fer Proano Cuenca #################################################################FEB 2024#############################################

##########################################################################################################
##################8:SLOAN NEUTRAL MODEL - Figure 4-5 Whole community, below and above ground ##############
##########################################################################################################

#patterns of biodiversity within microbial communities, particularly in the context of microbial biogeography. 
#https://github.com/seb369/landuse_comm_assembly/blob/master/Neutral_model.md

### Clear workspace ###
rm(list=ls())

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
write.csv(B, "Result_Prediction_ITSS.tsv")

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
p <- p + annotate("text", label="R^2*' = 0.312'",x=-5, y=0.95, size=5, parse=TRUE)

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
fungal.model  <- p  +
  theme(plot.margin = margin(25, 5, 5, 5, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_fungal.model  <- ggdraw() +
  draw_plot(fungal.model) +  # Add the plot
  draw_label("Fungi", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("A", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_fungal.model

ggsave(file = "Neutral_metacommunity_ITSF.tiff", dpi = 900, width = 16, height = 12, units = 'in')

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
p <- p + annotate("text", label="R^2*' = 0.139'",x=-4.8, y=0.95, size=5, parse=TRUE)

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
aboveITS.model  <- p  +
  theme(plot.margin = margin(25, 5, 5, 5, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_aboveITS.model  <- ggdraw() +
  draw_plot(aboveITS.model) +  # Add the plot
  draw_label("Above ground", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("B", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_aboveITS.model

ggsave(file = "Neutral_above_ground_ITSF.tiff", dpi = 900, width = 16, height = 12, units = 'in')

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
p <- p + annotate("text", label="R^2*' = 0.240'",x=-4.8, y=0.95, size=5, parse=TRUE)

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
belowITS.model  <- p  +
  theme(plot.margin = margin(25, 5, 5, 5, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_belowITS.model  <- ggdraw() +
  draw_plot(belowITS.model) +  # Add the plot
  draw_label("Below ground", fontface = "bold", size = 20, x = 0.5, y = 0.98) +  # Add the title
  draw_label("C", size = 20, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_belowITS.model

ggsave(file = "Neutral_below_ground_ITSF.tiff", dpi = 900, width = 16, height = 12, units = 'in')

#END#############################################################################################################
####################################Method 2

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
m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.1))
m.ci <- confint(m.fit, 'm', level=0.95)
m.sum <- summary(m.fit)
m.coef = coef(m.fit)

freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))

# Get table of model fit stats
fitstats <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                       Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                       Samples=nrow(OTU.table), Richness=length(p.list), 
                       Detect=d)

# Get confidence interval for predictions
freq.pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)

# Get table of predictions
pred.df <- data.frame(metacomm_RA=p.list, frequency=freq.pred, 
                      frequency_lowerCI=freq.pred.ci[,2], 
                      frequency_upperCI=freq.pred.ci[,3]) %>%
  unique()

# Get table of observed occupancy and abundance
obs.df <- C.no0 %>%
  dplyr::rename(metacomm_RA = p, frequency = freq)

fungal.model <- ggplot(data=obs.df) +
  geom_point(data=obs.df, aes(x=log10(metacomm_RA), y=frequency), 
             alpha=.4, size=3, color="#37AEC3") +
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency), color="red") + 
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_lowerCI), linetype=2, color="red") + 
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_upperCI), linetype=2, color="red") + 
  geom_text(data=fitstats, aes(label = paste("R^2 == ", round(Rsqr, 3))), 
            x=-5.75, y=0.96, size=5, parse=TRUE) +
  geom_text(data=fitstats, aes(label = paste("italic(m) ==", round(m, 3))), 
            x=-5.75, y=0.88, size=5, parse=TRUE) + 
  labs(x = expression(Log[10](Mean~relative~abundance)), y="Ocurrence frequency") +
  theme_bw(base_size = 32, base_family = "Times") +
  theme(legend.position = "none",axis.title = element_text(size=20),axis.text = element_text(size=20))+
  ylim(0,1)  # Set the x-axis limits
fungal.model 

# Increase the plot margins
fungal.model  <- fungal.model  +
  theme(plot.margin = margin(30, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_fungal.model  <- ggdraw() +
  draw_plot(fungal.model) +  # Add the plot
  draw_label("ITS Neutral Model", fontface = "bold", size = 18, x = 0.5, y = 0.98) +  # Add the title
  draw_label("A", size = 18, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_fungal.model

ggsave(file = "Neutral_metacommunity_ITSF.tiff", dpi = 900, width = 16, height = 12, units = 'in')


#######above#####

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
fitstats.a <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                         Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                         Samples=nrow(OTU.table), Richness=length(p.list), 
                         Detect=d)

# Get confidence interval for predictions
freq.pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)

# Get table of predictions
pred.df <- data.frame(metacomm_RA=p.list, frequency=freq.pred, 
                      frequency_lowerCI=freq.pred.ci[,2], 
                      frequency_upperCI=freq.pred.ci[,3]) %>%
  unique()

# Get table of observed occupancy and abundance
obs.df <- C.no0 %>%
  dplyr::rename(metacomm_RA = p, frequency = freq)

#plot

aboveITS.model = ggplot(data=obs.df) +
  geom_point(data=obs.df, aes(x=log10(metacomm_RA), y=frequency), 
             alpha=.4, size=3, color="#37AEC3") +
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency), color="red") + 
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_lowerCI), linetype=2, color="red") + 
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_upperCI), linetype=2, color="red") + 
  geom_text(data=fitstats.a, aes(label = paste("R^2 == ", round(Rsqr, 3))), 
            x=-5.25, y=0.96, size=5, parse=TRUE) +
  geom_text(data=fitstats.a, aes(label = paste("italic(m) ==", round(m, 3))), 
            x=-5.25, y=0.88, size=5, parse=TRUE) + 
  labs(x = expression(Log[10](Mean~relative~abundance)), y="Ocurrence frequency") +
  theme_bw(base_size = 32, base_family = "Times") +
  theme(legend.position = "none",axis.title = element_text(size=20),axis.text = element_text(size=20))+
  ylim(0,1)  # Set the x-axis limits

aboveITS.model  <- aboveITS.model +
  theme(plot.margin = margin(30, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_aboveITS.model  <- ggdraw() +
  draw_plot(aboveITS.model) +  # Add the plot
  draw_label("ITS Above ground", fontface = "bold", size = 18, x = 0.5, y = 0.98) +  # Add the title
  draw_label("A", size = 18, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_aboveITS.model

ggsave(file = "Neutral_above_ground_ITSF.tiff", dpi = 900, width = 16, height = 12, units = 'in')

#######below

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
fitstats.b <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2], 
                         Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N, 
                         Samples=nrow(OTU.table), Richness=length(p.list), 
                         Detect=d)

# Get confidence interval for predictions
freq.pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)

# Get table of predictions
pred.df <- data.frame(metacomm_RA=p.list, frequency=freq.pred, 
                      frequency_lowerCI=freq.pred.ci[,2], 
                      frequency_upperCI=freq.pred.ci[,3]) %>%
  unique()

# Get table of observed occupancy and abundance
obs.df <- C.no0 %>%
  dplyr::rename(metacomm_RA = p, frequency = freq)

#plot

belowITS.model <- ggplot(data=obs.df) +
  geom_point(data=obs.df, aes(x=log10(metacomm_RA), y=frequency), 
             alpha=.4, size=3, color="#37AEC3") +
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency), color="red") + 
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_lowerCI), linetype=2, color="red") + 
  geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_upperCI), linetype=2, color="red") + 
  geom_text(data=fitstats.b, aes(label = paste("R^2 == ", round(Rsqr, 3))), 
            x=-5.25, y=0.96, size=5, parse=TRUE) +
  geom_text(data=fitstats.b, aes(label = paste("italic(m) ==", round(m, 3))), 
            x=-5.25, y=0.88, size=5, parse=TRUE) + 
  labs(x = expression(Log[10](Mean~relative~abundance)), y="Ocurrence frequency") +
  theme_bw(base_size = 32, base_family = "Times") +
  theme(legend.position = "none",axis.title = element_text(size=20),axis.text = element_text(size=20))+
  ylim(0,1)  # Set the x-axis limits

belowITS.model  <- belowITS.model +
  theme(plot.margin = margin(30, 10, 10, 10, "pt"))  # Adjust the margins as needed

# Create a ggdraw object
plot_belowITS.model  <- ggdraw() +
  draw_plot(belowITS.model) +  # Add the plot
  draw_label("ITS Below ground", fontface = "bold", size = 18, x = 0.5, y = 0.98) +  # Add the title
  draw_label("B", size = 18, x = 0.025, y = 0.98, fontface = "bold")  # Add the tag "A"

# Print the plot
plot_belowITS.model

ggsave(file = "Neutral_below_ground_ITSF.tiff", dpi = 900, width = 16, height = 12, units = 'in')
