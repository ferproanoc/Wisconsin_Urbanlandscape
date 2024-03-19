#####################################################################################################
####################################1:PRE-PROCCESSING OF RAW DATA####################################
#####################################################################################################

#Set the working directory#
getwd()

### Clear workspace ###
rm(list=ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#Packages
BiocManager::install("dada2")
library(dada2)
packageVersion("dada2")
BiocManager::install("S4Vectors")
library(ShortRead)
library(Biostrings)

#Seth path:
path_16S <- "/Users/16s_raw" #directory containing the fastq files
list.files(path_16S) #list files in the directory

###Using only the Forward sequences since the reverse seqs bad quality and eliminates important information when used with the forward.

fnFs <- sort(list.files(path_16S, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path_ITS, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.namesF <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

#Identify primers
FWD <- "GTGCCAGCMGCCGCGGTAA"  ## 515F forward primer sequence 

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)

fnFs.filtN <- file.path(path_16S, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory

filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]))

#Quality - Inspect read quality profiles
plotQualityProfile(fnFs[1:6]) 
plotQualityProfile(fnFs[7:12]) 
plotQualityProfile(fnFs[13:18]) 
plotQualityProfile(fnFs[19:24]) 
plotQualityProfile(fnFs[25:30]) 
plotQualityProfile(fnFs[37:42]) 
plotQualityProfile(fnFs[43:47]) 

#Filter
filtFs <- file.path(path_16S, "Filtered", basename(fnFs))
#filtRs <- file.path(path.cut, "filtered", basename(cutRs))
names(filtFs) <- sample.namesF

system.time(out <- filterAndTrim(fnFs, filtFs, trimLeft =23, truncLen=c(245),
                                 maxN=0, maxEE=c(1), truncQ=2, rm.phix=TRUE,
                                 compress=TRUE, multithread=TRUE) )

head(out) #reads per sample 

#Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
#errR <- learnErrors(filtRs, multithread = TRUE)

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.namesF
#names(derepRs) <- sample.names

#Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#Sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

write.table(seqtab.nochim,"ASVs_16S_F.txt",sep=",")

saveRDS(seqtab.nochim, "16S_F_sqtb.asv.rds") #save the seqtab table

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.namesF
head(track)

write.table(track, "read-count-tracking16SF.tsv", quote=FALSE, sep="\t", col.names=NA) #tracking reads filtered

#####################################################################################################
####################################2:ASSIGN TAXONOMY################################################
#####################################################################################################

silva.ref <- "/Users/silva_nr99_v138.1_train_set.fa.gz" 

taxa_16S <- assignTaxonomy(seqtab.nochim, silva.ref, multithread = TRUE, tryRC = TRUE)
taxa_16S <- addSpecies(taxa_16S, "/Users/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa_16S # Removing sequence rownames for display only

rownames(taxa.print) <- NULL
head(taxa.print)

write.table(taxa_16S,"16S_Taxonomy_FRW.txt",sep=",") #taxonomy table without any filtering

#####################################################################################################
####################################3:EXTRACTING GOODS FROM DADA2####################################
#####################################################################################################

# Giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_16S_FRW.fa")

# Count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_16S_countsFORWARD.tsv", sep="\t", quote=F, col.names=NA)

#Taxa table:
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

colnames(taxa_16S) <- ranks
rownames(taxa_16S) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(taxa_16S, "ASVs_taxa_16S_FRW.tsv", sep = "\t", quote=F, col.names=NA)

#####################################################################################################
####################################4:REMOVE CONTAMINANTS############################################
#####################################################################################################

BiocManager::install("decontam")
library(decontam)
packageVersion("decontam")

colnames(asv_tab) 
vector_for_decontam <- c(rep(TRUE, 0), rep(FALSE, 47))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant)

##if any chloroplast or mitochondria contamination

is.chloroplast <- taxa_16S[,"Order"] %in% "Chloroplast"
seqtab.nochloro <- seqtab.nochim[,!is.chloroplast]
dim(seqtab.nochloro)

taxa.nochloro <- taxa_16S[!is.chloroplast,]
dim(taxa.nochloro)

is.mitochondria <- taxa.nochloro[,"Family"] %in% "Mitochondria"
seqtab.nomito <- seqtab.nochloro[,!is.mitochondria]
taxa.nomito <- taxa.nochloro[!is.mitochondria,]

seqtab_nochim_CM <- seqtab.nomito
dim(seqtab_nochim_CM)
taxonomy <- taxa.nomito
taxonomy

# Specify the directory path
output_directory <- "~/Documents/Koch_lab/1.R_projects/Urban_microbiome/Publication/R_analysis/"

# Create the directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)}

# Write the CSV file
write.csv(taxonomy, file.path(output_directory, "16S_noCM_taxonomy.tsv")) #taxonomy table without Mitochondria and Chloroplast

#####################################################################################################
####################################5:READ DATA AND PHYLOSEQ OBJECT####################################
#####################################################################################################

# Identify rows deleted during filtering for both chloroplast and mitochondria
deleted_rows <- taxa_16S[!(row.names(taxa_16S) %in% row.names(taxa.nochloro)) | !(row.names(taxa_16S) %in% row.names(taxa.nomito)), ]

# Write the deleted rows to a CSV file
write.csv(deleted_rows, file.path(output_directory, "deleted_taxa_rows.csv"))

# Display the deleted rows
deleted_rows

#Delete the indices from the ASV counts table

# Read the second TSV file into a data frame
df <- read.delim("ASVs_16S_countsFORWARD.tsv", stringsAsFactors = FALSE)

# Extract row indices from deleted_rows
deleted_indices <- as.data.frame(rownames(deleted_rows))
deleted_indices
# Extract numeric indices from row names by removing non-numeric characters
deleted_indices <- as.numeric(gsub("\\D", "", rownames(deleted_rows)))

# Print the result
print(deleted_indices)

# Convert the numeric indices to a comma-separated string with proper formatting
deleted_indices_str <- paste(deleted_indices, collapse = ", ")

# Print the result with commas between numbers
cat(deleted_indices_str)

# Define the indices of the rows you want to delete (replace with your actual indices)
indices_to_delete <- c(1, 2, 97, 317, 442, 455, 587, 603, 845, 890, 1095, 1227, 1254, 1322, 1364, 1476, 1581, 1620, 1766, 2031, 2090, 2153, 2191, 2250, 2300, 2338, 2361, 2392, 2408, 2499, 2612, 2759, 2805, 3004, 3015, 3120, 3183, 3202, 3297, 3305, 3431, 3433, 3453, 3492, 3496, 3541, 3709, 3711, 3747, 3764, 3794, 3818, 3930, 3993, 4046, 4220, 4359, 4392, 4403, 4410, 4454, 4507, 4509, 4534, 4547, 4566, 4621, 4683, 4857, 4918, 5006, 5027, 5029, 5051, 5168, 5252, 5301, 5334, 5338, 5342, 5426, 5427, 5435, 5501, 5616, 5730, 5788, 5829, 5837, 5936, 5974, 6002, 6100, 6135, 6160, 6186, 6187, 6206, 6250, 6259, 6326, 6472, 6589, 6640, 6782, 6789, 6908, 6981, 7088, 7138, 7243, 7244, 7265, 7323, 7428, 7450, 7453, 7481, 7483, 7489, 7510, 7515, 7517, 7518, 7692, 7725, 7882, 7887, 7891, 7901, 7933, 8117, 8122, 8130, 8171, 8175, 8405, 8406, 8425, 8482, 8595, 8646, 8654, 8709, 8771, 8909, 8919, 8920, 8939, 8956, 8959, 8997, 9212, 9216, 9220, 9232, 9254, 9255, 9331, 9549, 9552, 9561, 9627, 9635, 9647, 9675, 9693, 9703, 9839, 9899, 9999, 10017, 10055, 10180, 10231, 10303, 10304, 10313, 10320, 10321, 10346, 10350, 10354, 10422, 10428, 10548, 10685, 10699, 10721, 10767, 10771, 10772, 10778, 10833, 10846, 10931, 10932, 10946, 10947, 10949, 11051, 11086, 11116, 11208, 11300, 11301, 11346, 11367, 11456, 11480, 11481, 11482, 11483, 11495, 11496, 11503, 11504, 11505, 11508, 11566, 11594, 11603, 11634, 11749, 11961, 11975, 12027, 12047, 12057, 12102, 12108, 12113, 12114, 12133, 12144, 12146, 12154, 12256, 12260, 12276, 12409, 12497, 12541, 12560, 12626, 12686, 12711, 12713, 12820, 12832, 12834, 12854, 12863, 12864, 12872, 12877, 12895, 12936, 13077, 13087, 13092, 13105, 13190, 13317, 13330, 13434, 13449, 13456, 13460, 13480, 13482, 13493, 13494, 13498, 13511, 13512, 13535, 13543, 13560, 13686, 13692, 13698, 13700, 13705, 13716, 13804, 13805, 13886, 13887, 13891, 13899, 13902, 13986, 13987, 13993, 14033, 14035, 14059, 14087, 14163, 14290, 14323, 14348, 14355, 14405, 14406, 14426, 14438, 14465, 14468, 14469, 14498, 14528, 14546, 14563, 14675, 14678, 14685, 14686, 14695, 14699, 14707, 14709, 14730, 14732, 14757, 14758, 14760, 14761, 14766, 14769, 14774, 14775, 14776, 14802, 14850, 14969, 14977, 14992, 14994, 14997, 15082, 15084, 15101, 15140, 15143, 15172, 15194, 15197, 15301, 15350, 15352, 15358, 15389, 15431, 15433, 15463, 15483, 15497, 15498, 15587, 15620, 15637, 15640, 15643, 15656, 15657, 15658, 15659, 15709, 15741, 15764, 15767, 15791, 15793, 15843, 15892, 15923, 15961, 15963, 16002, 16003, 16004, 16005, 16026, 16027, 16028, 16032, 16033, 16058, 16059, 16062, 16070, 16075, 16078, 16079, 16080, 16081, 16082, 16085, 16098, 16117, 16118, 16150, 16189, 16195, 16233, 16236, 16240, 16245, 16321, 16322, 16325, 16332, 16335, 16336, 16346, 16384, 16389, 16393, 16447, 16450, 16453, 16456, 16460, 16463, 16533, 16633, 16636, 16640, 16641, 16684, 16686, 16729, 16735, 16766, 16768, 16769, 16772, 16773, 16811, 16813, 16814, 16815, 16818, 16889, 16892, 16949, 16954, 17017, 17018, 17019, 17031, 17035, 17086, 17131, 17133, 17184, 17186, 17210, 17240, 17241, 17242, 17243, 17247, 17268, 17271, 17273, 17299, 17330, 17347, 17348, 17349, 17374, 17420, 17431)
                      
# Remove rows with specified indices
df_clean <- df[-indices_to_delete, ]

# Print the cleaned data frame
print(df_clean)

# Write file with the cleaned data frame
write.table(df_clean, "ASVs_16s_countsFORWARD_USE.tsv", sep = "\t", quote=F, col.names=NA)

# Second, get rid of any row that have NA at the Phylum level

# Read the TSV file into a data frame
df <- read.delim("16S_taxonomy_use.tsv", stringsAsFactors = FALSE)

# Store the number of rows before deletion
rows_before <- nrow(df)
rows_before

# Subset the data frame to remove rows with "NA" in the Phylum column only
df_clean <- subset(df, !(Phylum == "NA"))

# Store the number of rows after deletion
rows_after <- nrow(df_clean)

# Calculate the number of rows deleted
rows_deleted <- rows_before - rows_after

# Identify the row indices that were deleted
deleted_indices <- setdiff(1:rows_before, which(!(df$Phylum == "NA")))

# Print the indices of the rows that were deleted
print(paste("Indices of deleted rows:", paste(deleted_indices, collapse = ", ")))

# Print the cleaned data frame
print(df_clean)

# Write file with the cleaned data frame
write.table(df_clean, "16S_ASVs_Taxonomy_F_USE.tsv", sep = "\t", quote=F, col.names=NA)

#Delete the indices from the ASV counts table

# Read the second TSV file into a data frame
df_second <- read.delim("ASVs_16S_countsFORWARD_USE.tsv", stringsAsFactors = FALSE)

# Define the indices of the rows you want to delete (replace with your actual indices)
indices_to_delete <- c(299, 968, 1717, 2059, 2103, 2182, 2365, 2394, 2469, 2520, 2585, 2606, 2760, 2915, 2922, 3228, 3547, 3668, 3707, 3795, 3824, 3929, 4021, 4058, 4432, 4502, 4526, 4706, 4726, 4896, 4904, 5162, 5182, 5237, 5386, 5749, 5751, 5774, 5805, 5961, 5975, 6174, 6249, 6432, 6551, 6803, 6842, 6843, 6868, 6916, 7006, 7012, 7226, 7243, 7267, 7344, 7418, 7456, 7485, 7486, 7561, 7601, 7702, 7789, 7931, 7961, 7987, 8062, 8096, 8119, 8218, 8288, 8304, 8309, 8315, 8388, 8550, 8577, 8688, 8771, 8807, 8857, 8903, 8999, 9064, 9124, 9158, 9174, 9197, 9247, 9252, 9265, 9267, 9302, 9348, 9392, 9424, 9440, 9456, 9459, 9472, 9481, 9500, 9509, 9511, 9536, 9550, 9625, 9648, 9738, 9806, 9835, 9903, 10048, 10123, 10146, 10155, 10235, 10364, 10390, 10449, 10452, 10501, 10569, 10586, 10638, 10639, 10651, 10670, 10742, 10785, 10787, 10861, 11044, 11144, 
                       11145, 11151, 11168, 11200, 11205, 11215, 11259, 11294, 11314, 11338, 11351, 11362, 11416, 11432, 11481, 11496, 11520, 11522, 11538, 11539, 11615, 11686, 11688, 11738, 11766, 11791, 11829, 11856, 11871, 11891, 11899, 11922, 11931, 11951, 11954, 11955, 11990, 11992, 12064, 12074, 12120, 12124, 12142, 12166, 12167, 12172, 12233, 12390, 12414, 12515, 12524, 12537, 12538, 12548, 12568, 12584, 12615, 12659, 12677, 12721, 12743, 12744, 12745, 12751, 12757, 12784, 12803, 12818, 12841, 12857, 12870, 12877, 12886, 12929, 12947, 12950, 12995, 13019, 13048, 13065, 13183, 13195, 13221, 13260, 13266, 13278, 13294, 13310, 13340, 13341, 13358, 13361, 13386, 13408, 13410, 13412, 13423, 13424, 13428, 13431, 13453, 13474, 13480, 13481, 13496, 13512, 13517, 13552, 13556, 13557, 13558, 13565, 13597, 13616, 13637, 13655, 13698, 13718, 13756, 
                       13783, 13793, 13794, 13805, 13824, 13829, 13838, 13848, 13869, 13899, 13931, 13939, 13955, 14001, 14025, 14056, 14077, 14105, 14123, 14149, 14168, 14179, 14194, 14214, 14240, 14250, 14251, 14264, 14267, 14288, 14289, 14322, 14323, 14331, 14337, 14338, 14340, 14366, 14387, 14405, 14407, 14411, 14427, 14446, 14450, 14465, 14469, 14471, 14474, 14478, 14503, 14515, 14524, 14537, 14538, 14554, 14558, 14566, 14598, 14601, 14603, 14631, 14632, 14648, 14676, 14683, 14691, 14695, 14696, 14698, 14712, 14714, 14728, 14744, 14745, 14758, 14802, 14821, 14822, 14825, 14827, 14829, 14833, 14836, 14842, 14860, 14862, 14863, 14872, 14880, 14921, 14996, 15008, 15026, 15031, 15063, 15077, 15089, 15097, 15141, 15155, 15168, 15179, 15180, 15265, 15275, 15276, 15277, 15279, 15281, 15284, 15293, 15303, 15304, 15320, 15322, 15339, 15341, 15360, 
                       15363, 15367, 15383, 15393, 15394, 15396, 15400, 15404, 15408, 15422, 15440, 15441, 15444, 15454, 15466, 15469, 15477, 15479, 15484, 15488, 15504, 15531, 15546, 15548, 15554, 15567, 15568, 15590, 15605, 15609, 15610, 15621, 15640, 15641, 15657, 15661, 15663, 15670, 15673, 15683, 15685, 15702, 15711, 15720, 15726, 15728, 15733, 15740, 15742, 15752, 15758, 15766, 15770, 15773, 15775, 15782, 15785, 15794, 15800, 15802, 15805, 15808, 15814, 15828, 15836, 15841, 15849, 15850, 15851, 15865, 15872, 15873, 15881, 15883, 15889, 15896, 15897, 15902, 15910, 15915, 15920, 15926, 15927, 15929, 15941, 15953, 15967, 15978, 15979, 15994, 15998, 16013, 16015, 16017, 16020, 16024, 16029, 16030, 16032, 16040, 16047, 16050, 16055, 16056, 16063, 16067, 16084, 16086, 16087, 16089, 16090, 16100, 16103, 16113, 16123, 16124, 16128, 16144, 16159, 
                       16163, 16164, 16165, 16183, 16184, 16194, 16196, 16198, 16203, 16209, 16211, 16216, 16255, 16256, 16257, 16263, 16269, 16271, 16287, 16300, 16301, 16322, 16326, 16329, 16333, 16342, 16345, 16352, 16363, 16377, 16380, 16399, 16400, 16402, 16404, 16411, 16415, 16423, 16437, 16444, 16445, 16452, 16453, 16456, 16458, 16467, 16474, 16478, 16490, 16503, 16504, 16505, 16514, 16517, 16521, 16532, 16547, 16559, 16567, 16578, 16581, 16588, 16599, 16614, 16616, 16630, 16631, 16632, 16636, 16647, 16648, 16649, 16650, 16652, 16670, 16680, 16687, 16708, 16727, 16728, 16732, 16733, 16734, 16746, 16751, 16759, 16779, 16780, 16787, 16794, 16803, 16811, 16814, 16833, 16837, 16843, 16844, 16858, 16897, 16905, 16916, 16917, 16928, 16953)

# Remove rows with specified indices
df_clean_second <- df_second[-indices_to_delete, ]

# Print the cleaned data frame
print(df_clean_second)

# Write file with the cleaned data frame
write.table(df_clean_second, "16S_ASVs_countsF_USE.tsv", sep = "\t", quote=F, col.names=NA)

###NOW THE ASVS COUNT TABLE AND THE TAXONOMY TABLE DONT INCLUDE ANY NA AT THE PHYLUM LEVEL AND ANY ASVS INDICE DELETED WAS ALSO DELETED IN THE TAXONOMY TABLE :AROUND 600 variants were filtered out!!!

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

#END 
#################################################################################################################################
#################################################################################################################################
