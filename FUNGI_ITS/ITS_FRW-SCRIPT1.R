#####################################################################################################
####################################1:PRE-PROCCESSING OF RAW DATA####################################
#####################################################################################################

### Clear workspace ###
rm(list=ls())

#Set the working directory#
getwd()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#Packages
BiocManager::install("dada2")
packageVersion("dada2")
BiocManager::install("S4Vectors")
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

#Seth path:
path_ITS <- "/Users/ITS_raw" #directory containing the fastq files
list.files(path_ITS) #list files in the directory

###Using only the Forward sequences since the reverse ones have bad quality and eliminates important information when used wit the forward.

fnFs <- sort(list.files(path_ITS, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path_ITS, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.namesF <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

#Identify primers
FWD <- "AACTTTYRRCAAYGGATCWCT"  ## 5.8S-Fun forward primer sequence 

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)

fnFs.filtN <- file.path(path_ITS, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory

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
plotQualityProfile(fnFs[43:48]) 
                              
#Filter
filtFs <- file.path(path_ITS, "Filtered", basename(fnFs))
#filtRs <- file.path(path.cut, "filtered", basename(cutRs))
names(filtFs) <- sample.namesF

system.time(out <- filterAndTrim(fnFs, filtFs, trimLeft = 23, truncLen=c(220),
                                 maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
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

write.table(seqtab.nochim,"ASVs_ITS_F.txt",sep=",")

saveRDS(seqtab.nochim, "ITS_F_sqtb.asv.rds") #save the seqtab table

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.namesF
head(track)

write.table(track, "read-count-trackingITSF.tsv", quote=FALSE, sep="\t", col.names=NA)

#####################################################################################################
####################################2:ASSIGN TAXONOMY################################################
#####################################################################################################

#FUNGI-UNITE REF DATABASE
unite.ref <- "/Users/sh_general_release_dynamic_all_25.07.2023_eukaryote.fasta" 

taxa_ITS <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa_ITS # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

head(taxa.print)

class(taxa_ITS)
class(taxa.print)
unname(taxa_ITS)

write.table(taxa_ITS,"ITS_Taxonomy_FRW.txt",sep=",") 

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
write(asv_fasta, "ASVs_ITS_FRW.fa")

# Count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_ITS_countsFORWARD.tsv", sep="\t", quote=F, col.names=NA)

#Taxa table:
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
ranks

colnames(taxa_ITS) <- ranks
rownames(taxa_ITS) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(taxa_ITS, "ASVs_taxa_ITS_FRW.tsv", sep = "\t", quote=F, col.names=NA)

taxonomy_table <- read.table("ASVs_taxa_ITS_FRW.tsv", sep="\t", header=T, row.names=1)
taxonomy_table

# Remove prefixes from each column
taxonomy_table[, -1] <- apply(taxonomy_table[, -1], 2, function(x) sub("^[a-z]__", "", x))

# Print the modified taxonomy data frame
print(taxonomy_table)

write.table(taxonomy_table, "ASVs_taxa_ITSFRW_NAs.tsv", sep = "\t", quote=F, col.names=NA) #some sequences have NAs at the Phylum level

#####################################################################################################
####################################4:REMOVE CONTAMINANTS############################################
#####################################################################################################

BiocManager::install("decontam")
library(decontam)
packageVersion("decontam")

colnames(asv_tab) 
vector_for_decontam <- c(rep(TRUE, 0), rep(FALSE, 48))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant)

#####################################################################################################
####################################5:READ DATA AND PHYLOSEQ OBJECT####################################
#####################################################################################################

# First, get rid of any row that have NA at the Phylum level

# Read the TSV file into a data frame
df <- read.delim("ASVs_taxa_ITSFRW_NAs.tsv", stringsAsFactors = FALSE)

# Store the number of rows before deletion
rows_before <- nrow(df)

# Subset the data frame to remove rows with "NA" in the Phylum column only
df_clean <- subset(df, !(Phylum == "NA"))

# Store the number of rows after deletion
rows_after <- nrow(df_clean)

# Calculate the number of rows deleted
rows_deleted <- rows_before - rows_after

# Print the number of rows deleted
print(paste("Rows deleted:", rows_deleted))

# Identify the row indices that were deleted
deleted_indices <- setdiff(1:rows_before, which(!(df$Phylum == "NA")))

# Print the indices of the rows that were deleted
print(paste("Indices of deleted rows:", paste(deleted_indices, collapse = ", ")))

# Print the cleaned data frame
print(df_clean)

# Write file with the cleaned data frame
write.table(df_clean, "ASVs_taxa_ITSFRW_USE.tsv", sep = "\t", quote=F, col.names=NA)

#Delete the indices from the ASV counts table

# Read the second TSV file into a data frame
df_second <- read.delim("ASVs_ITS_countsFORWARD.tsv", stringsAsFactors = FALSE)

# Define the indices of the rows you want to delete (replace with your actual indices)
indices_to_delete <- c(45, 166, 179, 180, 183, 196, 198, 203, 213, 240, 260, 293, 303, 319, 320, 335, 338, 342, 352, 364, 365, 372, 377, 397, 419, 428, 439, 446, 454, 457, 462, 470, 471, 473, 483, 492, 505, 512, 513, 518, 528, 530, 531, 539, 556, 565, 599, 600, 607, 618, 622, 627, 639, 649, 652, 655, 657, 670, 673, 674, 685, 691, 693, 695, 731, 732, 742, 746, 748, 753, 754, 779, 791, 793, 797, 804, 806, 808, 810, 815, 816, 818, 820, 821, 822, 828, 830, 832, 835, 840, 846, 857, 863, 869, 872, 875, 890, 899, 900, 909, 916, 924, 928, 938, 942, 943, 948, 952, 972, 982, 986, 991, 1014, 1018, 1024, 1030, 1036, 1037, 1038, 1054, 1060, 1062, 1066, 1070, 1071, 1072, 1077, 1087, 1090, 1096, 1100, 1114, 1120, 1122, 1123, 1128, 1133, 1140, 1141, 1152, 1156, 1160, 1164, 1168, 1169, 1175, 1182, 1189, 1202, 1204, 1206, 1215, 1222, 1233, 1238, 1239, 1246, 1255, 
                       1257, 1258, 1263, 1266, 1267, 1269, 1274, 1290, 1294, 1297, 1302, 1303, 1311, 1326, 1330, 1341, 1342, 1344, 1348, 1349, 1356, 1367, 1368, 1374, 1383, 1395, 1401, 1402, 1408, 1411, 1417, 1419, 1420, 1427, 1430, 1435, 1439, 1446, 1458, 1459, 1461, 1462, 1463, 1465, 1479, 1484, 1492, 1497, 1498, 1511, 1517, 1521, 1528, 1534, 1537, 1542, 1548, 1550, 1561, 1568, 1569, 1570, 1574, 1584, 1585, 1589, 1591, 1593, 1596, 1597, 1603, 1604, 1610, 1613, 1617, 1626, 1628, 1630, 1640, 1648, 1653, 1654, 1656, 1662, 1665, 1666, 1672, 1675, 1676, 1682, 1686, 1688, 1692, 1698, 1701, 1702, 1703, 1709, 1710, 1711, 1715, 1719, 1721, 1723, 1731, 1736, 1743, 1751, 1752, 1762, 1772, 1784, 1787, 1790, 1791, 1793, 1795, 1799, 1803, 1805, 1811, 1813, 1828, 1831, 1855, 1857, 1863, 1867, 1870, 1871, 1873, 1877, 1889, 1890, 1896, 1898, 1899, 1907, 1910, 
                       1917, 1923, 1927, 1928, 1934, 1937, 1938, 1939, 1942, 1957, 1962, 1966, 1968, 1970, 1971, 1972, 1979, 1987, 1988, 1994, 1996, 1999, 2003, 2011, 2014, 2016, 2019, 2023, 2025, 2028, 2029, 2031, 2033, 2035, 2036, 2039, 2044, 2047, 2052, 2056, 2061, 2068, 2069, 2073, 2077, 2083, 2085, 2086, 2087, 2091, 2096, 2104, 2109, 2111, 2112, 2119, 2121, 2122, 2123, 2126, 2143, 2146, 2147, 2149, 2155, 2159, 2160, 2167, 2168, 2175, 2186, 2190, 2197, 2201, 2205, 2225, 2228, 2238, 2243, 2245, 2249, 2251, 2275, 2276, 2278, 2280, 2282, 2283, 2288, 2289, 2290, 2292, 2295, 2311, 2315, 2316, 2322, 2331, 2336, 2339, 2341, 2342, 2347, 2349, 2350, 2351, 2354, 2365, 2369, 2377, 2380, 2385, 2386, 2387, 2391, 2394, 2404, 2407, 2409, 2414, 2422, 2423, 2434, 2449, 2457, 2459, 2474, 2479, 2485, 2487, 2494, 2495, 2500, 2501, 2503, 2505, 2508, 2509, 2511, 
                       2517, 2526, 2528, 2530, 2535, 2536, 2540, 2543, 2546, 2549, 2550, 2554, 2557, 2559, 2560, 2565, 2581, 2583, 2584, 2596, 2601, 2604, 2611, 2619, 2626, 2630, 2632, 2638, 2641, 2642, 2648, 2650, 2662, 2667, 2669, 2670, 2673, 2680, 2682, 2688, 2692, 2693, 2696, 2697, 2698, 2700, 2702, 2703, 2705, 2713, 2717, 2720, 2725, 2727, 2731, 2733, 2734, 2738, 2742, 2744, 2745, 2746, 2750, 2753, 2756, 2757, 2759, 2771, 2775, 2782, 2783, 2793, 2794, 2798, 2803, 2807, 2813, 2814, 2817, 2824, 2825, 2831, 2839, 2844, 2850, 2852, 2853, 2856, 2866, 2868, 2885, 2892, 2894, 2895, 2906, 2910, 2912, 2914, 2922, 2924, 2925, 2926, 2929, 2930, 2934, 2940, 2942, 2948, 2950, 2955, 2960, 2965, 2967, 2968, 2970, 2976, 2985, 2986, 2996, 3010, 3011, 3012, 3013, 3014, 3015, 3016, 3018, 3019, 3020, 3022, 3030, 3031, 3042, 3043, 3051, 3055, 3059, 3062, 3069, 
                       3075, 3079, 3080, 3084, 3088, 3094, 3097, 3100, 3108, 3113, 3117, 3118, 3125, 3135, 3136, 3138, 3143, 3145, 3148, 3152, 3161, 3173, 3182, 3183, 3186, 3188, 3191, 3193, 3194, 3199, 3206, 3212, 3215, 3224, 3234, 3236, 3237, 3238, 3240, 3242, 3248, 3249, 3254, 3255, 3257, 3258, 3259, 3261, 3262, 3263, 3272, 3274, 3276, 3289, 3291, 3294, 3302, 3305, 3309, 3310, 3312, 3321, 3322, 3326, 3332, 3337, 3339, 3340, 3343, 3346, 3355, 3359, 3360, 3365, 3370, 3371, 3373, 3377, 3379, 3380, 3381, 3385, 3387, 3399, 3400, 3402, 3403, 3404, 3405, 3411, 3436, 3437, 3445, 3449, 3452, 3457, 3458, 3460, 3461, 3462, 3473, 3477, 3478, 3485, 3487, 3488, 3493, 3496, 3501, 3508, 3511, 3512, 3523, 3531, 3537, 3538, 3540, 3542, 3543, 3544, 3551, 3555, 3564, 3567, 3570, 3573, 3575, 3579, 3581, 3582, 3587, 3588, 3591, 3594, 3602, 3610, 3617, 3619, 3620, 
                       3623, 3624, 3625, 3628, 3631, 3632, 3633, 3634, 3635, 3642, 3643, 3644, 3645, 3647, 3649, 3651, 3653, 3654, 3656, 3657, 3659, 3662, 3668, 3673, 3675, 3676, 3679, 3688, 3691, 3700, 3710, 3711, 3713, 3726, 3730, 3731, 3740, 3744, 3745, 3748, 3754, 3755, 3756, 3762, 3765, 3768, 3769, 3777, 3782, 3783, 3784, 3786, 3788, 3794, 3795, 3798, 3799, 3801, 3808, 3809, 3816, 3821, 3827, 3829, 3832, 3833, 3835, 3842, 3843, 3850, 3854, 3855, 3856, 3857, 3862, 3872, 3873, 3874, 3877, 3881, 3885, 3886, 3894, 3895, 3903, 3905, 3906, 3910, 3924, 3926, 3936, 3944, 3945, 3955, 3957, 3958, 3964, 3965, 3968, 3971, 3977, 3980, 3983, 3995, 3999, 4003, 4007, 4011, 4012, 4015, 4017, 4019, 4020, 4024, 4027, 4031, 4034, 4040, 4042, 4043, 4044, 4045, 4053, 4057, 4059, 4064, 4069, 4072, 4074, 4079, 4084, 4092, 4103, 4104, 4105, 4121, 4132, 4133, 4140, 
                       4141, 4142, 4143, 4145, 4147, 4162, 4163, 4168, 4170, 4182, 4184, 4186, 4187, 4190, 4191, 4194, 4197, 4200, 4206, 4211, 4214, 4223, 4228, 4231, 4233, 4234, 4242, 4245, 4246, 4248, 4253, 4260, 4269, 4274, 4276, 4281, 4283, 4289, 4291, 4293, 4298, 4299, 4304, 4305, 4307, 4309, 4310, 4311, 4313, 4314, 4315, 4317, 4323, 4325, 4326, 4327, 4330, 4331, 4335, 4337, 4342, 4345, 4346, 4352, 4354, 4356, 4357, 4360, 4361, 4369, 4373, 4374, 4376, 4378, 4380, 4386, 4387, 4388, 4398, 4399, 4401, 4402, 4405, 4410, 4412, 4421, 4426, 4430, 4432, 4437, 4441, 4443, 4460, 4461, 4466, 4470, 4472, 4477, 4479, 4480, 4482, 4493, 4502, 4503, 4506, 4510, 4518, 4529, 4530, 4536, 4537, 4538, 4541, 4542, 4544, 4547, 4548, 4551, 4552, 4553, 4564, 4565, 4567, 4570, 4580, 4591, 4593, 4594, 4604, 4606, 4607, 4611, 4612, 4613, 4620, 4623, 4628, 4633, 4640, 
                       4643, 4650, 4661, 4662, 4671, 4673, 4674, 4675, 4677, 4679, 4680, 4683, 4687, 4688, 4689, 4695, 4697, 4698, 4706, 4709, 4710, 4711, 4714, 4716, 4718, 4719, 4720, 4721, 4733, 4736, 4737, 4740, 4742, 4754, 4757, 4758, 4759, 4760, 4761, 4763, 4765, 4773, 4774, 4778, 4781, 4783, 4785, 4787, 4803, 4817, 4820, 4825, 4831, 4833, 4837, 4839, 4846, 4847, 4850, 4855, 4862, 4870, 4872, 4873, 4875, 4876, 4877, 4878, 4879, 4880, 4882, 4883, 4884, 4886, 4889, 4893, 4896, 4903, 4904, 4915, 4916, 4917, 4936, 4937, 4938, 4950, 4952, 4954, 4959, 4964, 4965, 4974, 4976, 4980, 4984, 4994, 4995, 4999, 5003, 5004, 5006, 5016, 5019, 5023, 5028, 5030, 5035, 5049, 5052, 5064, 5067, 5077, 5081, 5096, 5103, 5105, 5109, 5110, 5119, 5124, 5132, 5143, 5148, 5150, 5162, 5179, 5185, 5196, 5200, 5216, 5218, 5220, 5225, 5227, 5229, 5232, 5236, 5237, 5238, 
                       5240, 5243, 5255, 5256, 5257, 5258, 5260, 5264, 5267, 5271, 5272, 5273, 5279, 5280, 5281, 5284, 5286, 5290, 5291, 5311, 5312, 5313, 5314, 5319, 5320, 5321, 5324, 5325, 5329, 5330, 5333, 5334, 5344, 5346, 5347, 5352, 5353, 5362, 5365, 5366, 5374, 5378, 5381, 5382, 5383, 5388, 5390, 5391, 5393, 5395, 5396, 5404, 5405, 5411, 5414, 5419, 5420, 5422, 5424, 5427, 5432, 5443, 5448, 5452, 5453, 5456, 5459, 5464, 5469, 5475, 5479, 5480, 5482, 5484, 5489, 5494, 5495, 5496, 5498, 5503, 5508, 5524, 5525, 5529, 5531, 5539, 5542, 5545, 5549, 5551, 5552, 5556, 5557, 5565, 5567, 5569, 5572, 5577, 5582, 5585, 5592, 5593, 5596, 5611, 5621, 5622, 5628)  

# Remove rows with specified indices
df_clean_second <- df_second[-indices_to_delete, ]

# Print the cleaned data frame
print(df_clean_second)

# Write file with the cleaned data frame
write.table(df_clean_second, "ASVs_ITS_countsFORWARD_USE.tsv", sep = "\t", quote=F, col.names=NA)


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

#END 
#################################################################################################################################
#################################################################################################################################
