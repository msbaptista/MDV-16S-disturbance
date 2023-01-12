# Data from the manipulation experiments in 2015
# Sequenced at Waikato with IonTorrent in three different runs

# SRA accession number PRJNA527658

# Tutorial followed for the analysis
# https://benjjneb.github.io/dada2/tutorial.html


# Load library
library(dada2)


### JV_P1_MM Library ----


# Which directory
getwd()
#"C:/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/JV_P1_MM"

# Define as the path
path <- "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/JV_P1_MM" 

# read the names in the fastaq files 
fn <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# Inspect read quality profiles
plotQualityProfile(fn[1:2]) 

# IonExpress_77 and _91 have 672 reads and 296 reads - remove files 

# re-read the names in the fastaq files 
fn <- sort(list.files(path, pattern=".fastq", full.names = TRUE))


# Check if we still have primers ----

# Load the libraries to manipulate DNA strings ----
library(ShortRead)
packageVersion("ShortRead")

library(Biostrings)
packageVersion("Biostrings")


# Set up the primer sequences 
# where to check primer sequences - https://www.earthmicrobiome.org/protocols-and-standards/16s/

FWD <- "GTGYCAGCMGCCGCGGTAA"   ## 515F (Parada) 
REV <- "GGACTACNVGGGTWTCTAAT"  ## 806R (Apprill)


# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Primer orientations - forward / complement / reverse / reverse complement
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fns) {
  nhits <- vcountPattern(primer, sread(readFastq(fns)), fixed = FALSE)
  return(sum(nhits > 0))
}

# See if sequences have primers
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fns = fn[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fns = fn[[1]]))

#                  Forward Complement Reverse RevComp
# FWD.ForwardReads   55794          0       0       0
# REV.ForwardReads       0          0       0   17080



# Filter and trim ----

# Give files the sample names (not the IonTorrent Codes)
sample.names <- c("11M_C","15M_C","09M_C","17M_C","00B_C",
                  "07M","18M","23M","18M_C","20M_C",
                  "21M_C","13M_C","14M_C","13B","23M_C",
                  "11B","04M","15M","21M","07M_C",
                  "07B_C","11B_C","23B_C","14B","03B_C",
                  "11M","20M","15B")


# Place filtered files in filtered/ subdirectory
filt <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))


# Filter
out <- filterAndTrim(fn, filt, truncLen=265, trimLeft=19,
                     maxN=0, maxEE=2, truncQ=2,
                     compress=TRUE, multithread=TRUE)
out


# Check if trimming removed the primers ----


# Define as the path
path <- "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/JV_P1_MM/filtered"

# Read the names
fn.filt <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))

# Look at primer detection for the a number of samples
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fns = fn.filt[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fns = fn.filt[[1]]))

#                 Forward Complement Reverse RevComp
# FWD.ForwardReads       1          0       0       0
# REV.ForwardReads       0          0       0       1


# We have removed the primers 


# Learn the Error Rates ----
err <- learnErrors(filt, multithread=TRUE)

# Visualize the estimated error rates
plotErrors(err, nominalQ=TRUE)

# The plots don't look so good, but the quality plots did not look so good either


# Sample Inference ----
dada1 <- dada(filt, err=err, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

# Construct sequence table
seqtab1 <- makeSequenceTable(dada1)

dim(seqtab1)
# [1]  28 1282 

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab1)))
# 246 

saveRDS(seqtab1, "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/seqtab1.rds")


### JV_P2_MM Library ----

# Which directory
getwd()
#"C:/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/JV_P2_MM"

# Define as the path
path <- "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/JV_P2_MM" 

# read the names in the fastaq files 
fn <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# Inspect read quality profiles
plotQualityProfile(fn[5:6]) 


# Check if we still have primers ----

# See if sequences have primer - in this example sample #1, change to check others
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fns = fn[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fns = fn[[1]]))

#                  Forward Complement Reverse RevComp
# FWD.ForwardReads   85607          0       0       0
# REV.ForwardReads       0          0       0   44744



# Filter and trim ----

# Give files the sample names (not the IonTorrent Codes)
sample.names <- c("00B","04B","07B","23B",
                  "18B","03M","13M","02M_C","02B",
                  "03B","05B","06B","09B","21B",
                  "10B","22B","17B","02M","19B",
                  "09M","17M","05M","05M_C","14M",
                  "03M_C")

# Place filtered files in filtered/ subdirectory
filt <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))


# Filter
out <- filterAndTrim(fn, filt, truncLen=265, trimLeft=19,
                     maxN=0, maxEE=2, truncQ=2,
                     compress=TRUE, multithread=TRUE)
out


# Check if trimming removed the primers ---


# Define as the path
path <- "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/JV_P2_MM/filtered"

# Read the names
fn.filt <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))

# Look at primer detection for the a number of samples
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fns = fn.filt[[10]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fns = fn.filt[[10]]))

#                 Forward Complement Reverse RevComp
# FWD.ForwardReads       1          0       0       0
# REV.ForwardReads       0          0       0       9


# We have removed the primers 


# Learn the Error Rates ----
err <- learnErrors(filt, multithread=TRUE)

# Visualize the estimated error rates
plotErrors(err, nominalQ=TRUE)

# The plots don't look so good, but the quality plots did not look so good either


# Sample Inference ----
dada2 <- dada(filt, err=err, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

# Construct sequence table
seqtab2 <- makeSequenceTable(dada2)

dim(seqtab2)
# [1]  25 1933 

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))
# 246 

saveRDS(seqtab2, "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/seqtab2.rds")

# Double the size of previous run, but better part of the P1 library are cDNA whereas better part of P2 library are DNA


### Maria_R13 Library ----

# Which directory
getwd()
#"C:/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/Maria_R13"

# Define as the path
path <- "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/Maria_R13" 

# read the names in the fastaq files 
fn <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# Inspect read quality profiles
plotQualityProfile(fn[7:8]) 


# Check if we still have primers ----

# See if sequences have primer - in this example sample #1, change to check others
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fns = fn[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fns = fn[[1]]))

#                  Forward Complement Reverse RevComp
# FWD.ForwardReads   60549          0       0       0
# REV.ForwardReads       0          0       0   25823



# Filter and trim ----

# Give files the sample names (not the IonTorrent Codes)
sample.names <- c("0M","22M","10M_C","6M","10M",
                  "0M_C","6M_C","22M_C")


# Place filtered files in filtered/ subdirectory
filt <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))


# Filter
out <- filterAndTrim(fn, filt, truncLen=265, trimLeft=19,
                     maxN=0, maxEE=2, truncQ=2,
                     compress=TRUE, multithread=TRUE)
out


# Check if trimming removed the primers ---


# Define as the path
path <- "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/Maria_R13/filtered"

# Read the names
fn.filt <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))

# Look at primer detection for the a number of samples
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fns = fn.filt[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fns = fn.filt[[1]]))

#                 Forward Complement Reverse RevComp
# FWD.ForwardReads       1          0       0       0
# REV.ForwardReads       0          0       0       3


# We have removed the primers 


# Learn the Error Rates ----
err <- learnErrors(filt, multithread=TRUE)

# Visualize the estimated error rates
plotErrors(err, nominalQ=TRUE)

# The plots don't look so good, but the quality plots did not look so good either


# Sample Inference ----
dada3 <- dada(filt, err=err, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

# Construct sequence table
seqtab3 <- makeSequenceTable(dada3)

dim(seqtab3)
# [1]  8 1038 

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab3)))
# 246 

saveRDS(seqtab3, "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/seqtab3.rds")


### MERGE LIBRARIES ----


# Merge the runs - https://benjjneb.github.io/dada2/bigdata_paired.html
st.all <- mergeSequenceTables(seqtab1, seqtab2, seqtab3)

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab)
# [1]   61  2925

sum(seqtab)/sum(st.all)
# [1]  0.9899162 - the frequency of non-chimeric sequences
# chimeras account for ca 2% of the abundance

saveRDS(seqtab, "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/seqtab.rds")

# Create an ASV table
write.csv(seqtab, "asv_dry_valleys.csv")


# Create fasta table 
library(seqinr)

name.seq <- paste0("ASV", seq(ncol(seqtab)))
uniqueSeqs <- as.list(colnames(seqtab))
write.fasta(uniqueSeqs, name.seq, "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/uniqueSeqs_dry_valleys.fasta")


#### ASSIGN TAXONOMY ----


# DADA2 native implementation of the naive Bayesian classifier

### SILVA

taxa_SILVA <- assignTaxonomy(seqtab, "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
tax_SILVA <- addSpecies(taxa_SILVA, "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/silva_species_assignment_v138.fa.gz")
saveRDS(tax_SILVA, "tax_SILVA.rds") 
write.csv(tax_SILVA, "taxa_SILVA_dry_valleys.csv")


### GTDB

tax_GTDB <- assignTaxonomy(seqtab, "/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/GTDB_bac-arc_ssu_r86.fa", multithread=T) 
saveRDS(tax_GTDB, "tax_GTDB.rds") 
write.csv(tax_GTDB, "taxa_GTDB_dry_valleys.csv")

# IDTAXA 

library(DECIPHER)

dna <- DNAStringSet(getSequences(seqtab))

### SILVA

load("/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/SILVA_SSU_r138_2019.RData")
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab)
saveRDS(taxid, "tax_ID_SILVA.rds") 
write.csv(taxid, "taxa_ID_SILVA_dry_valleys.csv")

### GTDB

load("/Users/Mafalda/Documents/Artigo_Catarina_Dr_Valleys/Raw Sequences/GTDB_r95-mod_August2020.RData")
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab)
saveRDS(taxid, "tax_ID_GTDB.rds") 
write.csv(taxid, "taxa_ID_GTDB_dry_valleysb.csv")


#### TREE ----

library(phangorn)

### Align sequences
seqs <- getSequences(seqtab)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)


#Use only NNI rearrangements first
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 1))
saveRDS(fitGTR, "fitGTR.rds")

#This is only 5 times instead of (a maximum) of 100 rearrangements
fitGTR <- optim.pml(fitGTR, rearrangement = "stochastic", ratchet.par = list(iter = 5L, maxit = 5L, prop = 1/3), multicore = TRUE)
saveRDS(fitGTR, "fitGTR.rds")

# and in the end optimize all the other parameters again (they should not change much)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 1))

saveRDS(fitGTR, "fitGTR.rds")


#### HAND DOWN TO PHYLOSEQ ----

library(phyloseq)

seqtab<- readRDS("seqtab.rds")
tax_SILVA <- readRDS("tax_SILVA.rds")


# Trim Sample names
row.names(seqtab)
row.names(seqtab) <- sub("_filt.fastq.gz", "", row.names(seqtab))


# Import metadata

# Create phyloseq
dryvalleys <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
                       tax_table(tax_SILVA),
                       phy_tree(fitGTR$tree))



saveRDS(dryvalleys, file = "dryvalleys.rds")


# Rename ASVs sequence to something shorter

# Inspect the data
otu_table(dv)[1:10, 1:10]

# Rename
new.names <- paste0("ASV_", seq(ntaxa(dv))) # define new names ASV_1, ASV_2, ...
seqs <- taxa_names(dv) # store sequences
names(seqs) <- new.names # make map from ASV to full sequence
taxa_names(dv) <- new.names # rename 

# Re-inspect the data  check it is renamed
otu_table(dv)[1:10, 1:10]



# Now lets add the metadata to this phyloseq object. For that we need the samples names to what match the metadata

# Lets check the samples names on the asv table
sample_names(dv)


# Load the metadata and environmental data
envdata <-  read.csv("manipulation_dv_metadata.csv", row.names = 1, header = TRUE)

# Lets check the samples names on the metadata table
row.names(envdata)

# The names don't match. Lets create the sample names in the order that we want for the asv table 
sample.names <- c("11M_C", "15M_C", "09M_C", "17M_C", "00B_C", "07M",   "18M",   "23M",   "18M_C", "20M_C", "21M_C",
                  "13M_C", "14M_C", "13B",   "23M_C", "11B",   "04M",   "15M",   "21M",   "07M_C", "07B_C", "11B_C",
                  "23B_C", "14B",   "03B_C", "11M",   "20M",   "15B",   "00B",   "04B",   "07B",   "23B",   "18B",
                  "03M",   "13M",   "02M_C", "02B",   "03B",   "05B",   "06B",   "09B",   "21B",   "10B",   "22B",
                  "17B",   "02M",   "19B",   "09M",   "17M",   "05M",   "05M_C", "14M",   "03M_C", "00M",    "22M",
                  "10M_C", "06M",    "10M",   "00M_C", "06M_C",  "22M_C")


# Create a data frame
asvmat <- data.frame(otu_table(dv))

# Rename the samples
row.names(asvmat)<- paste0(sample.names)

# Create a otu_table
asvs <- otu_table(asvmat, taxa_are_rows = FALSE)

# Create a new phyloseq
dv2 <-  phyloseq(otu_table(asvmat, taxa_are_rows = FALSE),
                 tax_table(dv),
                 sample_data(envdata),
                 phy_tree(dv))

dv2

saveRDS(dv2, file="dryvalleys_manipulation.rds")

