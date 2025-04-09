# Load libraries
library(dada2)
library(DECIPHER)
library(tidyverse)
library(ShortRead)

# Set the paths
marver3_path <- "MarVer3"
mifish_path <- "MiFish-E"

# List FASTQ files (assuming paired-end, adapt accordingly if single-end)
marver3_fnFs <- sort(list.files(marver3_path, pattern="_R1_001.fastq.gz", full.names=TRUE))
marver3_fnRs <- sort(list.files(marver3_path, pattern="_R2_001.fastq.gz", full.names=TRUE))

mifish_fnFs <- sort(list.files(mifish_path, pattern="_R1_001.fastq.gz", full.names=TRUE))
mifish_fnRs <- sort(list.files(mifish_path, pattern="_R2_001.fastq.gz", full.names=TRUE))

############ MarVer3 ############  ############  ############  ############  ############ 

# Create filtered directories
filt_marver3_F <- file.path(marver3_path, "filtered_F")
filt_marver3_R <- file.path(marver3_path, "filtered_R")
# dir.create(filt_marver3_F); dir.create(filt_marver3_R)

# Quality Control
# plotQualityProfile(marver3_fnFs[1:2])  # Forward reads
# plotQualityProfile(marver3_fnRs[1:2])  # Reverse reads

# Filtering and trimming
# out_marver3 <- filterAndTrim(marver3_fnFs, filt_marver3_F,
#                              marver3_fnRs, filt_marver3_R,
#                               # adjust based on quality profile
#                              maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                              compress=TRUE, multithread=TRUE)
# head(out_marver3)


# Learn errors
errF <- learnErrors(filt_marver3_F, multithread=TRUE)
errR <- learnErrors(filt_marver3_R, multithread=TRUE)

# Infer ASVs
dadaFs <- dada(filt_marver3_F, err=errF, multithread=TRUE)
dadaRs <- dada(filt_marver3_R, err=errR, multithread=TRUE)

# Merge pairs
mergers <- mergePairs(dadaFs, filt_marver3_F, dadaRs, filt_marver3_R, verbose=TRUE)

# Create sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Load MIDORI or your custom elasmobranch-specific database
# Assume your reference fasta file is already downloaded as "MIDORI2_UNIQ_NUC.fasta"
# https://reference-midori.info/download.php#latest
dna_ref <- readDNAStringSet("MIDORI2_UNIQ_NUC_SP_GB264_CO1_DADA2.fasta")

# Reformat headers to DECIPHER's required format
new_names <- sapply(names(dna_ref), function(x){
  # Extract taxonomy part after "tax="
  tax_part <- sub(".*tax=", "", x)
  tax_part <- gsub(",", ";", tax_part) # commas to semicolons
  new_header <- paste0("Root;", tax_part)
  new_header
})

# Check if reformatting worked correctly
head(new_names)

# Assign new names to dna_ref
names(dna_ref) <- new_names

taxa_ref <- LearnTaxa(dna_ref, names(dna_ref)) # This needs A LOT of memory

saveRDS(taxa_ref, file = "MIDORI_CO1_DECIPHER_taxa_ref.rds")

# Ensure your ASVs are loaded
asv_seqs <- DNAStringSet(getSequences(seqtab.nochim))

# Classify using your trained taxa_ref
asv_taxonomy <- IdTaxa(asv_seqs, trainingSet = taxa_ref, strand = "both", processors = NULL)

taxa_df <- t(sapply(asv_taxonomy, function(x) {
  ranks <- x$taxon
  ranks[is.na(ranks)] <- "unassigned"
  ranks
}))

colnames(taxa_df) <- c("Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Inspect:
head(taxa_df)

asv_counts <- data.frame(seqtab.nochim)
asv_tax_final <- cbind(ASV_ID = paste0("ASV", seq_len(ncol(asv_counts))),
                       taxa_df,
                       t(asv_counts))

# Inspect your final combined table:
head(asv_tax_final)

write.csv(asv_tax_final, "MarVer3_ASV_taxonomy_final.csv", row.names = FALSE)
