# metabar
This repository utilizes the DADA2 pipeline to assign taxonomic identification to metabarcoded environmental DNA (eDNA) samples with two distinct primer pairs. All aspects can be run locally EXCEPT the taxonomic model training, which requires substantial memory. Steps to run the pipeline are outlined below.  

Current data tested: *Romeiro-MiS*
Primer sets: *MarVer3, MiFish-E*

## instructions
BEFORE
    - download reference database from [midori](https://reference-midori.info/download.php#latest)
    - I used the unique COI fasta for DADA2 (must be a fasta)
    - name of file: `MIDORI2_UNIQ_NUC_SP_GB264_CO1_DADA2.fasta`
    - Download these R packages `DECIPHER` and `BiocParallel` (might need to force) and `dada2` 

1. Remove corrupted files from MarVer3 sequences 
`removeCorrupt.R` identified sequence files smaller than 1kb-- these are empty files

2. Run pipeline to sort files, filter and trim (uncomment these lines if have not run already), and proceed with error quantification, and parsing taxonomic reference database *MIDORI_CO1_DECIPHER_taxa_ref*
`metabarcoding.R` 

3. Need to setup advanced memory usage for training taxonomic classifier
```
taxa_ref <- LearnTaxa(dna_ref, names(dna_ref)) # This needs A LOT of memory
```


# Source
Jeremy J - Virginia Tech
