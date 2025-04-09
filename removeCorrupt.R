# Set your path:
marver3_path <- "MarVer3"

# List all gzipped FASTQ files:
fastq_files <- list.files(marver3_path, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Check file sizes:
file_sizes <- file.info(fastq_files)$size

# Inspect sizes first:
data.frame(file=basename(fastq_files), size_bytes=file_sizes)

# Identify files smaller than 1 KB (1000 bytes):
small_files <- fastq_files[file_sizes < 1000]

# Inspect which files will be removed:
print(basename(small_files))

# Remove the small files after confirmation:
file.remove(small_files)
