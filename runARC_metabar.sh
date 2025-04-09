#!/bin/bash
#SBATCH -J metabarcoding_MarVer3       # Job name
#SBATCH --account=sharkpulse           # Use your lab account allocation
#SBATCH --partition=normal_q           # Queue/Partition name
#SBATCH --nodes=1                      # Request 1 node
#SBATCH --ntasks-per-node=1            # 1 task per node
#SBATCH --cpus-per-task=32             # Use 32 CPU cores
#SBATCH --mem=128G                     # Allocate 128GB RAM
#SBATCH --time=24:00:00                # Set max time 24 hours
#SBATCH --output=out/output.log        # Output log file

# Load required modules
module reset
module load R

cd /projects/sharkpulse/metabar/
Rscript metabarcoding.R

echo "job completed successfully!"