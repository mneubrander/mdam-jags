#!/bin/bash

# --- SLURM DIRECTIVES ---
#SBATCH --account=dctrl-mn141
#SBATCH --partition=common
#SBATCH --array=1-50  # <-- !! IMPORTANT: This MUST match 4th argument !!
#SBATCH --mem=2G          
#SBATCH -c 2 # <-- !! IMPORTANT: Should match N_CHAINS              
#SBATCH -o logs/item_script_%A_%a.out  
#SBATCH -e logs/item_script_%A_%a.err  

# --- JOB PARAMETERS ---
# Arguments are now read from the sbatch command line
# $1: CSV_PATH
# $2: POP_PARAMS
# $3: MISS_PARAMS
# $4: N_SAMP
# $5: L
# $6: N_CHAINS
# $7: N_BURNIN
# $8: N_ITER
# ---

# --- SCRIPT START ---
# Create logs directory if it doesn't exist
mkdir -p logs

# Load modules
module load JAGS/4.3.2-rhel9
module load R/4.4.0

echo "--- JOB START ---"
echo "Main Job ID: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
# Use $4 to show N_SAMP
echo "Running sample $SLURM_ARRAY_TASK_ID of $4"
echo "---"

# Call your *new* R script
# Pass the command-line args $1-$8, and the SLURM task ID as the 9th
Rscript jags-item-script-array.R \
  "$1" \
  "$2" \
  "$3" \
  "$4" \
  "$5" \
  "$6" \
  "$7" \
  "$8" \
  "$SLURM_ARRAY_TASK_ID" # This is the 9th argument

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo "R script failed with exit code: $exit_code"
else
    echo "--- JOB SUCCESS ---"
fi

exit $exit_code