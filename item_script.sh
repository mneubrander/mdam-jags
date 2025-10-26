#! /bin/bash
#SBATCH --account=dctrl-mn141
#SBATCH --partition=common
#SBATCH -o item_script.out
#SBATCH -e item_script.err
#SBATCH --mem=16G   # Keep 16G memory
#SBATCH -c 2

# Load modules
module load JAGS/4.3.2-rhel9
module load R/4.4.0

echo "--- JOB START ---"
echo "Job is running in directory:"
pwd
echo "---"

# @ - command line arguments
# We pass shell args ($@) to the R script
echo "Attempting to run with Rscript..."
Rscript jags_item-script.R "$@"
exit_code=$? # Capture the exit code immediately

# Check if the R script failed
if [ $exit_code -ne 0 ]; then
    echo "R script failed with exit code: $exit_code"
    echo "---"
    echo "Check item_script.err for the R error message."
    exit $exit_code
fi

echo "--- JOB SUCCESS ---"