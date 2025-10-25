#! /bin/bash
#SBATCH --account=dctrl-mn141
#SBATCH --partition=common
#SBATCH -o item_script.out
#SBATCH -e item_script.err
#SBATCH --mem-per-cpu=1G
#SBATCH -c 2

# Load R module (adjust for your cluster)
module load R/4.4.0

# @ - command line arguments
R CMD BATCH --args "$@" jags_item-script.R