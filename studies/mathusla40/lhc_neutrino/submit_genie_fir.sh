#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=rrg-mdiamond
#SBATCH --array=1-40
#SBATCH --mem=4G

# ---------------------------------------------------------------------------------------
# TEMPORARY DIRECTORIES and task ID in case running locally
if [ -z "$SLURM_TMPDIR" ]; then 
  echo "No SLURM_TMPDIR, running locally"
  SLURM_TMPDIR=$PATH_MG5_in
  SLURM_ARRAY_TASK_ID=0
fi

module load apptainer/1.3.5

singularity exec /home/tomren/data/bin/release/singularity.sif ./run_genie_fir.sh $SLURM_ARRAY_TASK_ID