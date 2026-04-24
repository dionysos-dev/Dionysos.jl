cd problems/BipedRobot

export DIONYSOS_ROBOT_MODEL=6D
export DIONYSOS_NCHUNKS=20

jid=$(sbatch --parsable execution/slurm_array/compute_array.sbatch)
sbatch --dependency=afterok:$jid execution/slurm_array/merge_chunks.sbatch