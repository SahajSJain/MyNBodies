#!/bin/bash
#SBATCH -A rmittal3_gpu
#SBATCH --job-name=ACCGPU    # Job name
#SBATCH --time=02:00:00   # Total time
#SBATCH --partition=a100    # Choose partition
#SBATCH --qos=qos_gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-type=all   # Email condition
#SBATCH --mail-user=sjain43@jhu.edu    # Email address
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
export ACC_DEVICE_TYPE=nvidia
export ACC_DEVICE_NUM=0
echo $CUDA_VISIBLE_DEVICES
time ./mynbodies_acc_out > acc.log
