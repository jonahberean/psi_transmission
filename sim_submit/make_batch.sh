#!/bin/sh
# bash script for making batch.sh files for PENTrack simulations, and submitting them to sbatchsacct

for var in 005s_lp01_sqrt 005s_lp10_sqrt 005s_lp01_pow2 005s_miro_sqrt 020s_lp01_sqrt 020s_lp10_sqrt 020s_lp01_pow2 020s_miro_sqrt 100s_lp01_sqrt 100s_lp10_sqrt 100s_lp01_pow2 100s_miro_sqrt
do
    filename=85mm_norm_${var}

    cat >batch_${filename}.sh <<-EOF
#!/bin/sh
            
###SBATCH --time=200
#SBATCH --mem=4G
#SBATCH --array=0-50
#SBATCH --output=/home/jberean/slurmOutput/${filename}.%A_%a.out
#SBATCH --error=/home/jberean/slurmErrors/${filename}.%A_%a.out
#PBS -l walltime=10:00:00

ID=\$SLURM_ARRAY_TASK_ID
JOB=\$SLURM_ARRAY_JOB_ID

/home/jberean/PENTrack/PENTrack \$JOB\$ID /home/jberean/scratch/config_${filename}.in /home/jberean/scratch/${filename}_out/
EOF

sbatch batch_${filename}.sh
done