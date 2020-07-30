#!/bin/bash
#SBATCH -J raxml_dist
#SBATCH -o /lustre/groups/cbi/Projects/covid19/nate_raxml/slurm_files/raxml_%A_%a.out
#SBATCH -e /lustre/groups/cbi/Projects/covid19/nate_raxml/slurm_files/raxml_%A_%a.err
#SBATCH -p defq
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stearrettnt@gwu.edu
#SBATCH --array=1-8

module load raxml

#Marcos says: estimate Raxml tree you can run quick bootstrap analysis but for bbest ML he would use at least 10 random additions. Keith wants me to run ModelTest-NG with RAXML to get best-fit model
# --thread $(nproc)
#-N/# the number of alternative runs on distinct starting trees (autoMRE computes a maximum of 1000 BS searches)

# -T is number of threads, should be number of alignment patterns (found in info.runID file)/500. Pegasus has 40 cores total. 
#-f a (conduct a rapid Bootstrap analysis and search for the best-scoring ML tree in one single program run) -T (num threads)
# -f x computes pairwise distances

#--- Start the timer
t1=$(date +"%s")

cd /lustre/groups/cbi/Projects/covid19/nate_raxml/

FILES=($(cat /lustre/groups/cbi/Projects/covid19/nate_raxml/raxml_names))

INPUT=${FILES[$SLURM_ARRAY_TASK_ID]}

echo ${FILES[$SLURM_ARRAY_TASK_ID]}

raxmlHPC \
 -f x \
 -# autoMRE \
 -m GTRGAMMAX \
 -s ${INPUT}.fasta \
 -n ${INPUT}_corona.tre \
 -p 42117 \
 -N 10 \
 -T 40 \
 -x 12345

#raxmlHPC ­f a ­x 12345 ­p 12345 ­# 100 ­m GTRCAT -s GISAID_SARS_BAT_MERS_alignvf.fasta -n corona.tree

t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
