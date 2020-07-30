#!/bin/bash
#SBATCH -N 1
#SBATCH -t 14-00:00:00
#SBATCH -p defq
#SBATCH -o iqtree.out
#SBATCH -e iqtree.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rebeccaclement@gwu.edu

module load iq-tree

#--- Start the timer
t1=$(date +"%s")

# changed settings to go faster
iqtree -s GISAID_SARS_BAT_MERS_alignvf.fasta -B 1000 -alrt 1000

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
