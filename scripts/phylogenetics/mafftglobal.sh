#!/bin/bash
#SBATCH -N 1
#SBATCH -t 14-00:00:00
#SBATCH -p defq,short
#SBATCH -o mafft_coronaglobal.out
#SBATCH -e mafft_coronaglobal.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rebeccaclement@gwu.edu

module load mafft

#--- Start the timer
t1=$(date +"%s")

# changed settings to go faster
mafft --thread $(nproc) --retree 2 --maxiterate 0 --globalpair allseqs.fa > allseqs_aln.fa

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
