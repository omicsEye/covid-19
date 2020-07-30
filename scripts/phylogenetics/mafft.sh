#!/bin/bash
#SBATCH -N 1
#SBATCH -t 7-00:00:00
#SBATCH -p defq
#SBATCH -o mafft_corona_big.out
#SBATCH -e mafft_corona_big.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rebeccaclement@gwu.edu

module load mafft

#--- Start the timer
t1=$(date +"%s")

# changed settings to go faster
mafft --thread $(nproc) --retree 2 --maxiterate 2 --nofft allseqs.fa > allseqsalnnofft.fa

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
