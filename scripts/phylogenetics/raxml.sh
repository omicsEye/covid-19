#!/bin/bash
#SBATCH -J raxml0613a
#SBATCH -o raxml0613a.out
#SBATCH -e raxml0613a.err
#SBATCH -p defq,small-gpu
#SBATCH -n 16
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rebeccaclement@gwu.edu

module load raxml

#Marcos says: estimate Raxml tree you can run quick bootstrap analysis but for bbest ML he would use at least 10 random additions. Keith wants me to run ModelTest-NG with RAXML to get best-fit model
# --thread $(nproc)
#-N/# the number of alternative runs on distinct starting trees (autoMRE computes a maximum of 1000 BS searches)

# -T is number of threads, should be number of alignment patterns (found in info.runID file)/500. Pegasus has 40 cores total. 
#-f a (conduct a rapid Bootstrap analysis and search for the best-scoring ML tree in one single program run) -T (num threads)
#--- Start the timer
t1=$(date +"%s")


raxmlHPC \
 -f a \
 -m GTRGAMMAX \
 -s allseqs_dedup.fasta \
 -n corona0613a.tre \
 -p 42117 \
 -N 10 \
 -T 40 \
 -x 12345

#raxmlHPC ­f a ­x 12345 ­p 12345 ­# 100 ­m GTRCAT -s GISAID_SARS_BAT_MERS_alignvf.fasta -n corona.tree

t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
