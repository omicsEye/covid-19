#!/bin/bash
#SBATCH -J raxml5
#SBATCH -o raxml5_autoMRE.out
#SBATCH -e raxml5.err
#SBATCH -p defq
#SBATCH -n 16
#SBATCH -t 14-00:00:00
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
 -s easyOutput.fasta \
 -n corona5.tre \
 -p 42117 \
 -N autoMRE \
 -T 40 \
 -x 12345

#raxmlHPC ­f a ­x 12345 ­p 12345 ­# 100 ­m GTRCAT -s GISAID_SARS_BAT_MERS_alignvf.fasta -n corona.tree

t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
