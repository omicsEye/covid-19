#!/bin/bash
#SBATCH -N 1
#SBATCH -t 14-00:00:00
#SBATCH -p defq
#SBATCH -o ../nate_apr23/modeltest_nate.out
#SBATCH -e ../nate_apr23/modeltest_nate.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stearrettnt@gwu.edu

module load modeltest-ng

cd /lustre/groups/cbi/Projects/covid19/scripts/

#--- Start the timer
t1=$(date +"%s")

# changed settings to go faster
modeltest-ng-static -i ../easyOutput.fasta -h uigf -f ef 

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
