home_dir=~/Box/COVID19_Project
regions=("Genome" "S" "E" "M" 'N' 'ORF1a' 
         'ORF1b' 'ORF1ab' 'ORF3a' 
         'ORF6' 'ORF7a' 'ORF7b' 
         'ORF8' "ORF10" "2'-O-ribose methyltransferase" 
         "3C-like proteinase" "3'-to-5' exonuclease" 
         "endoRNAse" "helicase" "leader protein" 
         "nsp2" "nsp3" "nsp4" "nsp6" 
         "nsp7" "nsp8" "nsp9" "nsp10" 
         "nsp11" "RNA-dependent RNA polymerase" 
          ) 
#regions=("ORF10") 
          
dist_method="K80"
for region in "${regions[@]}"
do
echo ${region}
# omeClust -i "${home_dir}/data/Distances_500_GISAID_With_NCBI/Distances/${dist_method}/${region}.tsv" -o "/Users/rah/Documents/m2clust_TN93/${region}"  --metadata "${home_dir}/data/Distances_500_GISAID_With_NCBI/Distances/Metadata.txt" --resolution low --linkage complete
omeClust -i "/Users/rah/Documents/Distances_2000_GISAID/${dist_method}/${region}.tsv" -o "/Users/rah/Documents/omeClust2000_${dist_method}/${region}"  --metadata "${home_dir}/data/Distances_2000_GISAID/Distances/GISAID_2000_Metadata.txt" --resolution low --linkage complete
done

# omeClust -i /Users/rah/Documents/Distances_2000_GISAID/TN93/Genome.tsv -o /Users/rah/Documents/omeClust2000_TN93/Genome_low  --metadata ~/Box/COVID19_Project/data/Distances_2000_GISAID/Distances/GISAID_2000_Metadata.txt --resolution low --linkage complete
# /Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/site-packages/omeClust/utilities.py:149: RuntimeWarning: invalid value encountered in double_scalars
#   s = (b - a) / max(a, b)
# The number of major clusters:  5
# Lineage  is the most influential metadata in clusters
# There are 5 clusters
# Output is written in /Users/rah/Documents/omeClust2000_TN93/Genome_low

# omeClust -i /Users/rah/Documents/Distances_2000_GISAID/TN93/Genome.tsv -o /Users/rah/Documents/omeClust2000_TN93/Genome_high  --metadata ~/Box/COVID19_Project/data/Distances_2000_GISAID/Distances/GISAID_2000_Metadata.txt --resolution high --linkage complete
# /Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/site-packages/omeClust/utilities.py:149: RuntimeWarning: invalid value encountered in double_scalars
#   s = (b - a) / max(a, b)
# /Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/site-packages/omeClust/utilities.py:149: RuntimeWarning: invalid value encountered in double_scalars
#   s = (b - a) / max(a, b)
# The number of major clusters:  5
# Lineage  is the most influential metadata in clusters
# There are 123 clusters
# Output is written in /Users/rah/Documents/omeClust2000_TN93/Genome_high

# omeClust -i /Users/rah/Documents/Distances_2000_GISAID/TN93/Genome.tsv -o /Users/rah/Documents/omeClust2000_TN93/Genome_medium  --metadata ~/Box/COVID19_Project/data/Distances_2000_GISAID/Distances/GISAID_2000_Metadata.txt --resolution medium --linkage complete
# /Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/site-packages/omeClust/utilities.py:149: RuntimeWarning: invalid value encountered in double_scalars
#   s = (b - a) / max(a, b)
# The number of major clusters:  5
# Lineage  is the most influential metadata in clusters
# There are 10 clusters
# Output is written in /Users/rah/Documents/omeClust2000_TN93/Genome_medium

