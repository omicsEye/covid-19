home_dir=~/Box/COVID19_Project
regions=("Genome" "S" "E" "M" 'N' 'ORF1a' 
  'ORF1b' 'ORF1ab' 'ORF3a' 
  'ORF6' 'ORF7a' 'ORF7b' 
  'ORF8' 'ORF10' "2'-O-ribose methyltransferase" 
  "3C-like proteinase" "3'-to-5' exonuclease" 
  "endoRNAse" "helicase" "leader protein" 
  "nsp2" "nsp3" "nsp4" "nsp6" 
  "nsp7" "nsp8" "nsp9" "nsp10" 
  "nsp11" "RNA-dependent RNA polymerase" 
  ) 
dist_method="K80"
for region in "${regions[@]}"
do
  echo ${region}
	# omeClust -i "${home_dir}/data/Distances_500_GISAID_With_NCBI/Distances/${dist_method}/${region}.tsv" -o "/Users/rah/Documents/m2clust_TN93/${region}"  --metadata "${home_dir}/data/Distances_500_GISAID_With_NCBI/Distances/Metadata.txt" --resolution low --linkage complete
	omeClust -i "/Users/rah/Documents/Distances/${dist_method}/${region}.tsv" -o "/Users/rah/Documents/omeClust_${dist_method}/${region}"  --metadata "${home_dir}/data/Distances_500_GISAID_With_NCBI/Distances/Metadata.txt" --resolution low --linkage complete
done


## failed for TN93: E ORF3a ORF6 ORF7a ORF7b ORF8 ORF10 "leader protein" nsp2 nsp11 "RNA-dependent RNA polymerase" 

# failed for K80: ORF3a ORF6 ORF7a ORF7b ORF8 ORF10 nsp11 "RNA-dependent RNA polymerase""

# omeClust -i /Users/rah/Documents/Distances/TN93/Genome.tsv -o /Users/rah/Documents/omeClust_TN93/Genome_low  --metadata ~/Box/COVID19_Project/data/Distances_500_GISAID_With_NCBI/Distances/Metadata.txt --resolution low --linkage complete
# The number of major clusters:  4
# Organism  is the most influential metadata in clusters
# There are 7 clusters
# Output is written in /Users/rah/Documents/omeClust_TN93/Genome_low
# omeClust -i /Users/rah/Documents/Distances/TN93/Genome.tsv -o /Users/rah/Documents/omeClust_TN93/Genome_medium  --metadata ~/Box/COVID19_Project/data/Distances_500_GISAID_With_NCBI/Distances/Metadata.txt --resolution medium --linkage complete
# The number of major clusters:  4
# Organism  is the most influential metadata in clusters
# There are 7 clusters
# Output is written in /Users/rah/Documents/omeClust_TN93/Genome_medium

# omeClust -i /Users/rah/Documents/Distances/TN93/Genome.tsv -o /Users/rah/Documents/omeClust_TN93/Genome_high  --metadata ~/Box/COVID19_Project/data/Distances_500_GISAID_With_NCBI/Distances/Metadata.txt --resolution high --linkage complete
# The number of major clusters:  4
# Organism  is the most influential metadata in clusters
# There are 19 clusters
# Output is written in /Users/rah/Documents/omeClust_TN93/Genome_high



