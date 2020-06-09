
home_dir=~/Box/COVID19_Project
regions=("S" "E" "M" "N" "ORF1a" "ORF1b" "ORF1ab" "ORF3a" "ORF6" "ORF7a" "ORF7b" "ORF8" "ORF10") 
dist_method="TN93"
for region in "${regions[@]}"
do
	m2clust -i ${home_dir}/data/Distances_0508/${dist_method}/${region}_500.tsv -o ${home_dir}/analysis/m2clust_0508/${dist_method}/${region}_500 --metadata ${home_dir}/data/MetaData/metadata.tsv
done
echo "Done running m2clust"
