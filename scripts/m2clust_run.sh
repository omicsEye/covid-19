
regions=("E" "M" "N" "ORF1a" "ORF1b" "ORF1ab" "ORF3a" "ORF6" "ORF7a" "ORF7b" "ORF8" "ORF10") 
dist_method="TN93"
for region in "${regions[@]}"
do
	m2clust -i ~/Box/COVID19_Project/data/Distances_0508/${dist_method}/${region}_500.tsv -o ~/Box/COVID19_Project/analysis/m2clust_0508/${dist_method}/${region}_500 --size-to-plot 100 --metadata ~/Box/COVID19_Project/data/MetaData/metadata.tsv
done
echo "Done running m2clust"
