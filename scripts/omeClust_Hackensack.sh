
home_dir=~/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/hackensack
regions=(
  "Genome" "S" "E" "M" 'N' 'ORF1a' 
  'ORF1b' 'ORF1ab' 'ORF3a' 
  'ORF6' 'ORF7a' 'ORF7b' 
  'ORF8' 'ORF10' "2'-O-ribose methyltransferase" 
  "3C-like proteinase" "3'-to-5' exonuclease" 
  "endoRNAse" "helicase" "leader protein" 
  "nsp2" "nsp3" "nsp4" "nsp6" 
  "nsp7" "nsp8" "nsp9" "nsp10" 
  "nsp11" "RNA-dependent RNA polymerase") 
dist_method="TN93"
for region in "${regions[@]}"
do
	m2clust -i "${home_dir}/data/${dist_method}/${region}.tsv" -o "${home_dir}/analysis/m2clust/${region}_metadata"  --metadata "${home_dir}/data/HackensackMetaData.txt" --resolution low --linkage complete
done
echo "Done running m2clust"
### --metadata "${home_dir}/data/Targeted_Clinical_Data.txt"