
# This program calaucltes distance matrix from MSA files for different regions of CoV genome
library(ape)
library(vegan)
library(stringr)

setwd("~/Box/COVID19_Project")

regions <- c("S", "E", "M", 'N', 'ORF1a', 'ORF1b', 'ORF1ab', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF10') #'reference'
num_samples <- 500
similarity_method <- 'TN93' #'K80'

# Creat a dirctory for MaAsLin input files
output_path <- paste('data/Distances_0508/',similarity_method, sep = '')
dir.create(file.path(output_path), showWarnings = TRUE)
#region <- 'S'
for (region in regions) {
  seq <- read.FASTA(paste('data/Genes_GISAID_15000/', 'msa_0508_',region, '_gene.fasta', sep=''))
  #seq <- seq[2:length(seq)]
  sub_seq <- seq[sample(1:length(seq), num_samples, replace=F)]
  D <- dist.dna(sub_seq, model = similarity_method, gamma = F, variance = TRUE,
                pairwise.deletion = TRUE,
                 base.freq = NULL, as.matrix = TRUE)
  # remove columns with all NA
  D <- as.data.frame(D[,colSums(!is.na(D)) > 0])
  
  # remove rows with all NA
  D <- as.data.frame(D[apply(D, 1, function(y) !all(is.na(y))),])
  
  # Data with metadata
  metadata_names <- sapply( strsplit(names(D), "\\|"), "[", 2 )
  D <- D[!is.na(metadata_names), ]
  D <- D[, !is.na(metadata_names)]
  colnames(D) <- rownames(D) <- metadata_names[!is.na(metadata_names)]
  #D <- as.matrix(vegdist(D, method="bray"))
  #pheatmap::pheatmap(D)
  write.table( D, paste(output_path, '/', region, '_',num_samples,'.tsv', sep='') , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
}
