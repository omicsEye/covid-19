# This program calculates distance matrix from MSA files for different regions of CoV genome
library(ape)
library(vegan)
library(stringr)

setwd("C:/Users/tyson/Box/COVID19_Project/data/Distances_Hackensack")

regions <- c("S", "E", "M", 'N', 'ORF1a',
             'ORF1b', 'ORF1ab', 'ORF3a',
             'ORF6', 'ORF7a', 'ORF7b',
             'ORF8', 'ORF10', "2'-O-ribose methyltransferase",
             "3C-like proteinase", "3'-to-5' exonuclease",
             "endoRNAse", "helicase", "leader protein",
             "nsp2", "nsp3", "nsp4", "nsp6",
             "nsp7", "nsp8", "nsp9", "nsp10",
             "nsp11", "RNA-dependent RNA polymerase") #'reference'
num_samples <- 74
similarity_method <- 'TN93' #'K80'

# Creat a dirctory for MaAsLin input files
output_path <- paste('Distances/',similarity_method, sep = '')
#dir.create(file.path(output_path), showWarnings = TRUE)
#region <- 'S'
for (region in regions) {
  seq <- read.FASTA(file=paste(region, '.fasta', sep=''))
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
  metadata_names <- names(D)
  D <- D[!is.na(metadata_names), ]
  D <- D[, !is.na(metadata_names)]
  colnames(D) <- rownames(D) <- metadata_names[!is.na(metadata_names)]
  #D <- as.matrix(vegdist(D, method="bray"))
  #pheatmap::pheatmap(D)
  write.table( D, paste(output_path, '/', region, '_',num_samples,'.tsv', sep='') , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
}