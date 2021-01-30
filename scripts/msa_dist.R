
# This program calaucltes distance matrix from MSA files for different regions of CoV genome
library(ape)
library(vegan)
library(stringr)

setwd("~/Box/COVID19_Project")

regions <- c("Genome", "S","E","M","N","ORF1a",
             "ORF1b","ORF1ab","ORF3a",
             "ORF6","ORF7a","ORF7b",
             "ORF8","ORF10","2'-O-ribose methyltransferase",
             "3C-like proteinase","3'-to-5' exonuclease",
             "endoRNAse","helicase","leader protein",
             "nsp2","nsp3","nsp4","nsp6",
             "nsp7","nsp8","nsp9","nsp10",
             "nsp11","RNA-dependent RNA polymerase") #'reference'
#num_samples <- 500
similarity_method <- 'TN93' # 'TN93' #   'TN93' # 

# Creat a dirctory for omeCLust input files
output_path <- paste('/Users/rah/Documents/Distances/',similarity_method, sep = '')
#output_path <- paste('/Users/rah/Documents/Distances_2000_GISAID/',similarity_method, sep = '')
dir.create(file.path(output_path), showWarnings = TRUE)
region <- 'Genome'
for (region in regions) {
  print(region)
  seq <- read.FASTA(paste('data/Distances_500_GISAID_With_NCBI/',region, '.fasta', sep=''))
  #seq <- read.FASTA(paste('data/Distances_2000_GISAID/',region, '.fasta', sep=''))
  seq <- read.FASTA('/Users/rah/Box/COVID19_Project/data/Genes_GISAID_15000/msa_0508_reference.fasta')
  seq <- seq[2:length(seq)]
  #sub_seq <- seq[sample(1:length(seq), num_samples, replace=F)]
  #seq['EPI_ISL_463071 ORF10']
  #seq <- seq[!names(seq) %in% c(paste("NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1 complete genome ", 
  #                                    region, sep = ""))]
  seq <- seq[!names(seq) %in% c(paste("NC_045512 ", 
                                      region, sep = ""))]
  seq <- seq[!names(seq) %in% c("NC_045512")]
  
  D <- dist.dna(seq, model = similarity_method, gamma = T, variance = TRUE,
                pairwise.deletion = TRUE,
                base.freq = NULL, as.matrix = TRUE)
  # remove columns with all NA
  D <- as.data.frame(D[,colSums(!is.na(D)) > 0])
  
  # remove rows with all NA
  D <- as.data.frame(D[apply(D, 1, function(y) !all(is.na(y))),])
  
  # remove columns and rowa with less tahn 10 values NA
  D <- as.data.frame(D[,colSums(!is.na(D)) > dim(D)[1]*.1])
  D<- D[colnames(D),]

  # remove all columns and rows with any NA
  #D <- as.data.frame(D[,colSums(is.na(D)) == 0])
  #D<- D[colnames(D),]

  

  # Data with metadata
  #metadata_names <- sapply( strsplit(names(D), "\\|"), "[", 2 )
  #D <- D[!is.na(metadata_names), ]
  #D <- D[, !is.na(metadata_names)]
  #colnames(D) <- rownames(D) <- metadata_names[!is.na(metadata_names)]
  rownames(D) <- gsub(paste(" ",region, sep=""), "",rownames(D))
  colnames(D) <- gsub(paste(" ",region, sep=""), "",colnames(D))
  print(dim(D))
  class(D)
  isSymmetric(as.matrix(D),check.attributes=TRUE)
  #D <- as.matrix(vegdist(D, method="bray"))
  #pheatmap::pheatmap(D)
  write.table( D, paste(output_path, '/', region,'.tsv', sep='') , sep = "\t", eol = "\n", na = "", col.names = NA, quote= F, row.names = T)
}
