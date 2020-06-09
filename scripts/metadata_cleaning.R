setwd("~/Box/COVID19_Project")

regions <- c("S", "E", "M", 'N', 'ORF1a', 'ORF1b', 'ORF1ab', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF10') #'reference'

metadata <- read.delim(
  'data/MetaData/genome_metadata_formatted.tsv',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
metadata2 <- metadata[, c("country", "country_exposure", "sex", "age")]
metadata2[metadata2 == '?'] <- NA
metadata2[metadata2 == 'nan'] <- NA
metadata2$age <- cut(as.numeric(metadata2$age), breaks = seq(0,100,by=15), right = TRUE)
write.table( metadata2,'data/MetaData/metadata.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)