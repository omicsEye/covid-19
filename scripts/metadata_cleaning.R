setwd("~/Box/COVID19_Project")

metadata <- read.delim(
  'data/MetaData/genome_metadata_formatted.tsv',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# select the metadata of interest
metadata2 <- metadata[, c("country", "country_exposure", "sex", "age")]

# replace charecters to be properly NAs
metadata2[metadata2 == '?'] <- NA
metadata2[metadata2 == 'nan'] <- NA

# categorize age as m2clust only accepts categorical data 
metadata2$age <- cut(as.numeric(metadata2$age), breaks = seq(0,100,by=15), right = TRUE)

# write the file
write.table( metadata2,'data/MetaData/metadata.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)