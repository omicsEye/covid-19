#####metadata cleaning


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

# fix dates
metadata["EPI_ISL_436426", "date"] <- "2020-04-05"
metadata$date[metadata$date=="2020"] <- "2020-01-01"
metadata$date[metadata$date=="2020-01"] <- "2020-01-01"

#set a start date for COVID-19 pandamic
start <- "2019-12-01" 
Days_from_start <- as.Date(metadata$date)-as.Date(start)
metadata$Days <- Days_from_start

# select the metadata of interest
metadata2 <- metadata[, c("country", "country_exposure", "sex", "age", "Days")]

# replace charecters to be properly NAs
metadata2[metadata2 == '?'] <- NA
metadata2[metadata2 == 'nan'] <- NA

# categorize age as m2clust only accepts categorical data 
#metadata2$age <- cut(as.numeric(metadata2$age), breaks = seq(0,100,by=15), right = TRUE)


# add a column with number of availble metadata 
metadata2$n <- rowSums(!is.na(metadata2))

# filter to rows with complete metadata
comp_meta <- metadata2[complete.cases(metadata2), ]


# write the file
write.table(comp_meta,'data/MetaData/complete_metadata.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
write.table(metadata2,'data/MetaData/metadata_interest.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)





