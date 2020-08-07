#####metadata cleaning

#This script takes in the raw metadata file from GISAID, formats into a standard format, and then outputs 
# 3 files 1) formatted metadata, 2) subset of metadata that has complete metadata, 3) subset of metadata that has a portion of each clade
# first download raw reads from gisaid
# log in to https://www.epicov.org/epi3/frontend#49f1a9 after you have applied for an account
# Under the Epicov tab, click "Downloads" and click on the most recent msa file as well as "nextmeta" to download the metadata


setwd("~/Box/COVID19_Project")
# read in metadata
metadata <- read.delim(
  'data/MetaData/gisaid_0804/metadata_2020-08-04_11-57.tsv',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
# set rownames
rownames(metadata)<-metadata$gisaid_epi_isl

# select the metadata of interest
metadata <- metadata[, c("country", "country_exposure", "sex", "age","date","region","pangolin_lineage")]

# replace charecters to be properly NAs
metadata[metadata == '?'] <- NA
metadata[metadata == 'nan'] <- NA

#replace sex to be male, female or NA
table(metadata$sex)
metadata$sex<- ifelse(metadata$sex=='Female' | metadata$sex=='FEmale' | metadata$sex=='Woman','Female',
                       ifelse(metadata$sex=='Male','Male',NA))
table(metadata$sex)

# fix dates
metadata["EPI_ISL_436426", "date"] <- "2020-04-05"
metadata$date[metadata$date=="2020"] <- "2020-01-01"
metadata$date[metadata$date=="2020-01"] <- "2020-01-01"
metadata$date[metadata$date=="2020-02"] <- "2020-02-01"
metadata$date[metadata$date=="2020-03"] <- "2020-03-01"
metadata$date[metadata$date=="2020-04"] <- "2020-04-01"
metadata$date[metadata$date=="2020-05"] <- "2020-05-01"
metadata$date[metadata$date=="2020-06"] <- "2020-06-01"
metadata$date[metadata$date=="2020-07"] <- "2020-07-01"

#set a start date for COVID-19 pandamic
start <- "2019-12-01" 
Days_from_start <- as.Date(metadata$date)-as.Date(start)
metadata$Days <- Days_from_start

# categorize age as m2clust only accepts categorical data 
metadata$age <- cut(as.numeric(metadata$age), breaks = seq(0,100,by=15), right = TRUE)


# add a column with number of availble metadata 
metadata$n <- rowSums(!is.na(metadata))

# filter to rows with complete metadata
comp_meta <- metadata[complete.cases(metadata), ]

#make a subset of this that has sequences from each type of pangolin thing (using Tyson's python script)
#read in list of sequences: 
clade400 <- read.delim(
  'data/MetaData/RandomSelection_400.tsv',
  sep = '\t',header = TRUE,fill = T,comment.char = "" ,check.names = F
)
metaclade<-metadata[clade400$gisaid_epi_isl,]

# write the file. Change date to date of original tsv file
write.table(comp_meta,'data/MetaData/complete_metadata_0804.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
write.table(metadata,'data/MetaData/metadata_interest_0804.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
write.table(metaclade,'data/MetaData/metadata_clades_0804.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)




