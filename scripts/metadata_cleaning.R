##### metadata cleaning ----

#This script takes in the raw metadata file from GISAID, formats into a standard format, and then outputs 
# 3 files 1) formatted metadata, 2) subset of metadata that has complete metadata, 3) subset of metadata that has a portion of each clade
# It also reads in fasta from GISAID, renames it to them to EPI IDs and creates corresponding files for metadata
# first download raw reads from gisaid
# log in to https://www.epicov.org/epi3/frontend#49f1a9 after you have applied for an account
# Under the Epicov tab, click "Downloads" and click on the most recent msa file as well as "nextmeta" to download the metadata


setwd("~/Box/COVID19_Project")
# read in metadata
metadata <- read.delim(
  'data/gisaid_0804/metadata_2020-08-04_11-57.tsv',
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

dim(metadata) #76047     9


#### Read in GISAID fasta file ----
#### read in the msa ----
seq <- read.FASTA('data/gisaid_0804/msa_0804.fasta')
length(seq) #69613
names(seq)<-str_split_fixed(names(seq), "\\|",3)[,2]
seq_complete<-seq[names(seq)!=""] #length is still 69613

#Get only the metadata that goes with the formatted fasta
metadat_genome<-metadata[rownames(metadata)%in%names(seq_complete),] #69606 entries
# keep only the formatted fasta files that are in the metadata
names(seq_complete[names(seq_complete)%in%rownames(metadat_genome)==FALSE])
seq_complete<-seq_complete[names(seq_complete)%in%rownames(metadat_genome)==TRUE] #Now this has 69606 entries

# filter to rows with complete metadata
comp_meta <- metadat_genome[complete.cases(metadat_genome), ] #length 19819

#make a subset of this that has sequences from each type of pangolin thing (using Tyson's python script)
#read in list of sequences: 
clade400 <- read.delim(
  'data/MetaData/RandomSelection_400.tsv',
  sep = '\t',header = TRUE,fill = T,comment.char = "" ,check.names = F
)
metaclade<-metadat_genome[clade400$gisaid_epi_isl,] #dim 510

#### Write files ----
# write the metadata files. Change date to date of original tsv file
write.table(comp_meta,'data/MetaData/complete_metadata_0804.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
write.table(metadat_genome,'data/MetaData/metadata_interest_0804.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
write.table(metaclade,'data/MetaData/metadata_clades_0804.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

seq_comp<-seq_complete[names(seq_complete)%in%rownames(comp_meta)==TRUE] #19819
seq_clade<-seq_complete[names(seq_complete)%in%rownames(metaclade)==TRUE] #509

# write out formatted fasta file
write.FASTA(seq_complete,"data/gisaid_0804/msa_0804_all.fasta")
write.FASTA(seq_complete,"data/gisaid_0804/msa_0804_complete.fasta")
write.FASTA(seq_complete,"data/gisaid_0804/msa_0804_clades.fasta")
