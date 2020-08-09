 #### NCBI metadata extraction and cleaning and combining with GISAID ----
# This script takes in 4 fasta files that are downloaded from NCBI (bat, MERS, SARS, OUTgroup). It changes their name to be a unique name from NCBI, and extracts the metadata.  It also takes in GISAID data and metadata and concatenates the two together
# Rebecca Clement rebeccaclement@gwu.edu
# Aug 4, 2020

#load libraries
library(phangorn)
library(ape)
library(stringr)
#install.packages('rentrez')
library("rentrez")

setwd("~/Box/COVID19_Project")

#### read in complete metadata from GISAID ----
metadata<-read.delim("data/MetaData/complete_metadata_0804.tsv")
dim(metadata) #22759 entries, 12 columns

#add a column that says which clade it is
metadata$clade <- 'covid19'

 #### read in files from ncbi ----
#read in bat files downloaded from https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SLen_i=29000%20TO%2040000&Completeness_s=complete&VirusLineage_ss=Bat%20SARS%20coronavirus%20HKU3,%20taxid:442736&VirusLineage_ss=Bat%20SARS-like%20coronavirus,%20taxid:1508227&VirusLineage_ss=Bat%20coronavirus,%20taxid:1508220
batseq <- read.FASTA('data/NCBI/bat-sequences.fasta')
length(batseq) #20
names(batseq)<-str_split_fixed(names(batseq), " ",2)[,1]
#read in SARS files 
# These are now included in GISAID data
#covseq<-read.FASTA('data/NCBI/sars-cov-2_apr21.fasta') #942 sequences
#names(covseq)<-str_split_fixed(names(covseq), " ",2)[,1]
#read in COV19 files
sarseq<-read.FASTA('data/NCBI/sars-related-coronavirus_apr21.fasta') #1254
names(sarseq)<-str_split_fixed(names(sarseq), " ",2)[,1]
#read in mers files
merseq<-read.FASTA('data/NCBI/mers-cov-apr21.fasta') #528 sequences
names(merseq)<-str_split_fixed(names(merseq), " ",2)[,1]
#read in outgroup files
outseq<-read.FASTA('data/NCBI/outgoup.fasta') #43
names(outseq)<-str_split_fixed(names(outseq), "\\|",2)[,2]
names(outseq)<-str_split_fixed(names(outseq), "\\.",2)[,1]

# Combine all NCBI sequences
allNCBI<-c(batseq,merseq,sarseq,outseq)
write.dna(batseq, "data/NCBI/batseq.fa", format="fasta")
write.dna(merseq, "data/NCBI/merseq.fa", format="fasta")
write.dna(outseq, "data/NCBI/outseq.fa", format="fasta")
write.dna(sarseq, "data/NCBI/sarseq.fa", format="fasta")
write.dna(allNCBI, "data/NCBI/NCBIseq.fa", format="fasta")

#upload GISAID Seq, get a subset of it, and then write DNA file
covseq<-read.FASTA('data/gisaid_0804/msa_0804_complete.fasta') 
GISAID_seq_sub <- covseq[sample(1:length(covseq), 2000, replace=F)]
NCBI_GISAID<-c(batseq,merseq,sarseq,outseq,GISAID_seq_sub)
write.dna(NCBI_GISAID, "data/NCBI/NCBI_gisaid0804_2000.fa", format="fasta")

#get metadata of the gisaid_sub #Come back here and fix this Aug9
#metadata_sub<-subset(metadata,names(GISAID_seq_sub)%in%metadata$X)
write.csv(as.data.frame(metadata_sub), file="data/NCBI/metadata_sub0804.csv")

#make a metadata datafram that has all of the names
ncbinames<-c(names(batseq),names(merseq),names(sarseq),names(outseq)) #there are 1845 total here
#get rid of duplicates
ncbinames<-ncbinames[!duplicated(ncbinames)] #now its length is 1844
#### Entrez stuff ----
sum1<- entrez_summary(db="nucleotide",id=ncbinames[1:350]) #looks like 361 is the max
sum2<- entrez_summary(db="nucleotide",id=ncbinames[351:700])
sum3<- entrez_summary(db="nucleotide",id=ncbinames[701:1050])
sum4<- entrez_summary(db="nucleotide",id=ncbinames[1051:1400])
sum5<- entrez_summary(db="nucleotide",id=ncbinames[1401:1750])
sum6<- entrez_summary(db="nucleotide",id=ncbinames[1751:1844])
str(sum1) #to see list of metadata entries

metalist<- extract_from_esummary(sum1, c("caption","uid","taxid","subtype","subname","title","createdate","updatedate","extra","organism","strain","biosample"))
metalist.df<-as.data.frame(t(metalist))
#If you want to take one or more element from every summary record to a list
#for (onething in c(sum2,sum3,sum4,sum5)){
  metalist2<- as.data.frame(t(extract_from_esummary(sum2, c("caption","uid","taxid","subtype","subname","title","createdate","updatedate","extra","organism","strain","biosample"))))
  metalist3<- as.data.frame(t(extract_from_esummary(sum3, c("caption","uid","taxid","subtype","subname","title","createdate","updatedate","extra","organism","strain","biosample"))))
  metalist4<- as.data.frame(t(extract_from_esummary(sum4, c("caption","uid","taxid","subtype","subname","title","createdate","updatedate","extra","organism","strain","biosample"))))
  metalist5<- as.data.frame(t(extract_from_esummary(sum5, c("caption","uid","taxid","subtype","subname","title","createdate","updatedate","extra","organism","strain","biosample"))))
  metalist6<- as.data.frame(t(extract_from_esummary(sum6, c("caption","uid","taxid","subtype","subname","title","createdate","updatedate","extra","organism","strain","biosample"))))
  newmetalist<-rbind(metalist.df,metalist2,metalist3,metalist4,metalist5,metalist6)
#}
View(newmetalist)
df <- apply(newmetalist,2,as.character)
write.table(df, "data/MetaData/metadatacrapNCBI.txt" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
df2<-as.data.frame(df)

#### tyson's magic-parsing out weird metadata ----
#when I ran this, line 920 has the field "note" twice, so I am going to change the second instance to note2
#replace all instances of |note| with |note1|
df2$subtype <- str_replace_all(df2$subtype, "\\|note\\|", "\\|note1\\|")

allFields=c()
for (x in unique(df2$subtype)){
  myFields=strsplit(x, "\\|")
  allFields=append(allFields,myFields[[1]])
}
uniqueFields=unique(allFields)
uniqueFields=uniqueFields[-which(uniqueFields=="strain")]
loopTable=df2
for(x in uniqueFields){
  print(x)
  name=x
  loopTable=cbind(loopTable, x)
}
colnames(loopTable)=c(colnames(df2),
                      uniqueFields)
for(row in 1:nrow(loopTable)){
  thisRowsFields=strsplit(loopTable[row,"subtype"], "\\|")[[1]]
  thisRowsValues=strsplit(loopTable[row,"subname"], "\\|")[[1]]
  for(field in uniqueFields){
    if(field %in% thisRowsFields){
      loopTable[row, field]=thisRowsValues[which(thisRowsFields==field)]
    }
    else{
      loopTable[row, field]=NA
    }
  }
}
finalTable=loopTable
write.csv(loopTable, "data/MetaData/NCBI_Metadata_parsed.csv")
loopTable<-read.csv("data/MetaData/NCBI_Metadata_parsed.csv")

#### formatting metadata ----
#now we want to subset this to only the data we have for GISAID ("country", "country_exposure", "sex", "age","date","virus","host","region")
ncbimeta<-loopTable[,c("caption","organism","host","country","collection_date")]
#change organism to be categorical
ncbimeta$clade <- ifelse(ncbimeta$organism=="Severe acute respiratory syndrome coronavirus 2",'covid19', 
                         ifelse(ncbimeta$organism=="Middle East respiratory syndrome-related coronavirus","MERS",
                                ifelse(ncbimeta$organism=="Severe acute respiratory syndrome-related coronavirus","SARS_related",
                                       ifelse(grepl("Bat", ncbimeta$organism, fixed = TRUE),'BAT',
                                              ifelse(grepl("SARS", ncbimeta$organism, fixed = TRUE),"SARS_related","outgroup" )))))
                       
#read in NCBI metadata
#ncbimeta<-read.csv('data/MetaData/NCBI_Metadata_parsed.csv')

#get country to be better
ncbimeta$country<-str_split_fixed(ncbimeta$country, ":",2)[,1]

#fix dates
#change dates to months
ncbimeta$month<-str_split_fixed(ncbimeta$collection_date, "\\-",3)[,2]
ncbimeta$month<-month.abb[as.numeric(ncbimeta$month)]

#make row names caption
rownames(ncbimeta)<-ncbimeta$caption
#add age, sex, region columns
ncbimeta$age<-NA
ncbimeta$sex<-NA
ncbimeta$region<-NA
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

#### Combine data frames ----
ncbigisaid<-rbind(ncbimeta[,c("clade","country","month","age","sex","region")],metadata2[,c("clade","country","month","age","sex","region")])

#export
write.table(ncbigisaid,'data/MetaData/metadata_ncbi_gisaid_0613.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

# write the file
write.table( metadata2,'data/MetaData/metadata_gisaid_0613.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
