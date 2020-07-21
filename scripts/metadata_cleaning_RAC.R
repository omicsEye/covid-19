 
#from raw data to metadata and ready for analysis
# Rebecca Clement rebeccaclement@gwu.edu
# June 29, 2020

# first download raw reads from gisaid
# log in to https://www.epicov.org/epi3/frontend#49f1a9 after you have applied for an account
# Under the Epicov tab, click "Downloads" and click on the most recent msa file as well as "nextmeta" to download the metadata

#load libraries
library(phangorn)
library(ape)
library(stringr)
#install.packages('rentrez')
library("rentrez")

setwd("~/Box/COVID19_Project")

#### read in the msa ----
seq <- read.FASTA('data/msa_0613/msa_0613.fasta')
length(seq) #41748
names(seq)<-str_split_fixed(names(seq), "\\|",3)[,2]
seq_complete<-seq[names(seq)!=""] #length is now 15745

#### format metadata file ----
metadat<-read.delim("data/msa_0613/metadata_2020-06-16_11-55.tsv")
dim(metadat) #46474 entries
#Get only the metadata that goes with the formatted fasta
metadat_genome<-metadat[metadat$gisaid_epi_isl%in%names(seq_complete),] #41721 entries
# keep only the formatted fasta files that are in the metadata
names(seq_complete[names(seq_complete)%in%metadat$gisaid_epi_isl==FALSE])
genome_seq<-seq_complete[names(seq_complete)%in%metadat$gisaid_epi_isl==TRUE]
length(genome_seq) #41721
# write out formatted fasta file
write.FASTA(genome_seq,"data/msa_0613/msa_0613_genome_formatted.fasta")

#set names so that they are the EPI ones
rownames(metadat_genome)<-metadat_genome$gisaid_epi_isl

# select the metadata of interest
metadata2 <- metadat_genome[, c("country", "country_exposure", "sex", "age","date","virus","host","region")]

# replace charecters to be properly NAs
metadata2[metadata2 == '?'] <- NA
metadata2[metadata2 == 'nan'] <- NA

#replace sex to be good
metadata2$sex[metadata2$sex == 'FEmale'] <- 'Female'
metadata2$sex[metadata2$sex == 'Woman'] <- 'Female'
metadata2$sex[metadata2$sex == 'unknwon'] <- 'Unknown'
table(metadata2$sex)

#change dates to months
metadata2$month<-str_split_fixed(metadata2$date, "\\-",3)[,2]
metadata2$month<-month.abb[as.numeric(metadata2$month)]


# categorize age as m2clust only accepts categorical data 
metadata2$age <- cut(as.numeric(metadata2$age), breaks = seq(0,100,by=15), right = TRUE)

#make a column that says which clade it is
metadata2$clade <- 'covid19'

# write the file
write.table( metadata2,'data/MetaData/metadata_gisaid_0613.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

 #### read in files from ncbi ----
#read in bat files downloaded from https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SLen_i=29000%20TO%2040000&Completeness_s=complete&VirusLineage_ss=Bat%20SARS%20coronavirus%20HKU3,%20taxid:442736&VirusLineage_ss=Bat%20SARS-like%20coronavirus,%20taxid:1508227&VirusLineage_ss=Bat%20coronavirus,%20taxid:1508220
batseq <- read.FASTA('data/NCBI/bat-sequences.fasta')
length(batseq) #20
names(batseq)<-str_split_fixed(names(batseq), " ",2)[,1]
#read in SARS files 
covseq<-read.FASTA('data/NCBI/sars-cov-2_apr21.fasta') #942 sequences
names(covseq)<-str_split_fixed(names(covseq), " ",2)[,1]
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


#### combine subset of gisaid with these other sequences things ----
#subset gisaid to 1000
gisaid_sub <- genome_seq[sample(1:length(genome_seq), 1000, replace=F)]
allseqs <- cbind(gisaid_sub,batseq,covseq,sarseq,merseq,outseq)
write.dna(gisaid_sub, "~/Documents/Research/Crandall/Coronavirus/NCBI/gisaid.fa", format="fasta")
write.dna(batseq, "~/Documents/Research/Crandall/Coronavirus/NCBI/batseq.fa", format="fasta")
write.dna(covseq, "~/Documents/Research/Crandall/Coronavirus/NCBI/covseq.fa", format="fasta")
write.dna(merseq, "~/Documents/Research/Crandall/Coronavirus/NCBI/merseq.fa", format="fasta")
write.dna(outseq, "~/Documents/Research/Crandall/Coronavirus/NCBI/outseq.fa", format="fasta")
write.dna(sarseq, "~/Documents/Research/Crandall/Coronavirus/NCBI/sarseq.fa", format="fasta")
# in command line, combine these sequences
# cat gisaid.fa batseq.fa covseq.fa merseq.fa outseq.fa sarseq.fa > allseqs.fa 
# When you count, grep -c "^>" allseqs.fa Should add up to 3787. If it doesn't, open up in geneious and then export as fasta
# align using mafft


#get metadata of the gisaid_sub
metadata_sub<-metadata2[names(gisaid_sub),]
write.csv(as.data.frame(metadata_sub), file="~/Documents/Research/Crandall/Coronavirus/NCBI/metadata_sub.csv")

#make a metadata datafram that has all of the names
ncbinames<-c(names(batseq),names(covseq),names(merseq),names(sarseq),names(outseq)) #there are 2787 total here
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

#### Combine data frames ----
ncbigisaid<-rbind(ncbimeta[,c("clade","country","month","age","sex","region")],metadata2[,c("clade","country","month","age","sex","region")])

#export
write.table(ncbigisaid,'data/MetaData/metadata_ncbi_gisaid_0613.tsv' , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
