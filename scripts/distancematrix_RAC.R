# making distance matrix for coronavirus data
# Rebecca Clement rebeccaclement@gwu.edu
#25 May 2020 updated June 29

#Ape doesn't include the models of evolution that we want to use so instead we use phangorn
library(phangorn)
library(ape)
library(stringr)

setwd("Documents/Research/Crandall/Coronavirus/Genes_GISAID/")
# read in metadata file
metadat<-read.delim("metadata_2020-05-08.tsv")
dim(metadat)
#### Fix metadata ----

#get an alignment that has the right names
seq <- read.FASTA('msa_0508_reference.fasta')
length(seq) #15930
genome_seq <- seq[2:length(seq)]
length(genome_seq)
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
genome_seq<-genome_seq[names(genome_seq)!=""] #length is now 15745
#Get only the metadata that goes with the formatted fasta
metadat_genome<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq),]
table(metadat$gisaid_epi_isl%in%names(genome_seq))
object<-names(genome_seq)%in%metadat$gisaid_epi_isl
names(genome_seq[object==FALSE])
dim(metadat_genome) #15721
genome_seq<-genome_seq[object==TRUE]

write.FASTA(genome_seq,"msa_0508_genome_formatted.fasta")
#set names so that they are the EPI ones
rownames(metadat_genome)<-metadat_genome$gisaid_epi_isl
#write in m2clust format
write.table(metadat_genome, "genome_metadata_formatted.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)


#from Ali
# E gene ----

# read MSA file
seq <- read.FASTA('msa_0508_E_gene.fasta')
length(seq)
#remove the first alignment (reference)
genome_seq <- seq[2:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""]
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl]
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)]
#calculate the distance 'model' must be one of: "RAW" "JC69" "K80" "F81" "K81" "F84" "T92" "TN93" "GG95" "LOGDET" "BH87" "PARALIN" "N" "TS" "TV" "INDEL" "INDELBLOCK"
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
# write in a format suitable for m2clust 
write.table( D, "TN93/E_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_E_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_E_sub_2000)<-metadat_E_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_E_sub_2000<-metadat_E_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_E_sub_2000, "TN93/E_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#Run the following on terminal: m2clust -i E_gene_2000_dist_TN93.tsv --metadata E_gene_metadata.tsv -o E_gene_out




# ORF1a gene ---- 
#replace all E_ with ORF1a_

# read MSA file
seq <- read.FASTA('msa_0508_ORF1a_gene.fasta')
length(seq) #15929
#remove the first alignment (reference)
genome_seq <- seq[2:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""] #should be 15744 now
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl] #should be 15720
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)] #length is 2000
#calculate the distance 
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
# write in a format suitable for m2clust 
write.table( D, "TN93/ORF1a_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_ORF1a_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_ORF1a_sub_2000)<-metadat_ORF1a_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_ORF1a_sub_2000<-metadat_ORF1a_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_ORF1a_sub_2000, "TN93/ORF1a_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
#run the folloing on terminal: m2clust -i ORF1a_gene_2000_dist_TN93.tsv --metadata ORF1a_gene_metadata.tsv -o ORF1a_gene_out


# ORF1b gene ---- 
#replace all E_ with ORF1b_

# read MSA file
seq <- read.FASTA('msa_0508_ORF1b_gene.fasta')
length(seq) #15929
#remove the first alignment (reference)
genome_seq <- seq[2:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""] #should be 15744 now
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl] #should be 15720
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)] #length is 2000
#calculate the distance 
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
# write in a format suitable for m2clust 
write.table( D, "TN93/ORF1b_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_ORF1b_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_ORF1b_sub_2000)<-metadat_ORF1b_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_ORF1b_sub_2000<-metadat_ORF1b_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_ORF1b_sub_2000, "TN93/ORF1b_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
#run the folloing on terminal: m2clust -i ORF1b_gene_2000_dist_TN93.tsv --metadata ORF1b_gene_metadata.tsv -o ORF1b_gene_out

# S gene ---- 
#replace all E_ with S_

# read MSA file
seq <- read.FASTA('msa_0508_S_gene.fasta')
length(seq) #15929
#remove the first alignment (reference)
genome_seq <- seq[2:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""] #should be 15744 now
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl] #should be 15720
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)] #length is 2000
#calculate the distance 
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
# write in a format suitable for m2clust 
write.table( D, "TN93/S_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_S_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_S_sub_2000)<-metadat_S_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_S_sub_2000<-metadat_S_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_S_sub_2000, "TN93/S_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
#run the folloing on terminal: m2clust -i S_gene_2000_dist_TN93.tsv --metadata S_gene_metadata.tsv -o S_gene_out

# ORF3a_E gene ---- 
#replace all E_ with ORF3a_E_

# read MSA file
seq <- read.FASTA('msa_0508_ORF3a_E_gene.fasta')
length(seq) #15929
#remove the first alignment (reference)
genome_seq <- seq[2:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""] #should be 15744 now
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl] #should be 15720
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)] #length is 2000
#calculate the distance 
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
# write in a format suitable for m2clust 
write.table( D, "TN93/ORF3a_E_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_ORF3a_E_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_ORF3a_E_sub_2000)<-metadat_ORF3a_E_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_ORF3a_E_sub_2000<-metadat_ORF3a_E_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_ORF3a_E_sub_2000, "TN93/ORF3a_E_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
#run the folloing on terminal: m2clust -i ORF3a_E_gene_2000_dist_TN93.tsv --metadata ORF3a_E_gene_metadata.tsv -o ORF3a_E_gene_out


# M gene ---- 
#replace all E_ with M_

# read MSA file
seq <- read.FASTA('msa_0508_M_gene.fasta')
length(seq) #15929
#remove the first alignment (reference)
genome_seq <- seq[1:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""] #should be 15744 now
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl] #should be 15720
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)] #length is 2000
#calculate the distance 
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
# write in a format suitable for m2clust 
write.table( D, "TN93/M_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_M_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_M_sub_2000)<-metadat_M_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_M_sub_2000<-metadat_M_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_M_sub_2000, "TN93/M_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
#run the folloing on terminal: m2clust -i M_gene_2000_dist_TN93.tsv --metadata M_gene_metadata.tsv -o M_gene_out


# ORF6_7a gene ---- 
#replace all E_ with ORF6_7a_

# read MSA file
seq <- read.FASTA('msa_0508_ORF6_7a_gene.fasta')
length(seq) #15929
#remove the first alignment (reference)
genome_seq <- seq[2:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""] #should be 15744 now
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl] #should be 15720
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)] #length is 2000
#calculate the distance 
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
# write in a format suitable for m2clust 
write.table( D, "TN93/ORF6_7a_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_ORF6_7a_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_ORF6_7a_sub_2000)<-metadat_ORF6_7a_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_ORF6_7a_sub_2000<-metadat_ORF6_7a_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_ORF6_7a_sub_2000, "TN93/ORF6_7a_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
#run the folloing on terminal: m2clust -i ORF6_7a_gene_2000_dist_TN93.tsv --metadata ORF6_7a_gene_metadata.tsv -o ORF6_7a_gene_out

# ORF8 gene ---- 
#replace all E_ with ORF8_

# read MSA file
seq <- read.FASTA('msa_0508_ORF8_gene.fasta')
length(seq) #15929
#remove the first alignment (reference)
genome_seq <- seq[2:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""] #should be 15744 now
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl] #should be 15720
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)] #length is 2000
#calculate the distance 
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
# write in a format suitable for m2clust 
write.table( D, "TN93/ORF8_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_ORF8_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_ORF8_sub_2000)<-metadat_ORF8_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_ORF8_sub_2000<-metadat_ORF8_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_ORF8_sub_2000, "TN93/ORF8_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
#run the folloing on terminal: m2clust -i ORF8_gene_2000_dist_TN93.tsv --metadata ORF8_gene_metadata.tsv -o ORF8_gene_out

# ORF10 gene ---- 
#replace all E_ with ORF10_
#There are NAs in this because some sequences are shorter than others

# read MSA file
seq <- read.FASTA('msa_0508_ORF10_gene.fasta')
#I wonder if I can remove al sequences that are less than 15923
length(seq) #15923. fewer bcause it's at the end
#remove the first alignment (reference)
genome_seq <- seq[2:length(seq)]
# change names so it's the part that starts with EPI
names(genome_seq)<-str_split_fixed(names(genome_seq), "\\|",3)[,2]
# remove all alignments that don't have "EPI" in their name"
genome_seq<-genome_seq[names(genome_seq)!=""] #should be 15744 now
#Also remove all alignments that are not in the metadata file
genome_seq<-genome_seq[names(genome_seq)%in%metadat$gisaid_epi_isl] #should be 15720
# get a subset of 2000 random samples from the file
genome_seq_sub <- genome_seq[sample(1:length(genome_seq), 2000, replace=F)] #length is 2000
#calculate the distance 
D <- dist.dna(genome_seq_sub, model = "TN93", gamma = T, variance = TRUE,
              pairwise.deletion = TRUE,
              base.freq = NULL, as.matrix = TRUE)
#replace na's with 0s
#newD<-na.omit(D)
#D[is.na(D)]=0
# I kept getting the error that there were the wrong lengths, so I opened up table in excel and deleted all rows and columns with NA and then it worked

# write in a format suitable for m2clust 
write.table( D, "TN93/ORF10_gene_2000_dist_TN93.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)

#get a subset of metadata that is in the 2000 subset
metadat_ORF10_sub_2000<-metadat[metadat$gisaid_epi_isl%in%names(genome_seq_sub),]
#set names so that they are the EPI ones
rownames(metadat_ORF10_sub_2000)<-metadat_ORF10_sub_2000$gisaid_epi_isl
#subset to only categorical metadata
metadat_ORF10_sub_2000<-metadat_ORF10_sub_2000[c("sex","region","country","division","host","age","country_exposure","region_exposure")]
#write in m2clust format
write.table(metadat_ORF10_sub_2000, "TN93/ORF10_gene_metadata.tsv" , sep = "\t", eol = "\n", col.names = NA, quote= F, row.names = T)
#run the folloing on terminal: m2clust -i ORF10_gene_2000_dist_TN93.tsv --metadata ORF10_gene_metadata.tsv -o ORF10_gene_out

