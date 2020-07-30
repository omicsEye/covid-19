# This script takes an msa file (in which the ref SARS-CoV-2 genome appears) as
# input. It then generates and writes a csv-formateed table that contains 
# the nucleotide sequences across all of the genomes for each of the
# defined regions of the SARS-CoV-2 genome. In also
# generates and writes a multifasta file with the same.
# In this way, homologous regions 
# of coronavirus genomes can be fetched and formatted for further downstream
# analysis.

#Load necessary libraries
#BiocManager::install("genbankr") #works with R version 3.6.3
library(genbankr)
library(Biostrings)
library(stringr)
#Change to your working directory
setwd('~/Box/COVID19_Project/data/Nucleotide/')

#Use genbankr methods to load SARS-CoV-2 genbank file
id=GBAccession("NC_045512.2") #Accession number for SARS-CoV-2 Ref
gb=readGenBank(id)

#Read the genome from the genbank file
fullGenome=toString(gb@sequence$`Severe acute respiratory syndrome coronavirus2`)

#SARS-CoV-2 has genes but also has "other features"
#These include "mat_peptides" that code for nsps
#We load both genes and other features and combine into a table
myOtherFeatures=as.data.frame(gb@other_features@elementMetadata@listData)
myOtherFeatures$start=gb@other_features@ranges@start
myOtherFeatures$end=myOtherFeatures$start + gb@other_features@ranges@width - 1
myGenes=as.data.frame(gb@genes@elementMetadata@listData)
myGenes$start=gb@genes@ranges@start
myGenes$end=myGenes$start + gb@genes@ranges@width -1 
#Select columns that we need
myGenes=myGenes[,c(2,12,13)]
myOtherFeatures=myOtherFeatures[,c(5,10,11)]
# I want to add a gene that is just ORF1a and one that is just ORF1b
# first add levels
levels(myGenes$gene)<-c(levels(myGenes$gene),"ORF1a","ORF1b")
myGenes<-rbind(myGenes,c('ORF1a',266,13469))
myGenes<-rbind(myGenes,c('ORF1b',13470,21555))

#Change column names so they match
colnames(myGenes)=c("locus","start","end")
colnames(myOtherFeatures)=c("locus","start","end")
#Remove any duplicated features
myOtherFeatures=myOtherFeatures[which(!duplicated(myOtherFeatures$locus)),]
#Combine genes and other_features into a single "locusTable"
locusTable=as.data.frame(rbind(myGenes,myOtherFeatures))
#Remove any NA values
index=which(is.na(locusTable$locus))
locusTable=locusTable[-index,]

#Find the subsequences in the genome that correspond to our loci
colOfSequences=c()
for(x in 1:length(locusTable$locus)){
  thisSeq=substr(fullGenome,locusTable$start[x],locusTable$end[x])
  colOfSequences=append(colOfSequences, thisSeq)
}
locusTable$Sequences=colOfSequences

#Next we load our alignment and put the sequences into a
#dataframe that can be easily parsed
con = file("../hcov_mers_bat_covid_0421/subsample_test.fasta", "r")
#We will make a vector of sequence names ie MERS, BatCoronavirusAncestor1, etc.
sequenceNames=c()
#We will make a vector of nucleotide sequences
sequenceContents=c()
while ( TRUE ) {
  line = readLines(con, n = 1) #Read lines one at a time
  if ( length(line) == 0 ) {
    break #This might be the end of the file
  }else if(str_detect(line, "^>")){ #If it's a header
    name=line
  }else if (!str_detect(line, "^>") && !(length(line)==0)){ #if it's a sequence
    thisSequence=line
    #Add this header to names vector
    sequenceNames=append(sequenceNames,name)
    #Add this sequence to sequence vector
    sequenceContents=append(sequenceContents,thisSequence)
  }
}
close(con)
#Bind the vectors into a table
alignmentTable=as.data.frame(cbind(sequenceNames, sequenceContents))

#Next we will see what region of the msa corresponds to each of our loci
#across each of our genomes, fetching the subsequences
#and adding them to a table
#Inititialize some variables that we will use in the loop
loopCompleted=FALSE
tableExists=FALSE
listOfRows=c()
#For each of our loci
for(x in locusTable$Sequences){
  #Perform local alignment to the aligned SARS-CoV-2 sequence
  #Use low gap penalty so that we can account for gaps introduced by msa
  localAlign <-
    pairwiseAlignment(x,
                      #alignmentTable$sequenceContents[which(alignmentTable$sequenceNames==">NC_045512_SARSCoV2 |Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1| complete genome")], # I am changing this so that the alignment does not have to include this sequence
                      alignmentTable$sequenceContents[1],
                      gapOpening=-0.1, 
                      gapExtension=-0.01,
                      scoreOnly=FALSE,
                      type = "local")
  listOfGenomes=c()
  #Error check whether our sequence from the alignment is the same as the
  #query but with gaps
  if(x!=str_remove_all(substring(alignmentTable$sequenceContents[1], #Becca changed this to hopefully be more universal
                  localAlign@subject@range@start,
                  localAlign@subject@range@start + localAlign@subject@range@width-1), "-")){
    print(paste("There was an error with the alignment step for the locus",
                locusTable$locus[which(locusTable$Sequences==x)], sep=" "))
  }
  #We have the alignment for SARS-CoV-2 so now we find the sequences
  #At that locus for all of the other genomes
  for(q in 1:length(alignmentTable$sequenceNames)){
    nameOfLocus=alignmentTable$sequenceNames[q]
    sequenceToAdd=
    listOfGenomes=append(listOfGenomes, substring(alignmentTable$sequenceContents[q],
                                                  localAlign@subject@range@start,
                                                  localAlign@subject@range@start + localAlign@subject@range@width-1))
    if(!loopCompleted){ #Make a list of the genome names
      listOfRows=append(listOfRows, str_remove_all(alignmentTable$sequenceNames[q], ","))
    }
  }
  loopCompleted=TRUE
  if(tableExists){ #If the table exist, bind to it
    myTable=as.data.frame(cbind(myTable,listOfGenomes))
    colnames(myTable)[length(colnames(myTable))]=locusTable$locus[which(locusTable$Sequences==x)]
  }else{ #if the table doesn't exist yet, make it
    tableExists=TRUE
    myTable=as.data.frame(listOfGenomes)
    colnames(myTable)[length(colnames(myTable))]=locusTable$locus[which(locusTable$Sequences==x)]
  }
}
#Set the row names to the list of genomes we made earlier (line 105)
row.names(myTable)=listOfRows
# set column names
colnames(myTable)=locusTable$locus
#Write results to a table
write.table(cbind(row.names(myTable),myTable), file="becca/SequenceTable.csv", quote=FALSE, sep=",", row.names = FALSE)

#Write results to a multifasta file
for(y in 1:length(row.names(myTable))){
  for(z in 1:length(colnames(myTable))){
    header=paste(row.names(myTable)[y], colnames(myTable)[z], sep = " ")
    mySequence=paste("\n", myTable[y,z], "\n", sep="")
    myFileName=paste('becca/',colnames(myTable)[z], ".fasta", sep="")
    cat(header, file=myFileName, append = TRUE)
    cat(mySequence, file=myFileName, append = TRUE)
  }
}
