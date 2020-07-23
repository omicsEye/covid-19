# This script takes 3 msa files (in which the ref SARS-CoV-2 genome appears) as
# input. The three files represent the three possible reading frames. It also
# requires the amino acid sequences for the SARS-CoV-2 proteins of interest.
# It then generates and writes a csv-formateed table that contains 
# the nucleotide sequences across all of the genomes for each of the
# defined regions of the SARS-CoV-2 genome. In also
# generates and writes a multifasta file with the same.
# In this way, homologous regions 
# of coronavirus genomes can be fetched and formatted for further downstream
# analysis.

#Load necessary libraries
library(ggplot2)
library(stringr)
library(Biostrings)
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins")
#Your alignment files go here; the naming scheme should be the same
filepaths=c("COVID_genome_translation_rf1 alignment.fasta",
            "COVID_genome_translation_rf2 alignment.fasta",
            "COVID_genome_translation_rf3 alignment.fasta")

#Create objects for each genome sequence across the the three alignments
#First initialize some variables
nameList=c()
rfList=c()
index=1
for(x in filepaths){
  splitOne=str_split(x," ")[[1]][1] #Remove the "alignment.fasta" at the end
  thisRF=str_split(splitOne, "_")[[1]][[4]] #Fetch which reading frame
  con = file(x, "r") 
  while ( TRUE ) { #Iterate through file
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(str_detect(line, "^>")){ #If header then add name to headerList
      nameList=append(nameList,paste(str_split(line, " ")[[1]][1], thisRF, sep="_"))
    }else{
      rfList=append(rfList, line) #If sequence then add to rfList
      index=index+1
    }
  }
  close(con)
}
#Take sequence names and sequences and bind into dataframe 
rfTable=as.data.frame(cbind(nameList,rfList))

#Next we iterate through directory of amino acid sequences for 25 proteins
#First, initialize some variables
tableExists=FALSE
listOfRows=c()
loopCompleted=FALSE
listOfFiles=list.files("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/Reference_Protein_Sequences", full.names = TRUE)
for(x in listOfFiles){
  if(str_detect(x, "fasta")){ #Make sure it's a .fasta file
    con = file(x, "r")
    while ( TRUE ) {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }else if(str_detect(line, "^>")){ #If it's a header, get the name
        name=line
      }else if (!str_detect(line, "^>") && !(length(line)==0)){
        thisSequence=line #If it's a sequence then save it and align later
      }
    }
    close(con)
    #See which rf gets the best alignment by aligning to each of the 3 rfs
    #We align to the SARS-CoV-2 genome with low gap penalty
    localAlign1 <-
      pairwiseAlignment(AAString(thisSequence),
                        AAString(rfTable$rfList[which(rfTable$nameList==">NC_045512_SARSCoV2_translation_rf1")]),
                        gapOpening=-0.1, 
                        gapExtension=-0.01,
                        scoreOnly=FALSE,
                        type = "local")
    localAlign2 <-
      pairwiseAlignment(AAString(thisSequence),
                        AAString(rfTable$rfList[which(rfTable$nameList==">NC_045512_SARSCoV2_translation_rf2")]),
                        gapOpening=-0.1, 
                        gapExtension=-0.01,
                        scoreOnly=FALSE,
                        type = "local")
    localAlign3 <-
      pairwiseAlignment(AAString(thisSequence),
                        AAString(rfTable$rfList[which(rfTable$nameList==">NC_045512_SARSCoV2_translation_rf3")]),
                        gapOpening=-0.1, 
                        gapExtension=-0.01,
                        scoreOnly=FALSE,
                        type = "local")
    #Calculate the maximum score
    theMax=max(localAlign1@score,localAlign2@score, localAlign3@score)
    #Whichever alignment gave the max score becomes our alignment
    if(localAlign1@score==theMax){
      localAlign=localAlign1
      bestAlignName=">NC_045512_SARSCoV2_translation_rf1"
    }
    if(localAlign2@score==theMax){
      localAlign=localAlign2
      bestAlignName=">NC_045512_SARSCoV2_translation_rf2"
    }
    if(localAlign3@score==theMax){
      localAlign=localAlign3
      bestAlignName=">NC_045512_SARSCoV2_translation_rf3"
    }
    #Error check whether our sequence from the alignment is the same as the
    #query but with gaps
    if(thisSequence!=str_remove_all(substring(rfTable$rfList[which(rfTable$nameList==bestAlignName)],
                                   localAlign@subject@range@start,
                                   localAlign@subject@range@start + localAlign@subject@range@width-1), "-")){
      #If there is an error, write the two sequences to a file for the user to view
      fileName=paste(str_remove_all(str_split(x, "/")[[1]][length(str_split(x, "/")[[1]])], ".fasta"), "Error.txt", sep="_")
      print(paste("There was an error with the alignment step for the locus", x, "Error report written to", fileName, sep=" "))
      cat("Reference sequence", file=fileName)
      cat("\n", file=fileName, append = TRUE)
      cat(thisSequence, file=fileName, append = TRUE)
      cat("\n", file=fileName, append = TRUE)
      cat("\n", file=fileName, append = TRUE)
      cat("Alignment sequence", file=fileName, append = TRUE)
      cat("\n", file=fileName, append = TRUE)
      cat(str_remove_all(substring(rfTable$rfList[which(rfTable$nameList==bestAlignName)],
                                   localAlign@subject@range@start,
                                   localAlign@subject@range@start + localAlign@subject@range@width-1), "-"), file=fileName, append = TRUE)
    }
    #Now that we know the locus of the alignment we fetch those regions
    #from all genomes
    listOfGenomes=c()
    for(q in rfTable$nameList){ #Iterate through all genomes across the 3 rfs
      nameOfProtein=str_remove_all(str_split(x, "/")[[1]][length(str_split(x, "/")[[1]])], ".fasta")
      listOfGenomes=append(listOfGenomes, substring(rfTable$rfList[which(rfTable$nameList==q)],
                                                    localAlign@subject@range@start,
                                                    localAlign@subject@range@start + localAlign@subject@range@width-1))
      if(!loopCompleted){#First time through the list save the names
        listOfRows=append(listOfRows, q)
      }
    }
    loopCompleted=TRUE
    if(tableExists){ #If table exists, add to it
      myTable=as.data.frame(cbind(myTable,listOfGenomes))
      colnames(myTable)[length(colnames(myTable))]=nameOfProtein
    }else{ #If table doesn't exist, make it
      tableExists=TRUE
      myTable=as.data.frame(listOfGenomes)
      colnames(myTable)[length(colnames(myTable))]=nameOfProtein
    }
  }
}
row.names(myTable)=listOfRows #Add row names to table

finalTable=myTable
f=0
#Remove sequences with stop codons
for(w in 1:length(colnames(finalTable))){
  for (y in 1:length(row.names(finalTable))){
    if(str_detect(finalTable[y,w], "\\*")){ #If there is a stop codon
      finalTable[y,w]=" "
    }else{
      f=f+1
    }
  }
}

#Reorder the proteins so they are grouped together
myReorder=c()
for(i in 1:18){
  myReorder=append(myReorder, c(i, i+18, i+36))
}
finalTable=finalTable[myReorder,]
#Now we format the dataframe and write it to a .csv file
#First we initialize some variables
constructRow=c()
constructCol=c()
tableExists=FALSE
addRow=FALSE
p=1
rowInd=1
while(p < (length(row.names(finalTable)))){ #Iterate through all rows and cols
  for(l in 1:length(colnames(finalTable))){
    thisUnit=finalTable[c(p, p+1,p+2), l] #Look at all three reading frames
    if(length(which(thisUnit!=" "))>0){
      #If two equally good sequence alignments were found, add the one with 
      #fewer gaps
      if(length(which(thisUnit!=" "))>1){
        firstCollapsed=(thisUnit[which(thisUnit!=" ")])[1]
        secondCollapsed=(thisUnit[which(thisUnit!=" ")])[2]
        firstGap=str_count(firstCollapsed, "-")
        secondGap=str_count(secondCollapsed, "-")
        gapMin=min(firstGap,secondGap)
        if(firstGap==gapMin){
          collapsed=(thisUnit[which(thisUnit!=" ")])[1]
        }
        if(secondGap==gapMin){
          collapsed=(thisUnit[which(thisUnit!=" ")])[2]
        }
        constructRow=append(constructRow, collapsed)
      }else{ #If there is only one sequence, add it
        collapsed=thisUnit[which(thisUnit!=" ")]
        constructRow=append(constructRow, collapsed)
      }
    }else{ #If there is no sequence, leave it blank
      collapsed=" "
      constructRow=append(constructRow, collapsed)
    }
  }
  if(!tableExists){ #if table doesn't exist, make it
    collapsedTable=as.data.frame(constructRow)
    tableExists=TRUE
  }else{ #if table does exist, add to it
    collapsedTable=as.data.frame(cbind(collapsedTable, constructRow))
  }
  constructRow=c()
  p=p+3 #Look forward to the next three rows
  rowInd=rowInd+1
}
collapsedTable=t(collapsedTable) #Transpose
colnames(collapsedTable)=colnames(finalTable) #Add column names
collapsedTableRowNames=c()
#Get the column names by looking at the table we made earlier
for(i in 1:(length(row.names(myTable))/3)){
  print(i)
  thisRowName=str_remove_all(row.names(myTable)[i], "_rf1")
  collapsedTableRowNames=append(collapsedTableRowNames,thisRowName)
}
row.names(collapsedTable)=collapsedTableRowNames
write.csv(collapsedTable, file="AllAminoAcidRegionsWithSequences.csv")

#Write results to a multifasta file
for(y in 1:length(row.names(collapsedTable))){
  for(z in 1:length(colnames(collapsedTable))){
    header=paste(row.names(collapsedTable)[y],
                 colnames(collapsedTable)[z], sep = " ")
    mySequence=paste("\n", collapsedTable[y,z], "\n", sep="")
    cat(header, file="18CoronavirusesMSA-AminoAcid.fasta", append = TRUE)
    cat(mySequence, file="18CoronavirusesMSA-AminoAcid.fasta", append = TRUE)
  }
}