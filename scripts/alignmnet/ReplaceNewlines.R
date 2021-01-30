# This script will take all of the fasta files in a directory
# and remove any newlines in the sequence portion of the files
# Note: Only designed for fasta, not multifasta
library(readtext)
library(stringr)
#Get a vector of all of the fasta files in a chosen directory
listOfFiles=list.files("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/Reference_Protein_Sequences", full.names = TRUE)

# Iterate through the files; first line is header, all subsequent lines get 
# newlines removed so they're on the same line
for(x in listOfFiles){
  if(str_detect(x, "fasta")){ #Make sure it's a fasta file
    myTXT=readtext(x) #Read in the file
    text=str_split(myTXT$text, "\n") #Split by newline character
    thisText="" #Create empty string to hold our new file contents
    for(q in 1:length(text[[1]])){ #Iterate through lines
      thisText=paste(thisText, text[[1]][q], sep="") #Add without newlines
      if(q==1){
        thisText=paste(thisText, "\n", sep="") #Only add newline after header
      }
    }
    cat(thisText, file=x) #Write new file contents to old file
  }
}
