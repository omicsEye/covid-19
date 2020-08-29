Steps for the analyses that we will do: <br/>
<details> 
  <summary> 0) Download metadata and MSA from GISAID</summary>
  a) First [register for an account](https://platform.gisaid.org/epi3/cfrontend#335368). This may take several days.
  
  b) Once you have an account, sign in [here](https://www.epicov.org/epi3/frontend#a3eb) with your username and password.
  
  c) From the EpiCov tab, click on `Downloads` and select the down arrow for the multiple sequence alignment (ex. MSA_0728), and the metadata (nextmeta)
  
  d) Extract from tar.xz file with `tar -xf file.tar.xz`
</details> <br/>
<details>
  <summary> 1) Download NCBI fasta files for each outgroup </summary>
  
   a) [Bat Coronavirus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SLen_i=29000%20TO%2040000&Completeness_s=complete&VirusLineage_ss=Bat%20SARS%20coronavirus%20HKU3,%20taxid:442736&VirusLineage_ss=Bat%20SARS-like%20coronavirus,%20taxid:1508227&VirusLineage_ss=Bat%20coronavirus,%20taxid:1508220)
  
   b) [MERS](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&Completeness_s=complete&VirusLineage_ss=Middle%20East%20respiratory%20syndrome-related%20coronavirus%20(MERS-CoV),%20taxid:1335626&SLen_i=28000%20TO%2035000)
    
   c) [SARS-related Coronavirus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome-related%20coronavirus,%20taxid:694009&CollectionDate_dr=2003-01-01T00:00:00.00Z%20TO%202019-10-29T23:59:59.00Z&SLen_i=25000%20TO%2035000)
    
   d) [Outgroups](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&SLen_i=28000%20TO%2035000&Completeness_s=complete&VirusLineage_ss=Transmissible%20gastroenteritis%20virus,%20taxid:11149)
   
  </details> <br/>
  
  2) Change names of GISAID fasta files to EPI IDs and [format metadata from GISAID](scripts/metadata_cleaning.R) 
  
  3) Retrieve and format [metadata from NCBI](scripts/NCBI_extract_metadata.R), combine files from NCBI and GISAID
  
  4) Align GISAID and NCBI sequences using mafft
  
  5) separate both alignments by locus.
  
  6) make trees using raxml for 
  
    a) GISAID + NCBI (whole genome)
    
    b) GISAID + NCBI (each gene)
  
  7) Visualize trees and metadata with [iTOL](scripts/iTOL/itol.md)
  
  8) Make distance matrices for a) GISAID genome alignment, b) GISAID individual gene alignments, c) GISAID+NCBI genome alignment, d) GISAID+NCBI individual alignments
  
  9) Run m2clust on everything.

