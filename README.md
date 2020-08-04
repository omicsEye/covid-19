Steps for the analyses that we will do: <br/>
<details> 
  <summary> 0) Download metadata and MSA from GISAID</summary>
  a) First [register for an account](https://platform.gisaid.org/epi3/cfrontend#335368). This may take several days.
  
  b) Once you have an account, sign in [here](https://www.epicov.org/epi3/frontend#a3eb) with your username and password.
  
  c) From the EpiCov tab, click on `Downloads` and select the down arrow for the multiple sequence alignment (ex. MSA_0728), and the metadata (nextmeta)
</details> <br/>
<details>
  <summary> 1) Download NCBI fasta files</summary>
    a) Get these 4 different things
  </details> <br/>
  
  2) Format names of fasta files from GISAID and NCBI and combine portion of GISAID with NCBI into one fasta
  
  3) Format metadata from GISAID, retrieve metadata from NCBI, combine into one fasta
  
  4) Align GISAID and NCBI sequences using mafft
  
  5) separate both alignments by locus.
  
  6) make trees using raxml for 
  
    a) GISAID + NCBI (whole genome)
    
    b) GISAID + NCBI (each gene)
  
  7) Visualize trees and metadata with [iTOL](scripts/iTOL/iTOL.md)
  
  8) Make distance matrices for a) GISAID genome alignment, b) GISAID individual gene alignments, c) GISAID+NCBI genome alignment, d) GISAID+NCBI individual alignments
  
  9) Run m2clust on everything.

