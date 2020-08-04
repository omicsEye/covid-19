## Using iTOL to create trees for Covid19

iTOL provides a classy way to display trees with metadata. Here is what we did. (note: apparently iTOL is planning on charging for use in 

0) Make an iTOL account. [Log in](https://itol.embl.de/login.cgi) and click "Upload tree files". Upload rooted tree in accepted format (Newick, Nexus, PhyloXML or Jplace)

    - If your tree isn't rooted, you can root it to a taa using figtree. Just upload the tree to figtree, click on the taxa or group of taxa you want to be the root, and click "Reroot" in the top left corner. Then, save the tree as a nexus.

1) Set display mode as circular or normal. I prefer circular with ignoring branch lengths.

2) Use [this script](itol.r) to prepare csvs for the scripts

3) Drag and drop scripts onto the tree
