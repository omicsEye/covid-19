library(ape)
setwd("~/Box/COVID19_Project")
metadata <- read.delim(
  'data/Distances_2000_GISAID/Distances/GISAID_2000_Metadata.txt',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
meta <- colnames(metadata)
regions <- c("Genome",
           'ORF1ab',
           "3C-like proteinase",
           'ORF1a',
           "leader protein","nsp2", "nsp3", "nsp4", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp11",
           'ORF1b',
           "RNA-dependent RNA polymerase",
           "helicase",
           "3'-to-5' exonuclease",
           "endoRNAse",
           "2'-O-ribose methyltransferase",
           "S",
           'ORF3a',
           "E", "M",
           'ORF6',
           'ORF7a', 
           'ORF7b',
           'ORF8',
           'N',  
           'ORF10'
)
stats_table <-
  setNames(
    data.frame(matrix(
      ncol = dim(metadata)[2], nrow = length(regions)
    )),
    colnames(metadata)
  )
rownames(stats_table) <- regions

similarity_method <- 'TN93' #'K80' #

# Creat a dirctory for MaAsLin input files
region = "S"
input_path <- '/Users/rah/Documents/omeClust2000_TN93/'
for (region in regions) {
  omeClust_output <- paste(input_path, region,'/clusters.txt', sep='')
  tryCatch({
    omeClust_overview <- read.table(
      omeClust_output,
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      row.names = 1
    )
    stats_table[region, meta] <- omeClust_overview[1,meta]
  },
  error = function(e) {
    print(region)
    message(e)
    #stats_table[region, meta] <- c(NA,NA, NA)
  })
  
}
write.table( stats_table,'/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/omeClust_NCBI_metadata.txt',
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

stats_table <- stats_table[complete.cases(stats_table),] #
stats_table <- t(stats_table) 
pdf("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/omeClust_NCBI_metadata.pdf", 4, 6)#5.2, 2.1)
p <- heatmap3::heatmap3(stats_table, scale = 'col')
print(p)
dev.off()

p <- gplots::heatmap.2(stats_table)
#stats_table <- stats_table *100
Scaled_stats_table <- scale(stats_table)
max_value <- ceiling(max(Scaled_stats_table))
min_value <- ceiling(min(Scaled_stats_table))
range_value <- max(c(abs(max_value),abs(min_value)))
breaks <- seq(range_value, range_value, by = .1)
color = colorRampPalette(c("WhiteSmoke", "orange", "darkorange"))
p <- pheatmap::pheatmap( Scaled_stats_table,
                        cellwidth = 5,
                        cellheight = 5,
                        # changed to 3
                        #main = title,
                        fontsize = 6,
                        #kmeans_k = NA,
                        #border = TRUE,
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        #scale = "column",
                        cluster_rows = FALSE,
                        cluster_cols = TRUE,
                        clustering_distance_rows = "euclidean",
                        clustering_distance_cols = "euclidean",
                        legend = TRUE,
                        #border_color = 'grey93',
                        #color = color(range_value),
                        #breaks = breaks,
                        treeheight_row = 20,
                        treeheight_col = 6)
                        #display_numbers = matrix(ifelse(
                       #   stats_table > 0.75, "",  ""), nrow(stats_table))
ggsave(filename='/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/manuscript/Figures/fig3/omeClust_GISAID_metadata_genome_order.pdf', plot=p, width = 90, height = 80, units = "mm", dpi = 350) #width = 183, height = 126.8,
 