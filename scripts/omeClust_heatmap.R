library(ape)

metadata2 <- read.table(
    'data/MetaData/metadata.tsv',
    sep = '\t',
    header = TRUE,
    fill = FALSE,
    comment.char = "" ,
    check.names = FALSE,
    row.names = 1
  )
meta <- colnames(metadata2)
stats_table <-
  setNames(
    data.frame(matrix(
      ncol = dim(metadata2)[2], nrow = length(regions)
    )),
    colnames(metadata2)
  )
rownames(stats_table) <- regions
num_samples <- 500
similarity_method <- 'K80' #'TN93'

# Creat a dirctory for MaAsLin input files
input_path <- paste('analysis/m2clust_0508/',similarity_method, sep = '')
region <- S
for (region in regions) {
  for (meta in colnames(metadata2)){
  m2clust_output <- paste(input_path, '/', 'E', '_',500,'/clusters.txt', sep='')
  # m2clust_overview <- read.table(
  #   m2clust_output,
  #   sep = '\t',
  #   header = TRUE,
  #   fill = FALSE,
  #   comment.char = "" ,
  #   check.names = FALSE,
  #   row.names = 1
  # )
  stats_table[region, meta] <- runif(1)
  }
}
stats_table <- t(stats_table)
max_value <- ceiling(max(stats_table))
min_value <- ceiling(min(stats_table))
range_value <- max(c(abs(max_value),abs(min_value)))
breaks <- seq(range_value, range_value, by = 0.1)
color = colorRampPalette(c("WhiteSmoke", "orange", "darkorange"))

p <- pheatmap::pheatmap(stats_table,
                        cellwidth = 7,
                        cellheight = 7,
                        # changed to 3
                        #main = title,
                        fontsize = 7,
                        kmeans_k = NA,
                        border = TRUE,
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        scale = "none",
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        clustering_distance_rows = "euclidean",
                        clustering_distance_cols = "euclidean",
                        legend = TRUE,
                        border_color = 'grey93',
                        color = color(range_value*100),
                        #breaks = breaks,
                        treeheight_row = 0,
                        treeheight_col = 0,
                        display_numbers = matrix(ifelse(
                          stats_table > 0.75, "*",  ""), nrow(stats_table))
)
ggsave(filename='manuscript/Figures/fig2/m2clust_hetampa.pdf', plot=p, width = 60, height = 35, units = "mm", dpi = 300) #width = 183, height = 126.8,
 