

library(dplyr)
library(reshape2)
setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/Coronavirus")
# You can effectively remove scientific notation in printing with this code:
options(scipen = 999)

genome_path <-  'data/ORF3a_dist_k80.tsv'
within_corr_threshold <- .9 # E threshold .85
genome <- NA
genome_data <- NA
genome <- read.table(
  genome_path,
  header = T,
  row.names = 1,
  sep = "\t",
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE
)
#rownames(D1) <- abbreviate(rownames(D1), minlength=20)
#colnames(D1) <- abbreviate(colnames(D1), minlength=20)
genome <- 1.00 - abs(genome)
genome[lower.tri(genome, diag = TRUE)] <- NA
genome_data <-
  melt(genome %>% mutate(from = rownames(.)), id.vars = 'from')
genome_data <-
  setNames(genome_data, c('from', 'to', 'weight'))
genome_data <-
  genome_data[abs(genome_data$weight) > within_corr_threshold, ]

genome_data <- genome_data[genome_data$from != genome_data$to,]
#x_y_data <- x_y_data[complete.cases(x_y_data),]


genome_data <- genome_data[!is.na(genome_data$weight), ]
x_y_data <- genome_data

#then call viz function script
#-------------vizNetwrok
library(visNetwork)

label = c(unlist(unique(union(x_y_data$from, x_y_data$to))))
id = seq(length(unlist(unique(union(x_y_data$from, x_y_data$to)))))
nodes = data.frame(id, label)
colnames(nodes) <- c("id", "label")

#id has to be the same like from and to columns in edges
nodes$id <- nodes$label
nodes$title <- nodes$label
#Edges
edges <- x_y_data[,1:3]# as.data.frame(lesmis[1])
colnames(edges) <- c("from", "to", "width")

#Create graph for Louvain
graph <- igraph::graph_from_data_frame(edges, directed = FALSE)


#Louvain Comunity Detection
cluster <- igraph::cluster_louvain(graph)

cluster_df <- as.data.frame(as.list(igraph::membership(cluster)))
cluster_df <- as.data.frame(t(cluster_df))
cluster_df$label <- label #rownames(cluster_df)

#Create group column
nodes <- left_join(nodes, cluster_df, by = "label")
colnames(nodes)[4] <- "group"

loaded_data <- load_data(input = 'data/Metabolites_netome_format.xlsx',
                         type = 'all', sheet = 1, ID = 'Metabolite')
features_info <- loaded_data$feature_metadata
nodes <- nodes[nodes$label %in% features_info$Metabolite,]
edges <- edges[edges$from %in% features_info$Metabolite & edges$to %in% features_info$Metabolite ,]
nodes$pathway <- features_info[match(nodes$label, features_info$Metabolite),"SUPER_PATHWAY"] 
pathways_of_interest = c("Aminosugar Metabolism",
                         "Diacylglycerol",
                         "Sphingomyelins"
)

nodes$pathway[!nodes$pathway %in% pathways_of_interest] <- NA
nodes$group <- nodes$pathway
#nodes$pathway[is.na(nodes$pathway)] <- features_info[match(nodes$label, features_info$Metabolite),"SUPER_PATHWAY"]
#vizgraph <- visNetwork(nodes, edges)
#vizgraph
vizgraph2 <- visNetwork(nodes, edges, width = "100%") %>%
  visIgraphLayout() %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#0085AF",
      border = "#013848",
      highlight = "#000000"
    ),
    shadow = list(enabled = F, size = 10)
  ) %>%
  visEdges(
    shadow = FALSE,
    color = list(color = "#0085AF", highlight = "#C62F4B")
  ) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
             selectedBy = "group", nodesIdSelection = F) %>% 
  #visLegend(width = 0.1, position = "right", main = "community") %>%
  addFontAwesome() %>%
  addExport(pdf = TRUE)  %>%
  visLayout(randomSeed = 90) %>%
  visExport(type = "pdf", name = "export-network-community",
            float = "left", label = "Save network", background = "white", style= "") 
vizgraph2
write.table(
  nodes,
  'analysis/ORF3a-community-network.txt',
  sep = "\t",
  eol = "\n",
  quote = F,
  col.names = NA,
  row.names = T
)
# same as
visSave(vizgraph2, , file = "./genome_network.pdf", background = "white")
dev.off()