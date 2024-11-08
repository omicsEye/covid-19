
library(dplyr)
library(reshape2)

# regions graph

C_2000 <- read.delim(
  '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/C_2000.txt', 
  sep = '\t',
  header = TRUE,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
#C_2000 <- as.data.frame(t(C_2000))
C_2000 <- C_2000 + min (C_2000, na.rm = T) * -1.0
D_2000 <- 1.0 - C_2000
for (i in 1: dim(D_2000)[1]){
  for (j in i:dim(D_2000)[2]){
    if (i==j)
      D_2000[i,j] = 0
    else
      D_2000[i,j] = D_2000[j,i]
  }
}
write.table(D_2000,
  '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/D_2000.txt', 
  sep = "\t",
  eol = "\n",
  quote = F,
  col.names = NA,
  row.names = T
)
# check fo rlower diagnal
#C_2000[lower.tri(C_2000, diag = TRUE)] <- NA
C_2000 <-
  melt(C_2000 %>% mutate(from = rownames(.)), id.vars = 'from')


# label dataframe 
C_2000 <-
  setNames(C_2000, c('from', 'to', 'weight'))
C_2000<- C_2000[!is.na(C_2000$weight),]
C_2000 <-
  C_2000[C_2000$weight >= 0.25, ]

C_2000$weight <- C_2000$weight * 10.0
x_y_data <-  C_2000



library(qgraph)
data(big5)
data(big5groups)
source('~/Documents/omicsEye/omicsPattern/R/utils.R')

# Correlations:
Q <- qgraph(cor(big5), minimum = 0.25, cut = 0.4, vsize = 1.5, groups = big5groups, 
            legend = TRUE, borders = FALSE)
title("Big 5 correlations", line = 2.5)

Q <- qgraph(x_y_corr_sub$r, minimum = 0.5, cut = 0.9, vsize = 1,  
            legend = F, borders = T)

# Libraries --------------------------------------------------------------

library(igraph)
#install.packages('~/Downloads/geomnet', repos = NULL, type = 'source', force =T)
library(geomnet)


# Data Preparation -------------------------------------------------------

#Load dataset
data(lesmis)


#Edges
#edges <- as.data.frame(lesmis[1])
#colnames(edges) <- c("from", "to", "weight")
edges <- x_y_data
#Create graph for the algorithms
g <- igraph::graph_from_data_frame(edges, directed = FALSE)

# Community Detection ----------------------------------------------------

# Louvain
lc <- igraph::cluster_louvain(g)
lc_com <- igraph::membership(lc)
#igraph::communities(lc)
plot(lc, g)

# Infomap
imc <- igraph::cluster_infomap(g)
igraph::membership(imc)
igraph::communities(imc)
plot(lc, g)
#igraph::plot.igraph(g, lc)

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
cluster <- igraph::cluster_louvain(g)

cluster_df <- as.data.frame(as.list(igraph::membership(lc)))
cluster_df <- as.data.frame(t(cluster_df))
cluster_df$label <- rownames(cluster_df)

#Create group column
nodes <- left_join(nodes, cluster_df, by = "label")
colnames(nodes)[4] <- "group"

#vizgraph <- visNetwork(nodes, edges)
#vizgraph
vizgraph2 <- visNetwork(nodes, edges) %>%
  visOptions(highlightNearest = T, nodesIdSelection = TRUE) %>% visLayout(randomSeed = 123) %>% visIgraphLayout() %>%
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
  visOptions(highlightNearest = list(enabled = T, degree = 0, hover = T),
             selectedBy = "group", nodesIdSelection = F) %>%  #
  #visLegend(width = 0.1, position = "right", main = "community") %>%
  addFontAwesome() %>%
  addExport(pdf = TRUE)  %>%
  visLayout(randomSeed = 90) %>%
  visExport(type = "pdf", name = "export-network-community",
            float = "left", label = "Save network", background = "white", style= "") 
vizgraph2
write.table(
  nodes,
  'data/community-network.txt',
  sep = "\t",
  eol = "\n",
  quote = F,
  col.names = NA,
  row.names = T
)
# same as
visSave(vizgraph2, , file = "./metabolites_network2.pdf", background = "white")
dev.off()
