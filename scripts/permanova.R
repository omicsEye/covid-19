

#' PERMANOVA results "heatmap"
#'
#' @param R2 A matrix of R2 values from \code{adonis}.
#' @param P A matrix of P-values from \code{adonis}.
#' @param fontsize The desired font size of the % R2 values in the heatmap.
#' @param FDR If \code{T}, P-values are first BH-adjusted, and significance
#' stars are shown as *** < 0.001, ** < 0.01, * < 0.05.
#' @return ggplot object.
PERMANOVA_heatmap <- function(R2, P, fontsize=5, FDR=T, alpha=NA, beta=NA) {
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  
  df <- melt(R2)
  colnames(df) <- c("Feature", "Dataset", "R2")
  df$P <- melt(P)[,3]
  if (FDR) {
    df$P <- p.adjust(df$P, method="fdr")
  }
  df$varExpPct <- sprintf("%.1f%%", 100*df$R2)
  df$varExpPct[is.na(df$R2)] <- ""
  df$NAtext <- ""
  df$NAtext[is.na(df$R2)] <- "N/A"
  df$stars <- cut(df$P, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
  df$stars[is.na(df$R2)] <- ""
  
  # Try to make a reasonable color scheme that has contrast where needed
  colors <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100)
  R2Q <- quantile(R2, c(0.25, 0.75), na.rm=T)
  if (is.na(alpha) || is.na(beta)) {
    labhat <- optim(par=c(0, 0), method="Nelder-Mead",
                    fn=function(lab) sum((pbeta(c(0.25, 0.75), exp(lab[1]), exp(lab[2])) - R2Q)^2))
    abhat <- exp(labhat$par)
    colorvalues <- pbeta(seq(0, 1, length=length(colors)), abhat[1], abhat[2])
    
    # Label colors flip to white when the color is too dark
    df$lblcolor <- ifelse(qbeta(df$R2, abhat[1], abhat[2]) < 0.8, "black", "white")
    #cat(sprintf("Best-fit alpha = %g, beta = %g\n", abhat[1], abhat[2]))
  } else {
    colorvalues <- pbeta(seq(0, 1, length=length(colors)), alpha, beta)
    
    # Label colors flip to white when the color is too dark
    df$lblcolor <- ifelse(qbeta(df$R2, alpha, beta) < 0.8, "black", "white")
  }
  
  
  ggp <- ggplot(data=df, aes(x=Dataset, y=Feature)) +
    geom_tile(aes(fill=R2)) +
    geom_text(aes(label=varExpPct, color=lblcolor), size=fontsize/(14/5), nudge_y=-0.15) +
    geom_text(aes(label=NAtext), color="grey", size=fontsize/(14/5), nudge_y=-0.15) +
    geom_text(aes(label=stars, color=lblcolor), fontface="bold", size=1.25*fontsize/(14/5), nudge_y=0.12) +
    scale_fill_gradientn(colors=colors, values=colorvalues, limits=c(0, 1), na.value="white") +
    scale_color_manual(values=c(black="black", white="white")) +
    scale_x_discrete(expand=c(0,0)) + xlab(NULL) +
    scale_y_discrete(expand=c(0,0), position = "right", limits = rev(levels(df$Feature))) + ylab(NULL) +
    guides(color="none",
           fill=guide_colourbar(title=NULL, barheight=unit(40,"mm"), label.position = "left")) +
    theme_omicsEye() +
    theme(axis.text.x = element_text(angle=-17, hjust=0),
          panel.border=element_rect(fill=NA),
          legend.position = "left", axis.ticks.y = element_blank())
  
  return (ggp)
}
library(permute)
library(ade4)
library(ggplot2)
library(vegan)
library(stringr)
library(deepath)
setwd("~/Box/COVID19_Project")

####### read metadata
metadata <- read.delim(
  'data/Distances_2000_GISAID/Distances/GISAID_2000_Metadata.txt',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
omeClustClade_low <- read.delim(
  '/Users/rah/Documents/omeClust2000_TN93/Genome_low/feature_cluster_label.txt',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
omeClustClade_medium <- read.delim(
  '/Users/rah/Documents/omeClust2000_TN93/Genome_medium/feature_cluster_label.txt',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
omeClustClade_high <- read.delim(
  '/Users/rah/Documents/omeClust2000_TN93/Genome_high/feature_cluster_label.txt',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
metadata$omeClustClade_low <- NA
metadata$omeClustClade_medium <- NA
metadata$omeClustClade_high <- NA
metadata[row.names(omeClustClade_low), 'omeClustClade_low'] <- as.character(omeClustClade_low[,'Cluster'])
metadata[row.names(omeClustClade_medium), 'omeClustClade_medium'] <- as.character(omeClustClade_med[,'Cluster'])
metadata[row.names(omeClustClade_high), 'omeClustClade_high'] <- as.character(omeClustClade_high[,'Cluster'])

omics <- c("Genome",
           'ORF1ab',
           'ORF1a',
           "leader protein","nsp2", "nsp3", "nsp4", 
           "3C-like proteinase",
           "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "nsp11",
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

M_perm <- matrix(NA, nrow=length(omics), ncol=length(metadata))
rownames(M_perm) <- omics
colnames(M_perm) <- colnames(metadata)
Nperms <- 999
R_perm = P_perm = M_perm
library(vegan)
library(permute)
i <- 1
multivariable_test <- F
for (i in seq_along(omics)){
  print(omics[i])
  D <- read.delim(
    paste('/Users/rah/Documents/Distances_2000_GISAID/TN93/', omics[i],'.tsv', sep=''), 
    sep = '\t',
    header = TRUE,
    fill = F,
    comment.char = "" ,
    check.names = F,
    row.names = 1
  )
  intersect_samples <- intersect(rownames(D), rownames(metadata))
  metadata2 = metadata[intersect_samples,] 
  tryCatch({
    #D <- as.matrix(D)
    # Test statistic from non-permuted data
    D[is.na(D)] <- 0
    #if (!multivariable_test){ #do it one by one
    for (meta in c("Lineage", "omeClustClade_low", "omeClustClade_med", "omeClustClade_high")){
      ad <- adonis(D ~ . , data=metadata2[, meta, drop=F], permutations=Nperms)
      temp_R <- ad$aov.tab$R2
      names(temp_R) <- rownames(ad$aov.tab)
      n <- length(temp_R)
      temp_P <- ad$aov.tab$`Pr(>F)`
      R_perm[omics[i], meta] <- temp_R[1]
      P_perm[omics[i], meta] <- temp_P[1]
     # } 
    #}else{
      ad <- adonis(D ~ . , data=metadata2[, c("Days", "Region", "Country", "Country_exposure", "Sex", "Age"),
                                          drop=F], permutations=Nperms)
      temp_R <- ad$aov.tab$R2
      names(temp_R) <- rownames(ad$aov.tab)
      n <- length(temp_R)
      temp_P <- ad$aov.tab$`Pr(>F)` #   $`Pr(>F)`
      R_perm[omics[i],c("Days", "Region", "Country", "Country_exposure", "Sex", "Age")] <- temp_R[1:(n-2)]
      P_perm[omics[i],c("Days", "Region", "Country", "Country_exposure", "Sex", "Age")] <- temp_P[1:(n-2)]
    }
  },
  error = function(e) {
    print(names(omics)[i])
    print(paste('error:', e))
  })
}


## reorder
R_perm2 <- R_perm[,c("Country", "Country_exposure", "Sex", "Age", "Region", "Days" , 
  "Lineage", "omeClustClade_high","omeClustClade_med", "omeClustClade_low")]
P_perm2 <- P_perm[,c("Country", "Country_exposure", "Sex", "Age", "Region", "Days" , 
                     "Lineage", "omeClustClade_high","omeClustClade_med", "omeClustClade_low")]
metadata <-metadata[,c("Country", "Country_exposure", "Sex", "Age", "Region", "Days" , 
                         "Lineage", "omeClustClade_high","omeClustClade_medium", "omeClustClade_low")]
## renaming metadata
colnames(P_perm2) <- c("Country", "Country exposure", "Sex", "Age", "Region", "Days" , 
                       "Lineage", "omeClust clades (high)","omeClust clades (med)", "omeClust clades (low)")
colnames(R_perm2) <- c("Country", "Country exposure", "Sex", "Age", "Region", "Days" , 
                       "Lineage", "omeClust clades (high)","omeClust clades (med)", "omeClust clades (low)")
colnames(metadata) <- c("Country", "Country exposure", "Sex", "Age", "Region", "Days" , 
                       "Lineage", "omeClust clades (high)","omeClust clades (medium)", "omeClust clades (low)")
 
# writing to file
write.table( R_perm,'/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/R_perm_2000.txt',
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table( P_perm,'/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/P_perm_2000.txt',
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table( metadata,'/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/Metadata_2069.txt',
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)


pdf("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/overview_omnibus_tests_2000.pdf", 7.4, 1.75)#5.2, 2.1)
png("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/overview_omnibus_tests_2000.png", 7.4, 1.7, units ="in", res =1200)
p <- PERMANOVA_heatmap(t(R_perm2), t(P_perm2), FDR = F)#, alpha=1.9, beta=0.1)
print(p)
dev.off()

