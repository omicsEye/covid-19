# Code based on ade4:::mantel.rtest

library(permute)
library(ade4)
library(ggplot2)

Nperms <- 4999

mantelbootstrap <- function(m1, m2, nrepet) {
  s1 <- unclass(m1)
  s2 <- unclass(m2)
  permi <- matrix(1:nrepet, nrow = nrepet, ncol = 1)
  bss <- apply(permi, 1, function(i) {
    bsi <- sample(seq_along(s1), length(s1), replace=T)
    return (cor(s1[bsi], s2[bsi]))
  })
  return (bss)
}


# Plots
library(RColorBrewer)
manteltest_plot <- function(O, P, Ocil, Ociu, title, stars=F) {
  
  library(reshape2)
  library(viridis)
  df <- melt(O)
  colnames(df) <- c("V1", "V2", "obs")
  df$V2 <- factor(df$V2, levels=rev(levels(df$V2)))
  df$obs <- df$obs^2
  if (!missing(P)) {
    Padj <- matrix(p.adjust(P, method="fdr"), nrow=nrow(P))
    df$padj <- melt(P)$value
    print(Padj)
  }
  if (!missing(Ocil)) {
    df$cil <- melt(Ocil)$value^2
    df$cil[melt(Ocil)$value < 0] <- 0
    df$ciu <- melt(Ociu)$value^2
  }
  
  # Increase the contrast of the color scale where it matters
  alpha <- 1.9
  beta <- 0.1
  colorvalues <- pbeta(seq(0, 1, length=101), alpha, beta)
  df$lblcolor <- ifelse(qbeta(df$obs, alpha, beta) < 0.8, "black", "white")
  
  ggp <- ggplot(data=df, aes(x=V1, y=V2)) +
    geom_tile(aes(fill=obs)) +
    scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
                         values=colorvalues,
                         na.value="white", limits=c(0, 1), name="Variance Explained")
  #ggp
  nudge_y <- 0
  if (!missing(P)) {
    ggp <- ggp + if (stars) {
      nudge_y <- -0.14
      geom_text(aes(label=ifelse(padj<=0.001, "***", ifelse(padj<=0.01, "**", ifelse(padj<=0.05, "*", ""))),
                    color=lblcolor),
                size=2, nudge_y=0.175, fontface="bold")
    } else {
      geom_text(aes(label=ifelse(is.na(padj), "", sprintf("FDR p\n%.2g", padj)),
                    color=lblcolor), size=1.7)
    }
  }
  ggp <- ggp + if (!missing(Ocil)) {
    geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%\n[%.1f%% - %.1f%%]", 100*obs, 100*cil, 100*ciu)),
                  color=lblcolor),
              size=1.7, nudge_y=nudge_y)
  } else {
    geom_text(aes(label=ifelse(is.na(obs), "", sprintf("%.1f%%", 100*obs)),
                  color=lblcolor),
              size=1.7, nudge_y=nudge_y)
  }
  ggp <- ggp +
    geom_text(aes(label=ifelse(V1!=V2, ifelse(is.na(obs), "", ""), "")), size=1.7, color="gray") +
    theme_nature() +
    theme(axis.text.x=element_text(angle=-17, hjust=0.05),
          legend.position = "left", axis.ticks.y = element_blank()) +
    scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0), position = "right") +
    scale_color_manual(values=c(white="white", black="black")) +
    guides(fill=guide_colourbar(title=NULL, barheight=unit(0.65,"npc"), label.position = "left"), color="none") +
    xlab(NULL) + ylab(NULL)
   if (!missing(title)) {
       ggp <- ggp + ggtitle(title)
   }
  
  return (ggp)
}


#' PERMANOVA results "heatmap"
#'
#' @param R2 A matrix of R2 values from \code{adonis}.
#' @param P A matrix of P-values from \code{adonis}.
#' @param fontsize The desired font size of the % R2 values in the heatmap.
#' @param FDR If \code{T}, P-values are first BH-adjusted, and significance
#' stars are shown as *** < 0.001, ** < 0.01, * < 0.05.
#' @return ggplot object.
PERMANOVA_heatmap <- function(R2, P, fontsize=6, FDR=T, alpha=NA, beta=NA) {
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
    theme_nature() +
    theme(axis.text.x = element_text(angle=-17, hjust=0),
          panel.border=element_rect(fill=NA),
          legend.position = "left", axis.ticks.y = element_blank())
  
  return (ggp)
}



#####metadata cleaning

#######
metadata <- read.delim(
  'data/MetaData/genome_metadata_formatted.tsv',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

# select the metadata of interest
metadata["EPI_ISL_436426", "date"] <- "2020-04-05"
metadata$date[metadata$date=="2020"] <- "2020-01-01"
metadata$date[metadata$date=="2020-01"] <- "2020-01-01"
#metadata$date <- as.Date(metadata$date, format="%yyyy-%mm-%dd")
start <- "2019-12-07" # as.Date("2019-12-07", format="%yyyy-%mm-%dd")
Transmission_date <- as.Date(metadata$date)-as.Date(start)
metadata$Days <- Transmission_date
transmission_data <- metadata$date
metadata2 <- metadata[, c("country", "country_exposure", "sex", "age", "Days")]

# replace charecters to be properly NAs
metadata2[metadata2 == '?'] <- NA
metadata2[metadata2 == 'nan'] <- NA

metadata_order = colnames(metadata2)
comp_meta <- metadata2[complete.cases(metadata2), ]


#####
sample_names <- rownames(comp_meta)




# This program calaucltes distance matrix from MSA files for different regions of CoV genome
library(ape)
library(vegan)
library(stringr)

setwd("~/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/hackensack")
regions <- list(
  "S", "E", "M", 'N', 'ORF1a',
  'ORF1b', 'ORF1ab', 'ORF3a',
  'ORF6', 'ORF7a', 'ORF7b',
  'ORF8', 'ORF10', "2'-O-ribose methyltransferase",
  "3C-like proteinase", "3'-to-5' exonuclease",
  "endoRNAse", "helicase", "leader protein",
  "nsp2", "nsp3", "nsp4", "nsp6",
  "nsp7", "nsp8", "nsp9", "nsp10",
  "nsp11", "RNA-dependent RNA polymerase") #'reference'
#num_samples <- 5000
#samples <- sample(1:15000, num_samples, replace=F)
#similarity_method <- 'TN93' #'K80'

# Creat a dirctory for MaAsLin input files
output_path <- 'analysis/'
dir.create(file.path(output_path), showWarnings = TRUE)

omics <- c(
  "Genome", "S", "E", "M", 'N', 'ORF1a',
  'ORF1b', 'ORF1ab', 'ORF3a',
  'ORF6', 'ORF7a', 'ORF7b',
  'ORF8', 'ORF10', "2'-O-ribose methyltransferase",
  "3C-like proteinase", "3'-to-5' exonuclease",
  "endoRNAse", "helicase", "leader protein",
  "nsp2", "nsp3", "nsp4", "nsp6",
  "nsp7", "nsp8", "nsp9", "nsp10",
  "nsp11", "RNA-dependent RNA polymerase"
)

M <- matrix(NA, nrow=length(omics), ncol=length(omics))
rownames(M) <- omics
colnames(M) <- omics
mt_inter_doa <- NA
mt_inter_doa <- list(C=M, Cil=M, Ciu=M, P=M)
Nperms <- 4999
#i =1
#j =2
#mt_inter_doa <- NULL

for (i in seq_along(omics)){
  for (j in  seq_along(omics)){
    if(i>=j)next 
    D1 <- read.delim(
      paste('data/TN93/', omics[i],'_74.tsv', sep=''), 
      sep = '\t',
      header = TRUE,
      fill = F,
      comment.char = "" ,
      check.names = F,
      row.names = 1
    )
    rownames(D1) = colnames(D1) = gsub(paste(" ",omics[i], sep=''), "",rownames(D1))
    D2 <- read.delim(
      paste('data/TN93/', omics[j],'_74.tsv', sep=''), 
      sep = '\t',
      header = TRUE,
      fill = F,
      comment.char = "" ,
      check.names = F,
      row.names = 1
    )
    rownames(D2) = colnames(D2) = gsub(paste(" ",omics[j], sep=''), "",rownames(D2))
    intersect_samples <- intersect(rownames(D1), rownames(D2))
    D1 = as.matrix(D1) #[intersect_samples, intersect_samples]
    D2 = as.matrix(D2) # [intersect_samples, intersect_samples]
    #if (all(D1 ==0 ) | all(D2 ==0 )){
    #  mt_inter_doa$C[i,j] <- 0
    #  mt_inter_doa$P[i,j] <- 0
    #}else 
    tryCatch({
      mt <- mantel.rtest(as.dist(D1), as.dist(D2), nrepet=Nperms)
      #mt$bootstraps <- mantelbootstrap(as.dist(D1), as.dist(D2), nrepet=Nperms)
      mt_inter_doa$C[i,j] <- mt$obs
      #mt_inter_doa$Cil[i,j] <- quantile(mt$bootstraps, 0.025, na.rm=T)
      #mt_inter_doa$Ciu[i,j] <- quantile(mt$bootstraps, 0.975, na.rm=T)
      mt_inter_doa$P[i,j] <- mt$pvalue
    },
    error = function(e) {
      print(names(omics)[i])
      print(names(omics)[j])
      print(paste('error:', e))
    })
  }
}
C <- t(mt_inter_doa$C)
P <- t(mt_inter_doa$P)
#C[is.na(C)] <- 0
#P[is.na(P)] <- ""
source("/Users/rah/Documents/omicsEye/hmp2_analysis/common/theme_nature.r")
mantel_plots <- manteltest_plot(t(mt_inter_doa$C), t(mt_inter_doa$P), stars = T)
#mantel_plots <- mantel_plots + theme_nature()
pdf("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/overview_mantel_hs_v3.png", 7.5, 4)#5.2, 2.1)
print(mantel_plots)
dev.off()

for (i in seq_along(omics)){
    D <- read.delim(
      paste('data/TN93/', omics[i],'_74.tsv', sep=''), 
      sep = '\t',
      header = TRUE,
      fill = F,
      comment.char = "" ,
      check.names = F,
      row.names = 1
    )
    rownames(D) = colnames(D) = gsub(paste(" ",omics[i], sep=''), "",rownames(D))
    rownames(D) = colnames(D) = gsub("NJ-Bio", "",rownames(D))
    rownames(D) = colnames(D) = gsub("-", "",rownames(D))
    rownames(D) = colnames(D) = gsub("Ampliseq", "",rownames(D))
    write.table(D,  paste('data/TN93/', omics[i],'.tsv', sep=''),
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
}

comm_samples <- read.delim('data/metadata_omnibus.txt', 
  sep = '\t',
  header = TRUE,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
colnames(D) %in% row.names(comm_samples)
M <- matrix(NA, nrow=length(omics), ncol=length(metadata_order))
rownames(M) <- names(omics)
colnames(M) <- metadata_order
Nperms <- 4999
R = P = M
library(vegan)
library(permute)

for (i in seq_along(omics)){
  print(names(omics)[i])
  D = regions[[names(omics)[i]]]
  tryCatch({
    D2 <- as.matrix(D)
    #D2 <- D2[,complete.cases(D2)]
    comm_samples <- intersect(rownames(comp_meta), rownames(D2))
    sub_comp_meta <- comp_meta[comm_samples, metadata_order, drop = F]
    #D <- as.dist(D)
    D2 <- D2[comm_samples, comm_samples]
    #D2 <- as.dist(D2)
    # Test statistic from non-permuted data
    D2[is.na(D2)] <- 0
    ad <- adonis(D2 ~ ., data=sub_comp_meta, permutations=Nperms)
    temp_R <- ad$aov.tab$R2
    names(temp_R) <- rownames(ad$aov.tab)
    n <- length(temp_R)
    temp_P <- ad$aov.tab$`Pr(>F)`
    R[names(omics)[i],] <- temp_R[1:(n-2)]
    P[names(omics)[i],] <- temp_P[1:(n-2)]
  },
  error = function(e) {
    print(names(omics)[i])
    print(paste('error:', e))
  })
}

pdf("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/overview_omnibus_tests_noFDR_3.pdf", 4.5, 1.75)#5.2, 2.1)
p <- PERMANOVA_heatmap(t(R), t(P), FDR = F, alpha=1.9, beta=0.1)
print(p)
dev.off()


