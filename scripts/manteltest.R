# Code based on ade4:::mantel.rtest

library(permute)
library(ade4)
library(ggplot2)

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
    df$padj <- melt(Padj)$value
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
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position = "left", axis.ticks.y = element_blank()) +
    scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0), position = "right") +
    scale_color_manual(values=c(white="white", black="black")) +
    guides(fill=guide_colourbar(title=NULL, barheight=unit(0.65,"npc"), label.position = "left"), color="none") +
    xlab(NULL) + ylab(NULL)
  # if (!missing(title)) {
  #     ggp <- ggp + ggtitle(title)
  # }
  
  return (ggp)
}



# Creat a dirctory for output files
output_path <- '/Users/rah/Documents/omeClust2000_TN93'
dir.create(file.path(output_path), showWarnings = TRUE)

omics <- c("Genome",
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

M <- matrix(NA, nrow=length(omics), ncol=length(omics))
rownames(M) <- omics
colnames(M) <- omics
mt_inter_doa <- NA
mt_inter_doa <- list(C=M, Cil=M, Ciu=M, P=M)
Nperms <- 1
i =1
j =2
#mt_inter_doa <- NULL

for (i in seq_along(omics)){
  print(omics[i])
  for (j in  seq_along(omics)){
    if(i>=j)next 
    D1 <- read.delim(
      paste('/Users/rah/Documents/Distances_2000_GISAID/TN93/', omics[i],'.tsv', sep=''), 
      sep = '\t',
      header = TRUE,
      fill = F,
      comment.char = "" ,
      check.names = F,
      row.names = 1
    )
    D2 <- read.delim(
      paste('/Users/rah/Documents/Distances_2000_GISAID/TN93/', omics[j],'.tsv', sep=''), 
      sep = '\t',
      header = TRUE,
      fill = F,
      comment.char = "" ,
      check.names = F,
      row.names = 1
    )
    intersect_samples <- intersect(rownames(D1), rownames(D2))
    D1 = D1[intersect_samples, intersect_samples] #as.matrix(D1) #
    D2 = D2[intersect_samples, intersect_samples] #as.matrix(D2) #
    D1[is.na(D1)] <- 0
    D2[is.na(D2)] <- 0
    
    #if (all(D1 ==0 ) | all(D2 ==0 )){
    #  mt_inter_doa$C[i,j] <- 0
    #  mt_inter_doa$P[i,j] <- 0
    #}else 
    tryCatch({
      mt <- mantel.rtest(as.dist(D1), as.dist(D2), nrepet=Nperms)
      #mt$bootstraps <- mantelbootstrap(as.dist(D1), as.dist(D2), nrepet=Nperms)
      mt_inter_doa$C[omics[i],omics[j]] <- mt$obs
      #mt_inter_doa$Cil[i,j] <- quantile(mt$bootstraps, 0.025, na.rm=T)
      #mt_inter_doa$Ciu[i,j] <- quantile(mt$bootstraps, 0.975, na.rm=T)
      mt_inter_doa$P[omics[i],omics[j]] <- mt$pvalue
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
write.table( C,'/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/C_2000.txt',
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table( P,'/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/P_2000.txt',
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
#C[is.na(C)] <- 0
#P[is.na(P)] <- ""
source("/Users/rah/Documents/omicsEye/hmp2_analysis/common/theme_nature.r")
#mantel_plots <- mantel_plots + theme_nature()
pdf("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/overview_mantel_C_2000_v3.pdf", 7.5, 4)#5.2, 2.1)
mantel_plots <- manteltest_plot(t(mt_inter_doa$C), t(mt_inter_doa$P), stars = T)
print(mantel_plots)
dev.off()




