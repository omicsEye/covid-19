
#require(gdata)
library(gdata)
library(pheatmap)
library(vegan)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(Hmisc)
library(Maaslin2)
library(readxl)
library(deepath)
library(dplyr)
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
library(ape)
region = "E"
 
recombination_data <- read.delim(
  '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/data/recomination_data.txt',
  sep = '\t',
  header = TRUE,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

library(Hmisc)
#recombination_data$Number[recombination_data$Number<0] <- NA
recombination_data <- recombination_data[, c("Site", "Number")]
recombination_data_rate <- as.data.frame(table(recombination_data$Site), drop = F)
colnames(recombination_data_rate) <- c("Site", "Number")
rownames(recombination_data_rate) <- recombination_data_rate$Site
recombination_data_rate$Size = NA
for (region in omics) {
  #if (!region %in% rownames(recombination_data_rate)) next
  print(region)
  seq <- read.FASTA(paste('~/Box/COVID19_Project/data/Filtered_2000_GISAID_Subset/',region, '.fasta', sep=''))
  if (region == "3C-like proteinase")
    recombination_data_rate["3C", "Size"] <- length(seq[[2]])
  else if (region == "endoRNAse")
    recombination_data_rate["endornase", "Size"] <- length(seq[[2]])
  else if (region == "3'-to-5' exonuclease")
    recombination_data_rate["exonuclease", "Size"] <- length(seq[[2]])
  else if (region == "leader protein")
    recombination_data_rate["leader", "Size"] <- length(seq[[2]])
  else if (region == "2'-O-ribose methyltransferase")
    recombination_data_rate["methyltransferase", "Size"] <- length(seq[[2]])
  else if (region == "RNA-dependent RNA polymerase")
    recombination_data_rate["RdRp", "Size"] <- length(seq[[2]])
  else
    recombination_data_rate[region, "Size"] <- length(seq[[2]])
}
recombination_data_rate$Rate <- recombination_data_rate$Number/recombination_data_rate$Size
recombination_data_rate <- recombination_data_rate[!is.na(recombination_data_rate$Site),]
library(viridis)
library(RColorBrewer)
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = sort(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
box_plot <-
  ggplot2::ggplot(
    data = recombination_data_rate, ggplot2::aes(x=reorder(Site, Number, median, na.rm = TRUE), Rate)) +
  ggplot2::geom_bar(stat = 'identity',
    #ggplot2::aes(fill = Site),
    fill = "Blue",
    outlier.alpha = 0.0,
    na.rm = TRUE,
    alpha = .5,
    show.legend = FALSE,
    lwd=.1
  ) +
  # ggplot2::geom_point(
  # #  #ggplot2::aes(fill = Site),
  #   ggplot2::aes(fill = "Blue"),
  #   #fill = "Blue",
  #   alpha = 0.75 ,
  #   size = .5,
  #   shape = 21,
  #   stroke = 0.05,
  #   color = 'black',
  #   position = ggplot2::position_jitterdodge()
  # ) +
  ggplot2::scale_fill_manual(values = col_vector) +  guides(colour = guide_legend(override.aes = list(size=3))) +
  geom_text(aes(x= Site,
                y = Rate + .03,
                label=Number),
            size = 1
  )
box_plot <- box_plot + theme_omicsEye() + ggplot2::xlab("") + ggplot2::ylab("Selection proportion") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
box_plot <- box_plot +  guides(fill = guide_legend(override.aes = list(size=5))) + theme(legend.position =  "none", panel.border = element_rect(colour = "black", fill=NA, size=.25)) #
box_plot
ggsave(filename = '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/analysis/Recombination_rate.pdf', plot = box_plot, height = 2,  width = 2.5 ,units = "in", dpi = 350)
dev.off()
