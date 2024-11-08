
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

recombination_data <- read.delim(
  '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/data/recomination_data.txt',
  sep = '\t',
  header = TRUE,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

library(viridis)
box_plot <-
  ggplot2::ggplot(
    data = recombination_data, ggplot2::aes(x=reorder(Site,Number,na.rm = TRUE), Number)) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = Site),
    outlier.alpha = 0.0,
    na.rm = TRUE,
    alpha = .5,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = Site),
    alpha = 0.75 ,
    size = 1,
    shape = 21,
    stroke = 0.15,
    color = 'black',
    position = ggplot2::position_jitterdodge()
  ) +
  ggplot2::scale_fill_brewer(palette = viridis(30), direction=1)
box_plot <- box_plot + theme_omicsEye() + ggplot2::xlab("") + ggplot2::ylab("Recombination") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
box_plot
ggsave(filename = '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/analysis/Recombination.pdf', plot = box_plot, height = 2,  width = 7.2 ,units = "in", dpi = 350)
dev.off()
