library(tidyr)
library(dplyr)
library(reshape2)
source('~/Documents/omicsEye/omicsPattern/R/utils.R')
#setting the working directory
setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/")

## read proteins
protiens_maaslin2 <- read.delim(
  "/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/protiens_maaslin2_GLM_log/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
score_data_severe <- protiens_maaslin2[protiens_maaslin2$metadata=="Group" & protiens_maaslin2$value=="Severe" ,]
rownames(score_data_severe) <- score_data_severe$feature

score_data_non_severe <- protiens_maaslin2[protiens_maaslin2$metadata=="Group" & protiens_maaslin2$value=="non-Severe" ,]
rownames(score_data_non_severe) <- score_data_non_severe$feature

score_data_non_covid <- protiens_maaslin2[protiens_maaslin2$metadata=="Group" & protiens_maaslin2$value=="non-COVID-19" ,]
rownames(score_data_non_covid) <- score_data_non_covid$feature

number_of_sig_to_keep <- 60

# use score_data_severe is reference
order_sig <- rownames(score_data_severe)[1:number_of_sig_to_keep]
score_data_severe <- score_data_severe[order_sig,]
score_data_severe<- score_data_severe[order(score_data_severe$coef),]
order_sig <- rownames(score_data_severe)
score_data_severe <- within(score_data_severe,
                           feature <- factor(feature,
                                             levels=order_sig))
score_data_non_severe <- score_data_non_severe[rownames(score_data_severe),]
score_data_non_severe <- within(score_data_non_severe,
                            feature <- factor(feature,
                                              levels=order_sig))


score_data_non_covid <- score_data_non_covid[rownames(score_data_severe),]
score_data_non_covid <- within(score_data_non_covid,
                                feature <- factor(feature,
                                                  levels=order_sig))

severe_temp_diff_bar <- diff_bar_plot(score_data_severe, threshold = 1.0, pvalue_col = "pval", 
                               fdr ="qval", orderby = NA, x_label = 'coef', y_label = '')
non_severe_temp_diff_bar <- diff_bar_plot(score_data_non_severe, threshold = 1.0, pvalue_col = "pval", 
                                     fdr ="qval", orderby = NA, x_label = 'coef', y_label = '')
non_covid_temp_diff_bar <- diff_bar_plot(score_data_non_covid, threshold = 1.0, pvalue_col = "pval", 
                                         fdr ="qval", orderby = NA, x_label = 'coef', y_label = '')

#
## read proteins
metabolites_maaslin2 <- read.delim(
  "/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19/Analysis/meatbolites_maaslin2_CPLM_log/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
metabolites_score_data_severe <- metabolites_maaslin2[metabolites_maaslin2$metadata=="Group" & metabolites_maaslin2$value=="Severe" ,]
rownames(metabolites_score_data_severe) <- metabolites_score_data_severe$feature

metabolites_score_data_non_severe <- metabolites_maaslin2[metabolites_maaslin2$metadata=="Group" & metabolites_maaslin2$value=="non-Severe" ,]
rownames(metabolites_score_data_non_severe) <- metabolites_score_data_non_severe$feature

metabolites_score_data_non_covid <- metabolites_maaslin2[metabolites_maaslin2$metadata=="Group" & metabolites_maaslin2$value=="non-COVID-19" ,]
rownames(metabolites_score_data_non_covid) <- metabolites_score_data_non_covid$feature

number_of_sig_to_keep <- 60

# use score_data_severe is reference
order_sig <- rownames(metabolites_score_data_severe)[1:number_of_sig_to_keep]
metabolites_score_data_severe <- metabolites_score_data_severe[order_sig,]
metabolites_score_data_severe<- metabolites_score_data_severe[order(metabolites_score_data_severe$coef),]
order_sig <- rownames(metabolites_score_data_severe)
metabolites_score_data_severe <- within(metabolites_score_data_severe,
                            feature <- factor(feature,
                                              levels=order_sig))
metabolites_score_data_non_severe <- metabolites_score_data_non_severe[rownames(metabolites_score_data_severe),]
metabolites_score_data_non_severe <- within(metabolites_score_data_non_severe,
                                feature <- factor(feature,
                                                  levels=order_sig))


metabolites_score_data_non_covid <- metabolites_score_data_non_covid[rownames(metabolites_score_data_severe),]
metabolites_score_data_non_covid <- within(metabolites_score_data_non_covid,
                               feature <- factor(feature,
                                                 levels=order_sig))

metabolites_severe_temp_diff_bar <- diff_bar_plot(metabolites_score_data_severe, threshold = 1.0, pvalue_col = "pval", 
                                      fdr ="qval", orderby = NA, x_label = 'coef', y_label = '')
metabolites_non_severe_temp_diff_bar <- diff_bar_plot(metabolites_score_data_non_severe, threshold = 1.0, pvalue_col = "pval", 
                                          fdr ="qval", orderby = NA, x_label = 'coef', y_label = '')
metabolites_non_covid_temp_diff_bar <- diff_bar_plot(metabolites_score_data_non_covid, threshold = 1.0, pvalue_col = "pval", 
                                         fdr ="qval", orderby = NA, x_label = 'coef', y_label = '')


## do plots


fig2 <- ggdraw() +
  draw_plot(severe_temp_diff_bar+ theme(axis.title=element_text(face = 'bold')),
            x = 0, y = 0, width = .4, height = .49) +
  
  draw_plot(metabolites_severe_temp_diff_bar + theme(axis.title.x = element_blank()),
            x = 0, y = .51, width = .5, height = .49) +
  
  draw_plot(non_severe_temp_diff_bar + theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title=element_text(face = 'bold')),
            x = .7, y = 0, width = .3, height = .49) +
  
  draw_plot(metabolites_non_severe_temp_diff_bar + theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.x = element_blank()),
            x = .5, y = .51, width = .25, height = .49) +
  draw_plot(non_covid_temp_diff_bar  + theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title=element_text(face = 'bold')),
            x = .4, y = 0, width = .3, height = .49) +
  draw_plot(metabolites_non_covid_temp_diff_bar + theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.x = element_blank()),
            x = .75, y = .51, width = .25, height = .49) +

  draw_plot_label((label = c("Metabolites", "Proteins", "Severe", "non-Severe", "non-COVID")),
                  size = 5,x = c(0, 0, .4, .6, .8), y = c(1, .6,.5, .5, .5))
fig2

ggsave(filename = '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-Omics/manuscript/figures/fig2/fig2.pdf', plot=fig2, width = 183, height = number_of_sig_to_keep*4, units = "mm", dpi = 350)
ggsave(filename = '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/COVID-Omics/manuscript/figures/fig2/fig2.png', plot=fig2, width = 183, height = number_of_sig_to_keep*4, units = "mm", dpi = 350)

