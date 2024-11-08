
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
#library(omicsPattern)
fdr_threshold <- 1.00
setwd("~/Dropbox/Ali-Docs/Research_docs/Projects/COVID-19")
source('~/Documents/omicsEye/omicsPattern/R/utils.R')
source('~/Documents/omicsEye/omicsPattern/R/pcl_utils.R')

metadata <- read.table(
  'data/protein_metabloutes_Shen_etal/clinical_info.txt',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
metabolites <- read.table(
  'data/protein_metabloutes_Shen_etal/metabolites.txt',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
metabolites <- read.delim(
  'data/protein_metabloutes_Shen_etal/metabolites.txt',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)


metadata_sub <- metadata[, c("Sex","Age", "BMI")]
metadata_sub$BMI <- as.numeric(metadata_sub$BMI)
Maaslin2(metabolites, #'data/protein_metabloutes_Shen_etal/metabolites.txt',
         metadata_sub, #'data/protein_metabloutes_Shen_etal/clinical_info.txt',
         'analysis/meatbolites_maaslin2_CPLM_log_residuals_Group',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #fixed_effects = c("Sex"),
         analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         transform = 'LOG',
         normalization = 'NONE')
          #transform = "AST")#,
         #analysis_method = 'SLM',
         #)#,
         #transform = 'NONE',
         #min_abundance = 0.0,
         #min_prevalence = 0.0)

veg_dist <- as.matrix(vegdist(t(metabolites), method="bray", na.rm = T))

write.table(veg_dist, 'data/protein_metabloutes_Shen_etal/metabolites_bray.txt',
           sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(metadata_sub, 'data/protein_metabloutes_Shen_etal/metabolites_metadata.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
pcoa_plots <- ordplots(as.data.frame(t(metabolites)), metadata_sub, output = 'analysis/', outputname = 'tsne', method = 'tsne')
pcoa_plots <- ordplots(as.data.frame(t(metabolites)), metadata_sub, output = 'analysis/', outputname = 'pcoa', method = 'pcoa')

stats_table <- test2groups(data = as.data.frame(t(metabolites)),
            metadata= metadata_sub,
            meta = 'Group',
            case_label = 'Severe',
            control_label = 'Healthy',
            test_type = 'wilcox.test',
            paired = F)
stats_table2 <- stats_table[!is.na(stats_table$P.Value),]
diff_bar <- diff_bar_plot(stats_table2, threshold = 0.00001, method = 'fdr', orderby = 'logFC', y_label = '') + theme_nature()
ggsave(filename = 'Analysis/bar_plot_meatbolites.pdf', plot = diff_bar,
       height = 4,  width = 4 ,units = "in", dpi = 350)
       dev.off()
#Proteins cleaning
proteins_training <- read.delim(
 'data/protein_metabloutes_Shen_etal/protein_training.txt',
 sep = '\t',
 header = TRUE,
 fill = F,
 comment.char = "" ,
 check.names = F,
 row.names = 1
)
proteins_test <- read.delim(
  'data/protein_metabloutes_Shen_etal/protein_test.txt',
  sep = '\t',
  header = TRUE,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)
protiens <- merge(proteins_training, proteins_test, by = 0)
rownames(protiens) <- protiens$Row.names
protiens <- protiens[-1]
# remove columns with all NA
#proteins_training <- as.data.frame(proteins_training[,colSums(!is.na(proteins_training)) > 0])

# remove rows with all NA
#proteins_training <- as.data.frame(proteins_training[apply(proteins_training, 1, function(y) !all(is.na(y))),])

Maaslin2(protiens, #'data/protein_metabloutes_Shen_etal/metabolites.txt',
         metadata_sub, #'data/protein_metabloutes_Shen_etal/clinical_info.txt',
         'analysis/protiens_maaslin2_GLM_log',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #fixed_effects = c("Sex"),
         analysis_method = 'LM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         transform = 'LOG',
         normalization = 'NONE')

write.table(protiens, 'data/protein_metabloutes_Shen_etal/protiens.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
#transform = "AST")#,
#analysis_method = 'SLM',
#)#,
#transform = 'NONE',
#min_abundance = 0.0,
#min_prevalence = 0.0)

veg_dist <- as.matrix(vegdist(t(protiens), method="bray", na.rm = T))

write.table(veg_dist, 'data/protein_metabloutes_Shen_etal/protiens_bray.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(metadata_sub, 'data/protein_metabloutes_Shen_etal/metabolites_metadata.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
pcoa_plots <- ordplots(as.data.frame(t(protiens)), metadata_sub, output = 'analysis/', outputname = 'tsne_protiens', method = 'tsne')
pcoa_plots <- ordplots(as.data.frame(t(protiens)), metadata_sub, output = 'analysis/', outputname = 'pcoa_protiens', method = 'pcoa')

stats_table <- test2groups(data = as.data.frame(t(protiens)),
                           metadata= metadata_sub,
                           meta = 'Group',
                           case_label = 'Severe',
                           control_label = 'Healthy',
                           test_type = 'wilcox.test',
                           paired = F)
stats_table2 <- stats_table[!is.na(stats_table$P.Value),]
diff_bar <- diff_bar_plot(stats_table2, threshold = 0.00001, method = 'fdr', orderby = 'logFC', y_label = '') + theme_nature()
ggsave(filename = 'Analysis/bar_plot_meatbolites.pdf', plot = diff_bar,
       height = 4,  width = 4 ,units = "in", dpi = 350)
dev.off()

metadata_sub <- metadata[, c("Group"  ,"Sex","Age", "BMI")]
metadata_sub$Group <- as.character(metadata_sub$Group)
metadata_sub$Group[metadata_sub$Group == "non-Severe"] <- "1-non-Severe"
metadata_sub$Group <- as.factor(metadata_sub$Group)
metadata_sub$BMI <- as.numeric(metadata_sub$BMI)
class(metadata_sub$Age)
Maaslin2(metabolites, #'data/protein_metabloutes_Shen_etal/metabolites.txt',
         metadata_sub, #'data/protein_metabloutes_Shen_etal/clinical_info.txt',
         'analysis/meatbolites_maaslin2_CPLM_log_noSevere_ref',
         min_abundance = 0,
         min_prevalence = 0,
         max_significance = 0.1,
         #fixed_effects = c("Sex"),
         analysis_method = 'CPLM',
         #random_effects = c("GROUP NUMBER"),
         standardize = FALSE,
         transform = 'LOG',
         normalization = 'NONE')



tryCatch({
  path_results <- deepath::deepath(
                 input_data = stats_table,
                 input_metadata = NA,
                 meta <- 'TREATMENT',
                 case_label <- case_group,
                 control_label <- control_group,
                 output = output,
                 score_col = 'logFC',
                 pval_threshold = .05,
                 fdr_threshold = NA,
                 Pathway.Subject = NA,#'Metabolic',
                 do_plot = TRUE,
                 mapper_file = features_info,
                 method = 'wilcox',
                 min_member = 2,
                 pathway_col = "SUPER PATHWAY",
                 feature_col = "Metabolite")
}, error = function(err) {
  #invisible(dev.off())
  logging::logerror("Unable to do deepath or its plots!!!")
  logging::logerror(err)
})





