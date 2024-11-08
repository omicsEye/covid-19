library(g3viz)
library(vroom)
library(Rfast)
library(ggplot2)
setwd("~/Documents/")
table=vroom("ncbivcf0.00005MAF.tsv", skip = 3)
table_mut <- table[,11:dim(table)[2]]
table_mut[!(table_mut==1) && !(table_mut==0) ] <- 1

myData=as.data.frame(table$POS)

myData$count=rowSums(table_mut)
colnames(myData)=c("Position","Count")
#myData=myData[-c(1,2,3),]
#myData=myData[-217,]
#myData=myData[-214,]
#png(filename="Lollipop.pdf",height=2000,width=10000,res=800)
#myData[myData<200]<-NA
myData=vroom("myData.csv")
lollipop_plot <- ggplot(myData, aes(x=Position, y=log(Count))) +
  geom_segment( aes(x=Position, xend=Position, y=0, yend=log(Count)), size = 0.005) +
  geom_point( fill = 'darkolivegreen4',
              color = 'black',
              alpha = .5,
              shape = 21,
              size = .65,
              stroke = 0.05) +
  #theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("Base pair position with mutation") +
  ylab("Number of genomes with mutation (log)") +
  geom_hline(yintercept=log(8000), color="red", size = 0.05) +
  omicsArt::theme_omicsEye()
lollipop_plot
ggsave(filename = 'Lollipop_tyson_log.png', plot = lollipop_plot,
       height = 1.75,  width = 7.2 ,units = "in", dpi = 350)
dev.off()
