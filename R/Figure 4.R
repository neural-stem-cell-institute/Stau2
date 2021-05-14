
####Figure 4

library(dplyr)
library(tibble)
library(foreach)
library(hypeR)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpmisc)
library(stringr)
library(biomaRt)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(RColorBrewer)
library(GenomicAlignments)
library(foreach)
library(akmedoids)
library(GO.db)

######4a, b
#####get genelists for top stable GO BP categories cell cycle, chromosome organization, macromolecule localization from HypeR output (Fig.3 code)
#####get genelists for top dynamic GO BP categories neurogenesis, differentiation, and cell projection organization from HypeR output (Fig.3 code)

stable_cellcycle_scores <- inner_join(stable_cellcycle, scores, by= "Gene.Name")
stable_cellcycle_scores_sorted <- stable_cellcycle_scores[order(stable_cellcycle_scores$Summed.BF),] 
colnames (masterpeaks_trimmed_wm)[1]<- "Gene.Name"
stable_cellcycle_data <- inner_join(stable_cellcycle_scores_sorted[,c(1,6)], masterpeaks_trimmed_wm, 
                                    by= "Gene.Name")

cc<- stable_cellcycle_data[c(1,6:10),]  ####use genes unique to each category
cc2 <- melt(cc, id.vars=c("Gene.Name", "Summed.BF"))
p1<-ggplot(aes(variable, Gene.Name), data=cc2) 
p1<- p1+geom_point(aes (size=value,color=Summed.BF)) 
p1<- p1+ scale_color_gradient2(low="blue", mid="yellow", high= "orange")
p1<- p1 + theme(axis.text.x  = element_text(size=6.5))
p1<- p1 + theme(axis.text.x = element_text(angle = 90))
p1


stable_chromosome_scores <- inner_join(stable_chromosome, scores, by= "Gene.Name")
stable_chromosome_scores_sorted <- stable_chromosome_scores[order(stable_chromosome_scores$Summed.BF),] 
stable_chromosome_data <- inner_join(stable_chromosome_scores_sorted[,c(1,6)], masterpeaks_trimmed_wm, 
                                     by= "Gene.Name")

ch<- stable_chromosome_data[c(3,4,7:10),]
ch2 <- melt(ch, id.vars=c("Gene.Name", "Summed.BF"))
p1<-ggplot(aes(variable, Gene.Name), data=ch2) 
p1<- p1+geom_point(aes (size=value,color=Summed.BF)) 
p1<- p1+ scale_color_gradient2(low="blue", mid="white", high= "yellow", midpoint=1.16)
p1<- p1 + theme(axis.text.x  = element_text(size=6.5))
p1<- p1 + theme(axis.text.x = element_text(angle = 90))
p1

stable_macromolecule<- data.frame(read.csv("Stable_GOBP_macromolecule_catabolic_process.csv"))
stable_macromolecule_scores <- inner_join(stable_macromolecule, scores, by= "Gene.Name")
stable_macromolecule_scores_sorted <- stable_macromolecule_scores[order(stable_macromolecule_scores$Summed.BF),] 
stable_macromolecule_data <- inner_join(stable_macromolecule_scores_sorted[,c(1,6)], 
                                        masterpeaks_trimmed_wm, by= "Gene.Name")

mc<- stable_macromolecule_data[c(1,4,8,10,12,13),]
mc2 <- melt(mc, id.vars=c("Gene.Name", "Summed.BF"))
p1<-ggplot(aes(variable, Gene.Name), data=mc2) 
p1<- p1+geom_point(aes (size=value, color=Summed.BF)) 
p1<- p1+ scale_color_gradient2(low="blue", mid="white", high= "yellow", midpoint=1.16)
p1<- p1 + theme(axis.text.x  = element_text(size=6.5))
p1<- p1 + theme(axis.text.x = element_text(angle = 90))
p1

dec_dyn_neurogenesis<- data.frame(read.csv("Dec_dyn_GOBP_neurogenesis.csv"))
dec_dyn_neurogenesis_scores <- inner_join(dec_dyn_neurogenesis, scores, by= "Gene.Name")
dec_dyn_neurogenesis_scores_sorted <- dec_dyn_neurogenesis_scores[order(-dec_dyn_neurogenesis_scores$Summed.BF),] 
dec_dyn_neurogenesis_data <- inner_join(dec_dyn_neurogenesis_scores_sorted[,c(1,6)], 
                                        masterpeaks_trimmed_wm, by= "Gene.Name")

ng<- dec_dyn_neurogenesis_data[c(1,3,5,7,8,9),]
ng2 <- melt(ng, id.vars=c("Gene.Name", "Summed.BF"))
p1<-ggplot(aes(variable, Gene.Name), data=ng2) 
p1<- p1+geom_point(aes (size=value, color=Summed.BF)) 
p1<- p1+ scale_color_gradient(low="orange", high= "red")
p1<- p1 + theme(axis.text.x  = element_text(size=6.5))
p1<- p1 + theme(axis.text.x = element_text(angle = 90))
p1


dec_dyn_differentiation<- data.frame(read.csv("Dec_dyn_GOBP_neuron_differentiation.csv"))
dec_dyn_differentiation_scores <- inner_join(dec_dyn_differentiation, scores, by= "Gene.Name")
dec_dyn_differentiation_scores_sorted <- dec_dyn_differentiation_scores[order(-dec_dyn_differentiation_scores$Summed.BF),] 
dec_dyn_differentiation_data <- inner_join(dec_dyn_differentiation_scores_sorted[,c(1,6)], 
                                           masterpeaks_trimmed_wm, by= "Gene.Name")

nd<- dec_dyn_differentiation_data[9:15,]
nd2 <- melt(nd, id.vars=c("Gene.Name", "Summed.BF"))
p1<-ggplot(aes(variable, Gene.Name), data=nd2) 
p1<- p1+geom_point(aes (size=value, color=Summed.BF))
p1<- p1+ scale_color_gradient(low="orange", high= "red")
p1<- p1 + theme(axis.text.x  = element_text(size=6.5))
p1<- p1 + theme(axis.text.x = element_text(angle = 90))
p1

dec_dyn_projection<- data.frame(read.csv("Dec_dyn_GOBP_cell_projection_organization.csv"))
dec_dyn_projection_scores <- inner_join(dec_dyn_projection, scores, by= "Gene.Name")
dec_dyn_projection_scores_sorted <- dec_dyn_projection_scores[order(-dec_dyn_projection_scores$Summed.BF),] 
dec_dyn_projection_data <- inner_join(dec_dyn_projection_scores_sorted[,c(1,6)], 
                                      masterpeaks_trimmed_wm, by= "Gene.Name")

pr<- dec_dyn_projection_data[c(15:20),]
pr2 <- melt(pr, id.vars=c("Gene.Name", "Summed.BF"))
p1<-ggplot(aes(variable, Gene.Name), data=pr2) 
p1<- p1+ scale_color_gradient(low="orange", high= "red")
p1<- p1+geom_point(aes (size=value, color=Summed.BF)) 
p1<- p1 + theme(axis.text.x  = element_text(size=6.5))
p1<- p1 + theme(axis.text.x = element_text(angle = 90))
p1



all<- rbind(cc2, ch2, mc2, ng2, nd2, pr2)
all_2<-with(all, all[order(variable),])
all_2$Gene.Name<- factor (all_2$Gene.Name, levels=all_2[1:36,]$Gene.Name)
p1<-ggplot(aes(variable, Gene.Name), data=all_2) 
p1<- p1+ scale_color_gradient2(low="blue", mid="yellow", high= "red", midpoint=5)
p1<- p1+geom_point(aes (size=value, color=Summed.BF)) 
p1<- p1 + theme(axis.text.x  = element_text(size=6.5))
p1<- p1 + theme(axis.text.x = element_text(angle = 90))
p1


#4c,d
##### use weighted means calculated in Fig.3 code for making heatmaps in Figs. 4 and 5


#on (pro-IPC) and off (non-IPC) gene lists from Aprea et al. 2013
#use the weighted means, masterpeaks_wm  file
library(dplyr)
#Rename column 1 of masterpeaks_wm as "Gene.Name"
colnames (on_genes)<-"Gene.Name"
on_data_wm<-inner_join(on_genes, wm, by="Gene.Name")
#add BFs from Table S2

z_on<- data.frame (t(apply (on_data_wm[,2:13], 1, function(x)(x-mean(x))/sd(x))), row.names = on_data_wm$Gene.Name)
z_on<-setDT(z_on, keep.rownames = "Gene.Name")[]
#add BFs from Table S2
z_on_bf<- data.frame (inner_join(z_on, scores, by="Gene.Name"))
rownames(z_on_bf)<- z_on_bf$Gene.Name
z_on_bf$E11.mean<- rowMeans(z_on_bf[2:4],)
z_on_bf$E13.mean<- rowMeans(z_on_bf[5:7],)
z_on_bf$E15.mean<- rowMeans(z_on_bf[8:10],)
z_on_bf$E17.mean<- rowMeans(z_on_bf[11:13],)
z_on_bf<- z_on_bf[which(z_on_bf$Summed.BF>1.163),] #dynamic+partially dynamic genes

library(superheat)
library(RColorBrewer)
###generate heatmap with clustering #use rowmeans columns 
superheat.z.on.means.bf.clustered<- superheat(z_on_bf[,19:22], n.clusters.rows = 2, 
 clustering.method="hierarchical",left.label.size=0.3, left.label.text.size = 2.5, bottom.label.text.size = 3, 
      pretty.order.rows = TRUE,left.label = "variable", heat.pal = rev(brewer.pal(11,"RdYlBu")), 
      yr = z_on_bf$Summed.BF,  yr.plot.type = "bar", yr.axis.name = "ln(BF)")
write.csv(superheat.z.on.means.bf.clustered$membership.rows, "genes_z_on_superheat.csv")

##Repeated for off (non-IPC) genes.


