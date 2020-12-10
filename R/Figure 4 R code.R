
####Figure 4
#4c, d


##### calculate weighted means for making heatmaps in Figs. 4 and 5
masterpeaks <- data.frame (read.csv("MasterPeaks.csv"), as.is=T) #MasterPeaks=Table S1
ids<-unique(as.character(masterpeaks$Gene.Name)) 
if(length(which(is.na(ids)))>=1) {ids<-ids[-(which(is.na(ids)))] }

timepoints<- c(1,2,3,4,5,6,7,8,9,10,11,12)# 1,2,3,etc correspond to E11 replicates E11.1, E11.2, E11.3 and so on
library(foreach)

wm<-foreach(i=1:length(ids),.combine="rbind") %do% {
  mat2<-masterpeaks[which(as.character(masterpeaks$Gene.Name)==ids[i]),]
  if(dim(mat2)[1]>1){
    wgt<-ecdf(0:max(mat2$Peak.Score)) 
    wgt<-sapply(mat2$Peak.Score,wgt) #apply the wgt from the ecdf function to the peak scores
  } else {wgt<-1}
  mat1<- mat2[,21:32] #dataframe containing E11-E17 gene values/counts
  if(dim(mat2)[1]>1){ 
    weight.mean<-foreach (j=1:12, .combine="c") %do% {
      mat3<- mat1[,which(timepoints==j)]
      sum (wgt*mat3)/sum (wgt) } } else {weight.mean<-foreach (j= 1:12, .combine="c") %do% {
        mat3<- mat1[,which(timepoints==j)]
        weight.mean<- as.numeric(mat3)}}  
  weight.mean
}
rownames (wm)<- ids
colnames (wm)<- c("E11.1", "E11.2", "E11.3", "E13.1", "E13.2", "E13.3", "E15.1", "E15.2", "E15.3", 
                  "E17.1", "E17.2", "E17.3")

#on (pro-IPC) and off (non-IPC) gene lists from Aprea et al. 2013
#use the weighted means wm file
library(dplyr)
#Rename column 1 of wm files as "Gene.Name"
colnames (on_genes)<-"Gene.Name"
on_data_wm<-inner_join(on_genes, wm, by="Gene.Name")
#add BFs from Table S2
on_data_wm_bf<- inner_join(on_data_wm, scores, by="Gene.Name")
on_dynamic<-on_data_wm_bf[which(exp(on_data_wm_bf$summed.ln(BF))>=3),] # dynamic and partially dynamic genes 
z_on<- t(apply (on_dynamic[,2:13], 1, function(x)(x-mean(x))/sd(x)))
rownames (z_on)<-on_dynamic$Gene.Name

###generate heatmap with clustering
library(gplots)
library(RColorBrewer)
heatmap.2 (as.matrix (z_on),
           scale="none",trace= "none", dendrogram="row" , 
           col=rev(brewer.pal(11,"RdYlBu")), 
           Colv=F,Rowv=T,labCol= c("E11.1", "E11.2", "E11.3", "E13.1", "E13.2", "E13.3", 
                "E15.1", "E15.2", "E15.3", "E17.1", "E17.2", "E17.3"), cexRow=0.75, keysize=1)
z_on_heatmap.2<- heatmap.2 (as.matrix (z_on),
           scale="none",trace= "none", dendrogram="row" , 
           col=rev(brewer.pal(11,"RdYlBu")), 
           Colv=F,Rowv=T,labCol= c("E11.1", "E11.2", "E11.3", "E13.1", "E13.2", "E13.3", 
          "E15.1", "E15.2", "E15.3", "E17.1", "E17.2", "E17.3"), cexRow=0.75, keysize=1)
genes_z_on<-data.frame(rownames(z_on)[z_on_heatmap.2$rowInd]) #gene names in order

##Repeated for off (non-IPC) genes.


