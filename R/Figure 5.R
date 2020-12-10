#VZ, SVZ, VZ-SVZ, CP gene lists from Fietz et al. 2012
#for all svz genes, combine vz_svz and svz genes
#generate weighted means file for svz genes with BFs column same as for Fig. 4


svz_wm<-inner_join(svz, wm_scores, by="Gene.Name")
#to find STAU2 cargo SVZ DB genes
library (GOfuncR)
svz_DB<- get_anno_genes(go_ids="GO:0006351",
                        genes=svz_wm$Gene.Name, database='Mus.musculus')
colnames (SVZ_DB)[2]<-"Gene.Name"

svz_DB_wm<- inner_join(svz_wm_DB, svz_wm, 
                       by="Gene.Name") #generate file columns for expression values/counts at timepoints

svz_DB_wm$E11_average<- rowMeans(svz_DB_wm[3:5],)#mean E11 values
#repeat for E13, E15, E17
svz_DB_wm_dynamic<-svz_DB_wm[which(exp(svz_DB_wm$summed.ln(BF))>=3),] # dynamic and partially dynamic genes 

svz_DB_wm_dyn_means<- svz_DB_wm_dynamic[,16:19]#subset mean expression values/counts columns and BFs

library(gplots)
library(RColorBrewer)
svz_heatmap<- t(apply (svz_DB_wm_dyn_means[1:4], 1, function(x)(x-mean(x))/sd(x))) #z-scores for E11, E13, E15, E17 mean expression values
rownames(svz_heatmap)<- svz_DB_wm_dynamic$Gene.Name
heatmap.2 (as.matrix (svz_heatmap),
           scale="row",trace= "none", dendrogram="row" , 
           col=rev(brewer.pal(11,"RdYlBu")), 
           Colv=F,Rowv=T,labCol= c("E11", "E13", "E15", "E17"), 
           cexRow=0.75, keysize=1)
svz_heatmap.2<-heatmap.2 (as.matrix (svz_heatmap),
                          scale="row",trace= "none", dendrogram="row" , 
                          col=rev(brewer.pal(11,"RdYlBu")), 
                          Colv=F,Rowv=T,labCol= c("E11", "E13", "E15", "E17"), 
                          cexRow=0.75, keysize=1)
genes_svz<-data.frame(rownames(svz_heatmap)[svz_heatmap.2$rowInd])

####Repeated for CP genes


#Figure 5D- layer markers

#download Belgard et al., 2011 (Table S2)- Layer Enrichment Probabilities of Known Genes

mydata<- data.frame (read.csv("MasterPeaks_trimmed.csv", as.is=T)) #Table S1
layermarkers<- data.frame (read.csv("adult_brain_markers_belgard.csv", as.is=T))
library(dplyr)
x<-layermarkers[-(which(layermarkers$X.1=="")),] #remove blank values
x<-x[-(which(duplicated(x$X.1))),] #remove duplicated values
x$l234_average<- rowMeans(x[10:11],) #average upper (L2,3,4) probabilties 
x$l56_average<- rowMeans(x[12:14],) #average deep (L5, 6, 6b) probabilties
y<-x[,16:17] #subset columns with avg l2/3/4 and average l5/6/6b probabilities
rownames(y)<-x$X.1


p<-ecdf(unlist(y)) 
z<-t(apply(y,1,function(x) p(x))) #apply ecdf function to normalize across rows
m<-apply(z,1,max) 
z<-z[which(m>.9),] #select genes with probabilities>0.9 for each layer group
colnames(z)<-colnames(y)
m<-foreach(i=1:2) %do% {
  rownames(z)[which(z[,i]>0.9)]
}


l2.3.4<-setdiff(m[[1]],m[[2]]) #values specific to l2.3.4
l5.6<-setdiff(m[[2]],m[[1]]) #values specific to l5, l6, l6b
lvec<-c(l2.3.4,l5.6)


y<-y[lvec,]
y<-as.matrix(y) #matrix containing selected layer marker genes

mydata2<- mydata[,21:32] #columns with counts/expression values for the timepoints
z<- t(apply (mydata2, 1, function(x)(x-mean(x))/sd(x))) #zscores,mean=0, sd=above/below 0
rownames (z)<-mydata$Gene.Name

yy<-y[intersect(rownames(z),rownames(y)),] #subset out the layer marker genes (y) from my data

gene.name<- data.frame (rownames (z)) 
s<-cbind(gene.name, z) #add gene names to z
colnames (s) [1]<-"Gene.Name" 

zz<-foreach(i=1:length(rownames(yy)),.combine='rbind') %do% {
  p<-s[s$Gene.Name==rownames(yy)[i],]
  p<-as.matrix(p[,2:13])
  if(dim(p)[1]>1) {p<-apply(p,2,sum)} #genes with multiple peaks have dim>1
  p
} 


rownames(zz)<-rownames(yy)


#apply weights from layer markers file to layer marker genes in my data
yszz<-foreach(i=1:dim(zz)[2],.combine='rbind') %do% {
  z1<-zz[,i]
  z2<-foreach(b=1:2,.combine='c') %dopar% {
    wghts<-yy[intersect(rownames(yy),m[[b]]),b]
    z3<-sum(z1[intersect(rownames(yy),m[[b]])]*wghts)/sum(wghts)
  }
}
colnames(yszz)<-c("L2/3/4","L5/6")
rownames(yszz)<-colnames(zz)

# for plotting
yszz_df<-data.frame(yszz)
yszz_df$timepoint<- c(11,11,11,13,13,13, 15,15,15,17,17,17)

library(reshape2)
ysz_plot_2<- melt(yszz_df) 

library(ggplot2)
p<-ggplot(aes(timepoint,value),data=yszz_plot_2)
p<-p + geom_smooth(aes(color=Layer),size=2,fill="gray30")
p<-p + theme_light() + scale_color_manual(values=c("orange", "purple"))
p

#####5E: See Fig. S3 code, gene lists for DL cargo and not cargo from Tables S11a,b
