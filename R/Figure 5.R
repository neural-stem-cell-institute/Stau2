########Figure 5 and Fig. S5 code
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

#Figure 5A
#VZ, SVZ, VZ-SVZ, CP gene lists from Fietz et al. 2012
# combine vz_svz and svz genes- this gives SVZ genelist
#generate weighted means file with BFs column same as for Fig. 4



#to find STAU2 cargo SVZ DB genes
library (GOfuncR)
svz_DB<- get_anno_genes(go_ids="GO:0006351",
                        genes=SVZ$Gene.Name, database='Mus.musculus')
colnames (svz_DB)[2]<-"Gene.Name"

svz_DB_wm<- inner_join(svz_DB, masterpeaks_wm, 
                       by="Gene.Name") #generate file columns for expression values/counts at timepoints


svz_DB_wm$E11_average<- rowMeans(svz_DB_wm[3:5],)#mean E11 values
#repeat for E13, E15, E17

library(superheat)
library(RColorBrewer)
#use row means columns
z_SVZ<- data.frame (t(apply (svz_DB_scores[,15:18], 1, function(x)(x-mean(x))/sd(x))), 
                    row.names = svz_DB_scores$Gene.Name) 
z_SVZ<-setDT(z_SVZ, keep.rownames = "Gene.Name")[]
z_SVZ_bf<- data.frame (inner_join(z_SVZ, scores, by="Gene.Name"))
rownames(z_SVZ_bf)<- z_SVZ_bf$Gene.Name
z_SVZ_bf<- z_SVZ_bf[which(z_SVZ_bf$Summed.BF>=1.163),] #dynamic+partially dynamic genes

superheat.z.SVZ.means.bf.clustered<- superheat(z_SVZ_bf[,2:5], n.clusters.rows = 2, 
clustering.method="hierarchical",left.label.size=0.2, left.label.text.size = 2.5, bottom.label.text.size = 3, 
   pretty.order.rows = TRUE,left.label = "variable", heat.pal = rev(brewer.pal(11,"RdYlBu")), 
   yr = z_SVZ_bf$Summed.BF,  yr.plot.type = "bar", yr.axis.name = "ln(BF)")

write.csv(superheat.z.SVZ.means.bf.clustered$membership.rows, "genes_z_SVZ_superheat.csv")

#Repeat for CP genes in Fig. S5




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
library(foreach)
m<-foreach(i=1:2) %do% {
  rownames(z)[which(z[,i]>0.9)]
}


l2.3.4<-setdiff(m[[1]],m[[2]]) #values specific to l2.3.4
l5.6<-setdiff(m[[2]],m[[1]]) #values specific to l5, l6, l6b
lvec<-c(l2.3.4,l5.6)


y<-y[lvec,]
y<-as.matrix(y) #matrix containing selected layer marker genes

mydata2<- mydata[,9:20] #columns with counts/expression values for the timepoints
z<- t(apply (mydata2, 1, function(x)(x-mean(x))/sd(x))) #zscores,mean=0, sd=above/below 0
rownames (z)<-mydata$gene_symbol

yy<-y[intersect(rownames(z),rownames(y)),] #subset out the layer marker genes (y) from my data

gene.name<- data.frame (rownames (z)) 
s<-cbind(gene.name, z) #add gene names to z
colnames (s) [1]<-"Gene.Name" 

zz<-foreach(i=1:length(rownames(yy)),.combine='rbind') %do% {
  p<-s[rownames(s)==rownames(yy)[i],]
  p<-as.matrix(p[,2:13])
  if(dim(p)[1]>1) {p<-apply(p,2,sum)} #genes with multiple peaks have dim>1
  p
} 

#apply weights from layer markers file to layer marker genes in my data
yszz<-foreach(i=1:dim(zz)[2],.combine='rbind') %do% {
  z1<-zz[,i]
  z2<-foreach(b=1:2,.combine='c') %dopar% {
    wghts<-yy[intersect(rownames(yy),m[[b]]),b]
    z3<-sum(z1[intersect(rownames(yy),m[[b]])]*wghts, na.rm = T)/sum(wghts)
  }
}
colnames(yszz)<-c("L2/3/4","L5/6")
rownames(yszz)<-colnames(zz)

# for plotting
  library(reshape2)
yszz_df<-data.frame(yszz)
yszz_plot<- melt(yszz_df)
yszz_plot$timepoint<- c(11,11,11,13,13,13, 15,15,15,17,17,17)
colnames (yszz_plot)[1]<- "Layer"

library(ggplot2)
p<-ggplot(aes(timepoint,value),data=yszz_plot)
p<-p + geom_smooth(aes(color=Layer),size=2,fill="gray30", level=0.9)
p<-p + theme_light() + scale_color_manual(values=c("orange", "purple"))
p

#####5E: See Fig. S3 code, gene lists for DL cargo and not cargo from Tables S11a,b

REACTOME <- msigdb_gsets(species="Mus musculus", category="C2", subcategory="CP:REACTOME")
dl.data.path<-hypeR(l5.6_data, REACTOME, background=55487, fdr=0.01)
dl.notdata.path<-hypeR(l5.6_notdata, REACTOME, background=55487, fdr=0.01)


dl.data.path.2<- dl.data.path$data
dl.notdata.path.2<- dl.notdata.path$data
ul.data.path.2<- ul.data.path$data
ul.notdata.path.2<- ul.notdata.path$data


RP<-read.delim("ReactomePathways.txt",as.is=T,header=F)  #Table S13a
RP<-RP[grep("Mus",RP[,3]),] #Reactome pathways for mus musculus
rownames(RP)<-RP$V1
#use DL_cargo, not_cargo, UL cargo, not_cargo pathways
dl.data.pathx<-foreach(i=1:85,.combine='rbind') %do% { #85- number of pathways
  x<-dl.data.path.2[i,]
  x$Cluster<-i-1
  x[grep("REACTOME",x$label),] 
}

#remove the special characters from the label
x<-strsplit(dl.data.pathx$label,"REACTOME_")
x<-unlist(x)[seq(2,length(unlist(x)),2)]
x<-gsub("_", " ", x)
dl.data.pathx$label1<-x #add new column label1 without characters

#Edit label names in RP df to match label1 in dl.data.pathx df
RP$caps<-toupper(RP$V2)
RP$caps<-gsub("\\(","",RP$caps)
RP$caps<-gsub("\\)","",RP$caps)
RP$caps<-gsub("-"," ",RP$caps)
RP$caps<-gsub(",","",RP$caps)
RP$caps<-gsub("\\+","",RP$caps)
RP$caps<-gsub(":"," ",RP$caps)
RP$caps<-gsub("\\/"," ",RP$caps)
RP$caps<-trimws(RP$caps)

dl.data.pathx$PathID<-foreach(i=1:length(dl.data.pathx$label1),.combine="c") %dopar% {
  x<-RP[which(RP$caps==dl.data.pathx$label1[i]),1]
  if(length(x)==0) {x<-NA} else{x}
}
#file with parent Reactome pathway identifiers- Table S13b
RPH<-read.delim("ReactomePathwaysRelation.txt",as.is=T,header=F)
RPH<-RPH[grep("MMU",RPH[,1]),]

#subtract specific parent pathways (top tier) to get second tier pathway branches as needed#subtract specific parent pathways (top tier) to get second tier pathway branches as needed
RPH_subset<-RPH[-(grep("R-MMU-1640170", RPH[,1])),] #cell cycle
RPH_subset<-RPH_subset[-(grep("R-MMU-1266738", RPH_subset[,1])),] #developmental bio
RPH_subset<-RPH_subset[-(grep("R-MMU-74160", RPH_subset[,1])),] #Gene expression
RPH_subset<-RPH_subset[-(grep("R-MMU-1430728", RPH_subset[,1])),]#metabolism
RPH_subset<-RPH_subset[-(grep("R-MMU-162582", RPH_subset[,1])),] #signal transduction
#column V1 in RPH_subset are higher tier pathways, and V2 are the branches or lower tier 

#parent pathways identifiers for dynamic genes 
x<-foreach(i=1:length(dl.data.pathx$PathID),.combine='c') %do% {
  y<-RPH[which(RPH[,2]==dl.data.pathx$PathID[i]),1]
  if(length(y)==0) {y<-dl.data.pathx$PathID[i]} else {y[1]}
}
repeat{
  z1<-length(unique(x))
  x<-foreach(i=1:length(x),.combine='c') %do% {
    y<-RPH_subset[which(RPH_subset[,2]==x[i]),1]
    if(length(y)==0) {y<-x[i]} else {y[1]}
  }
  if(length(unique(x))==z1) {break}
}
dl.data.pathx$rootPath<-x #rootPath- parent pathway identifer
dl.data.pathx$rootName<-RP[x,]$V2 #rootName- parent pathway name
table(dl.data.pathx$rootName) #table showing parent pathway names and number of second tier pathways (sub-pathways)

#Repeat for DL not cargo pathways
dl.notdata.pathx<-foreach(i=1:59,.combine='rbind') %do% { #59- number of pathways
  x<-dl.notdata.path.2[i,]
  x$Cluster<-i-1
  x[grep("REACTOME",x$label),] 
}

#remove the special characters from the label
x<-strsplit(dl.notdata.pathx$label,"REACTOME_")
x<-unlist(x)[seq(2,length(unlist(x)),2)]
x<-gsub("_", " ", x)
dl.notdata.pathx$label1<-x #add new column label1 without characters

dl.notdata.pathx$PathID<-foreach(i=1:length(dl.notdata.pathx$label1),.combine="c") %dopar% {
  x<-RP[which(RP$caps==dl.notdata.pathx$label1[i]),1]
  if(length(x)==0) {x<-NA} else{x}
}
x<-foreach(i=1:length(dl.notdata.pathx$PathID),.combine='c') %do% {
  y<-RPH[which(RPH[,2]==dl.notdata.pathx$PathID[i]),1]
  if(length(y)==0) {y<-dl.notdata.pathx$PathID[i]} else {y[1]}
}
repeat{
  z1<-length(unique(x))
  x<-foreach(i=1:length(x),.combine='c') %do% {
    y<-RPH_subset[which(RPH_subset[,2]==x[i]),1]
    if(length(y)==0) {y<-x[i]} else {y[1]}
  }
  if(length(unique(x))==z1) {break}
}
dl.notdata.pathx$rootPath<-x #rootPath- parent pathway identifer
dl.notdata.pathx$rootName<-RP[x,]$V2 #rootName- parent pathway name
table(dl.notdata.pathx$rootName) #table showing parent pathway names and number of second tier pathways (sub-pathways)
#make and read in a file with parent pathways for DL cargo and not cargo genes (Table S11e)
dl_parent_pathways<-read.csv("DL_parent_paths_plotting.csv")
dl_parent_pathways$Pathway<- factor (dl_parent_pathways$Pathway, levels=unique(dl_parent_pathways$Pathway)) 

p<-ggplot(dl_parent_pathways[1:39,], aes(y=factor(Category), x=Pathway, size=percent.of.all.pathways, color=Category))
p<-p + geom_point()+coord_flip()
p<-p+theme(axis.text.x = element_text(size=8, angle=90, 
                                      hjust=1, vjust=0.5))
p

#S5 C,D
ul.data.path<-hypeR(l2.3.4_data, REACTOME, background=55487, fdr=0.01)
ul.notdata.path<-hypeR(l2.3.4_notdata, REACTOME, background=55487, fdr=0.01)
#See function hypeR, figure 2 code
hyp_plot2(ul.data.path,main="Top Signaling Pathways for UL cargo Genes")
hyp_plot2(ul.notdata.path,main="Top Signaling Pathways for UL not cargo Genes")
