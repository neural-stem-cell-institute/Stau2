#######Figure S3

library(ggplot2)
library(foreach)


RP<-read.delim("ReactomePathways.txt",as.is=T,header=F)  #Table S13a
RP<-RP[grep("Mus",RP[,3]),] #Reactome pathways for mus musculus
rownames(RP)<-RP$V1
#use dynamic, partially dynamic, and stable pathways generated in Fig. 3 code
t1.path.2<- t1.path$data
dynx<-foreach(i=1:223,.combine='rbind') %do% { #223- number of pathways
  x<-t1.path.2[i,]
  x$Cluster<-i-1
  x[grep("REACTOME",x$label),] 
}

####remove the special characters from the label
x<-strsplit(dynx$label,"REACTOME_")
x<-unlist(x)[seq(2,length(unlist(x)),2)]
x<-gsub("_", " ", x)
dynx$label1<-x #add new column label1 without characters

###Edit label names in RP df to match label1 in dynx df
RP$caps<-toupper(RP$V2)
RP$caps<-gsub("\\(","",RP$caps)
RP$caps<-gsub("\\)","",RP$caps)
RP$caps<-gsub("-"," ",RP$caps)
RP$caps<-gsub(",","",RP$caps)
RP$caps<-gsub("\\+","",RP$caps)
RP$caps<-gsub(":"," ",RP$caps)
RP$caps<-gsub("\\/"," ",RP$caps)
RP$caps<-trimws(RP$caps)

dynx$PathID<-foreach(i=1:length(dynx$label1),.combine="c") %dopar% {
  x<-RP[which(RP$caps==dynx$label1[i]),1]
  if(length(x)==0) {x<-NA} else{x}
}
########file with parent Reactome pathway identifiers- Table S13b
RPH<-read.delim("ReactomePathwaysRelation.txt",as.is=T,header=F)
RPH<-RPH[grep("MMU",RPH[,1]),]

#subtract specific parent pathways (top tier) to get second tier pathway branches as needed#subtract specific parent pathways (top tier) to get second tier pathway branches as needed
RPH_subset<-RPH[-(grep("R-MMU-1640170", RPH[,1])),] #cell cycle
RPH_subset<-RPH_subset[-(grep("R-MMU-1266738", RPH_subset[,1])),] #developmental bio
RPH_subset<-RPH_subset[-(grep("R-MMU-74160", RPH_subset[,1])),] #Gene expression
RPH_subset<-RPH_subset[-(grep("R-MMU-1430728", RPH_subset[,1])),]#metabolism
RPH_subset<-RPH_subset[-(grep("R-MMU-162582", RPH_subset[,1])),] #signal transduction
#######column V1 in RPH_subset are higher tier pathways, and V2 are the branches or lower tier 

#######parent pathways identifiers for dynamic genes 
x<-foreach(i=1:length(dynx$PathID),.combine='c') %do% {
  y<-RPH[which(RPH[,2]==dynx$PathID[i]),1]
  if(length(y)==0) {y<-dynx$PathID[i]} else {y[1]}
}
repeat{
  z1<-length(unique(x))
  x<-foreach(i=1:length(x),.combine='c') %do% {
    y<-RPH_subset[which(RPH_subset[,2]==x[i]),1]
    if(length(y)==0) {y<-x[i]} else {y[1]}
  }
  if(length(unique(x))==z1) {break}
}
dynx$rootPath<-x #rootPath- parent pathway identifer
dynx$rootName<-RP[x,]$V2 #rootName- parent pathway name
table(dynx$rootName) #table showing parent pathway names and number of second tier pathways (sub-pathways)

#Repeat for stable pathways
s1.path.2<- s1.path$data
stabx<-foreach(i=1:605,.combine='rbind') %do% { #605- number of pathways
  x<-s1.path.2[i,]
  x$Cluster<-i-1
  x[grep("REACTOME",x$label),] 
}

#####remove the special characters from the label
x<-strsplit(stabx$label,"REACTOME_")
x<-unlist(x)[seq(2,length(unlist(x)),2)]
x<-gsub("_", " ", x)
stabx$label1<-x #add new column label1 without characters

stabx$PathID<-foreach(i=1:length(stabx$label1),.combine="c") %dopar% {
  x<-RP[which(RP$caps==stabx$label1[i]),1]
  if(length(x)==0) {x<-NA} else{x}
}
x<-foreach(i=1:length(stabx$PathID),.combine='c') %do% {
  y<-RPH[which(RPH[,2]==stabx$PathID[i]),1]
  if(length(y)==0) {y<-stabx$PathID[i]} else {y[1]}
}
repeat{
  z1<-length(unique(x))
  x<-foreach(i=1:length(x),.combine='c') %do% {
    y<-RPH_subset[which(RPH_subset[,2]==x[i]),1]
    if(length(y)==0) {y<-x[i]} else {y[1]}
  }
  if(length(unique(x))==z1) {break}
}
stabx$rootPath<-x #rootPath- parent pathway identifer
stabx$rootName<-RP[x,]$V2 #rootName- parent pathway name
table(stabx$rootName) #table showing parent pathway names and number of second tier pathways (sub-pathways)
#repeat for partially dynamic
p1.path.2<- p1.path$data
par.dynx<-foreach(i=1:629,.combine='rbind') %do% { #629- number of pathways
  x<-p1.path.2[i,]
  x$Cluster<-i-1
  x[grep("REACTOME",x$label),] 
}

#####remove the special characters from the label
x<-strsplit(par.dynx$label,"REACTOME_")
x<-unlist(x)[seq(2,length(unlist(x)),2)]
x<-gsub("_", " ", x)
par.dynx$label1<-x #add new column label1 without characters

par.dynx$PathID<-foreach(i=1:length(par.dynx$label1),.combine="c") %dopar% {
  x<-RP[which(RP$caps==par.dynx$label1[i]),1]
  if(length(x)==0) {x<-NA} else{x}
}
x<-foreach(i=1:length(par.dynx$PathID),.combine='c') %do% {
  y<-RPH[which(RPH[,2]==par.dynx$PathID[i]),1]
  if(length(y)==0) {y<-par.dynx$PathID[i]} else {y[1]}
}
repeat{
  z1<-length(unique(x))
  x<-foreach(i=1:length(x),.combine='c') %do% {
    y<-RPH_subset[which(RPH_subset[,2]==x[i]),1]
    if(length(y)==0) {y<-x[i]} else {y[1]}
  }
  if(length(unique(x))==z1) {break}
}
par.dynx$rootPath<-x #rootPath- parent pathway identifer
par.dynx$rootName<-RP[x,]$V2 #rootName- parent pathway name
table(par.dynx$rootName) #table showing parent pathway names and number of second tier pathways (sub-pathways)


#######plotting
#######make and read in a file with pathways for all 3 categories (Table S7d)
parent_pathways$Pathway<- factor (parent_pathways$Pathway, levels=unique(parent_pathways$Pathway)) 

p<-ggplot(parent_pathways, aes(y=factor(Stability), x=Pathway, size=percent, color=Stability))
p<-p + geom_point()+coord_flip()
p<-p+theme(axis.text.x = element_text(size=8, angle=90, 
                                      hjust=1, vjust=0.5))
p
p
