#script to generate Dotplots of BMP pathway genes included in the manuscript "A whole-body atlas of BMP signaling activity in an adult sea anemone"
#to generate Seurat object used for the analysis see https://github.com/technau/nv2_Atlas/tree/main

# Setup ----
# load libraries 
library("Seurat") 
library("pals")
library("RColorBrewer")
library("ggplot2")
library("readxl")

#load workspace and data
setwd()
load ("Alldata.Nv2.updated.Robj")
data1 <- Alldata.Nv2

#drop embryonic data 
data1<-SetIdent(Alldata.Nv2,value='lifehistory')
data1<-subset(data1,idents='AdultSubset')
#drop empty embryo clusters
data1$ID.separate<-droplevels(data1$ID.separate) 
#tabulate cells left
x=as.data.frame(table(data1$ID.separate))
data1<-SetIdent(data1,value='ID.separate')
#filter out any clusters with less than 10 cells
data1<-subset(data1,idents=levels(data1)[which(x$Freq>10)]) 
#drop empty embryo clusters
data1$ID.separate<-droplevels(data1$ID.separate) 
data1<-SetIdent(data1,value = 'ID.separate')
levels(data1)
View(as.data.frame(table(data1$ID.separate)))


#load BMP gene list of interest
mygenes = c(
  "NV2.23832",
  "NV2.12730",
  "NV2.7393",
  "NV2.7394",
  "NV2.12016",
  "NV2.122",
  "NV2.14185",
  "NV2.13963",
  "NV2.24312",
  "NV2.12106",
  "NV2.13010",
  "NV2.13851",
  "NV2.11643",
  "NV2.15308",
  "NV2.10628",
  "NV2.2977",
  "NV2.2978",
  "NV2.6489",
  "NV2.6269",
  "NV2.19847"
)
GOI = unique(genes$name[match(mygenes,genes$geneID)])


#plot expression data
plot=
  DotPlot(data1,'RNA',
  features = unique(GOI),
  scale.by='size',
  col.min = 0,
  col.max = 3,
  group.by='IDs', #group.by="ID.separate"
  scale.max = 100,
  cols = c('lightgrey','#C6237F') )
    +RotatedAxis()
    +coord_flip()
    +theme(legend.position = 'bottom')
    +theme_light()
    +guides(x =  guide_axis(angle = 90)#+geom_point(shape = 21)
  )
plot
