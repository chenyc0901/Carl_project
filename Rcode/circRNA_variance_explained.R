library(edgeR)
library(tidyverse)
library(variancePartition)
rm(list=ls())
#load file
lib_mtx <- read.csv("./rawdata/library_info.csv", row.names = 1)
gene_mtx <- read.csv("./rawdata/gene_count_matrix.csv", row.names = 1, check.names=FALSE)
gene_mtx <- gene_mtx[,rownames(lib_mtx)]
bsj_mtx <- read.csv("./rawdata/circRNA_bsj_noGTEX.csv", row.names = 1, check.names=FALSE)
#convert readcount to cpm by edgeR
circ_DGE <- DGEList(counts = bsj_mtx)
cpm_circRNA <- cpm(circ_DGE) 
#create formula
form <- ~age+(1|gender)+(1|ethnicity_self_identify)+(1|vital_status)+bmi
#fit model
varPart <- fitExtractVarPartModel(cpm_circRNA,form,info)
vp <- sortCols(varPart)
#plotPercentBars(vp[1:10,])
#To calculate the median value of percent
medianpercent <- function(x){return(median(x)*100)}
text <- vp %>% apply(.,2,medianpercent) %>% as.data.frame() %>% rownames_to_column("id")
colnames(text) <-c("id" ,"value")
#draw plot
plotVarPart(vp)+
  geom_text(data=text,aes(y=50,label=(paste0(round(value,digit=2),"%")),x=id),color='blue')

