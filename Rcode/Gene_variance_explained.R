library(tidyverse)
library(cqn)

rm(list=ls())
#load file
  gene_info <- read_csv("rawdata/geneLength_GCcontent.csv")
  gene_count <- read_csv("rawdata/Gene_count.csv")
 
#merge data (matrix without identical id and all zero values)
#merged with information data and filter out all zero values
  gene_selected <-left_join(gene_info,gene_count,by="id") %>% 
                  as_tibble() %>% 
                  rowwise() %>% 
                  mutate(sum_=sum(c_across(starts_with("C")))) %>% 
                  filter(sum_!=0) %>%
                  dplyr::select(-sum_) 
#create gene matrix based on our gene selected (id to rowname)
  gene_matrix <- gene_selected %>% 
                 column_to_rownames(.,var="id") %>% dplyr::select(-c(1:9)) 
#create information of gene (the gene order and id must identical to gene matrix) and convert GC to GC percentage
  gene_info_matrix <- gene_selected %>%dplyr::select(c("GC_content","Gene_length")) %>% 
                      mutate(GC_percent=GC_content/100)
#fit cqn 
  cqn_subset <- cqn(gene_matrix,lengths=gene_info_matrix$Gene_length,x=gene_info_matrix$GC_percent,verbose=TRUE)

#cqnplot draw
  par(mfrow=c(1,2))
  cqnplot(cqn_subset, n = 1,xlab="Length")
  cqnplot(cqn_subset, n = 2,xlab="GC Content")
  