# This script uses some simple code to de-identify Seurat object and add metadata (like cell annotations). The required inputs are - Seurat obj, table with sample/orig.ident (e.g. patientID1,patientID2,etc.) and the ID to change it to (e.g. U1,U2,etc.) in columns named 'SampleID' and 'Updated_ID', respectively. Addiitionally a metadata spreadsheet is needed that contains seurat cluster numbers and their matching annotations.

# Author: Surbhi Sona (Ting Lab)

library(Seurat)
library(dplyr)
library(data.table)
library(openxlsx)
library(readxl)
library(tidyverse)

# Read original Seurat obj
obj_g<-readRDS('/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/2020_09_21_ureter10/2020_09_21_ureter10_clustered.rds')
obj_str<-readRDS('/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/Stromal subset/2021_02_19_ureter10_stromal/2021_02_19_ureter10_stromal_clustered.rds')
obj_uro<-readRDS('/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/Uro subset/2021_04_01_ureter10_uro_PC50_res0.2/2021_04_01_ureter10_uro_PC50_res0.2_clustered.rds')
obj_imm<-readRDS('/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/Immune subset/2021_02_25_ureter10_immune/2021_02_25_ureter10_immune_clustered.rds')
obj_trio<-readRDS('/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/2020_09_22_ureter3_PC30_res0.3/2020_09_22_ureter3_PC30_res0.3_clustered.rds')


# Set variables
file<-'/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/SampleID_Key.xlsx'
sample_tab<-read_xlsx(file)

obj_out<-'/home/sonas/ucsc/'

# New filenames
prefix_global<-'adult-ureter'
prefix_uro<-'adult-ureter-uro-subset'
prefix_str<-'adult-ureter-stromal-subset'
prefix_imm<-'adult-ureter-immune-subset'
prefix_trio<-'adult-ureter-trio'


# Function to de-identify Seurat object

de_identify_obj<-function(obj,sample_tab,out,fname,save=TRUE)
{
  indx<-match(obj$orig.ident,sample_tab$Sample)
  obj$orig.ident<-sample_tab$Updated_ID[indx]
  indx2<-match(obj$sample,sample_tab$Sample)
  obj$sample<-sample_tab$Updated_ID[indx2]
  if(save)
  {
    saveRDS(obj,file=paste0(out,fname))
  }else{
    return(obj)
  }
}


# function to add cell annotations to Seurat obj metadata

# Update Metadata

update_metacol<-function(obj,tab_path,tab_name,out,fprefix,save=TRUE)
{
  #Remove redundant seurat cluster column
  col<-grep('integrated',colnames(obj@meta.data))
  obj@meta.data[,col]<-NULL
  
  # Read table
  tab<-read_xlsx(paste0(tab_path,tab_name))
  
  # Add cell type info
  indx<- match(obj$seurat_clusters,tab$clusters)
  cell_type<-tab$cell_type[indx]
  obj<-AddMetaData(obj,metadata = cell_type,col.name = 'cell_type')
  
  if(save==TRUE)
  {
    saveRDS(obj,paste0(out,fprefix,'.rds'))
  }else
  {
    return(obj)
  }
  
  
}



## Deidentify 

de_identify_obj(obj=obj_g,sample_tab = sample_tab,out=obj_out,fname=paste0(prefix_global,'.rds'),save = TRUE)
de_identify_obj(obj=obj_uro,sample_tab = sample_tab,obj_out,fname=paste0(prefix_uro,'.rds'),save = TRUE)
de_identify_obj(obj=obj_str,sample_tab = sample_tab,obj_out,fname=paste0(prefix_str,'.rds'),save = TRUE)
de_identify_obj(obj=obj_imm,sample_tab = sample_tab,obj_out,fname=paste0(prefix_imm,'.rds'),save = TRUE)
de_identify_obj(obj=obj_trio,sample_tab = sample_tab,obj_out,fname=paste0(prefix_trio,'.rds'),save = TRUE)



## Add metadata
path<-'/home/sonas/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/UCSC/cluster_annotation_tables/'

# These files should contain 2 columns: 
#'clusters' that denote seurat cluster (e.g. 0,1,2) and 'cell_type' which contains matching cell annotations for each cluster (e.g. basal cells, umbrella cells, intetermidiate cells)

global<-'global.xlsx'
immune<-'immune.xlsx'
str<-'stromal.xlsx'
uro<-'uro.xlsx'
trio<-'trio.xlsx'

# Read in de-identified obj

obj_g2<-readRDS(paste0(obj_out,prefix_global,'.rds'))
obj_uro2<-readRDS(paste0(obj_out,prefix_uro,'.rds'))
obj_str2<-readRDS(paste0(obj_out,prefix_str,'.rds'))
obj_imm2<-readRDS(paste0(obj_out,prefix_imm,'.rds'))
obj_trio2<-readRDS(paste0(obj_out, prefix_trio,'.rds'))

# add cell annotation
update_metacol(obj=obj_g2,tab_path = path,tab_name=global,out=obj_out,fprefix = prefix_global, save=TRUE)
update_metacol(obj=obj_str2,tab_path = path,tab_name=str,out=obj_out,fprefix = prefix_str, save=TRUE)
update_metacol(obj=obj_imm2,tab_path = path,tab_name=immune,out=obj_out,fprefix = prefix_imm, save=TRUE)
update_metacol(obj = obj_uro2,tab_path = path,tab_name=uro,out=obj_out,fprefix = prefix_uro, save=TRUE)
update_metacol(obj = obj_trio2,tab_path = path,tab_name=trio,out=obj_out,fprefix = prefix_trio, save=TRUE)


