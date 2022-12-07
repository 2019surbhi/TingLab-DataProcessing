#!/usr/bin/env Rscript --vanilla

#Loading dependencies for working in command line and HPC environment
library(argparser)
library(parallel) 
library(tictoc)

#Loading dependencies for data wrangling and processing
library(dplyr) 
library(data.table) 

#'omics' packages
library(Seurat)
library(Signac)
library(GenomicRanges)
library(patchwork)
library(microbenchmark)

#Visualizations packages
library(ggplot2)
library(lattice)    
library(grid) 
library(gridExtra)

source('/bme/home/surbhis/Scripts/R_Scripts_SS/sc_atac_functions.R')
##### Creating options to add command line arguments #####

parser<-arg_parser(name="sc_atac_merge.R",description="Processing multiple atac seq data using Seurat pipeline")

parser<-add_argument(
  parser,
  arg="--input_dir",
  short = '-i',
  type="character",
  default="./",
  help="Location+name of directory containing matrix count data.Must have trailing '/'. Default: [./]"
)

parser<-add_argument(
  parser,
  arg="--output_dir",
  short = '-o',
  type="character",
  default="./data/cluster_plots/",
  help="location+name for output directory. Must contain trailing '/'. Default: [./data/cluster_plots/]"
)

parser<-add_argument(
  parser,
  arg="--samples",
  short = '-s',
  type="character",
  default="",
  help="Name of sample to be processed. Must be name of sample directory. For multiple samples, split by comma (,)"
)

parser<-add_argument(
  parser,
  arg="--run_name",
  short = '-r',
  type="character",
  default="",
  help="Title used to identify all plots that correspond to the current run instance. Default: None"
)

parser<-add_argument(
  parser,
  arg="--matrix_name",
  short = '-n',
  type="character",
  default="filtered_peak_bc_matrix",
  help="Name of the matrix file/directory with no file extension. Eg: 'raw_feature_bc_matrix'. Default: [filtered_feature_bc_matrix]"
)

parser<-add_argument(
  parser,
  arg="--annotation_file",
  short = '-a',
  type="character",
  default="/home/surbhis/ref_files_bc/Homo_sapiens.GRCh37.87.gtf",
  help="Enter the annotation file with path. Default:[gtf file in cellranger ref folder]"
)

parser<-add_argument(
  parser,
  arg="--meta_file",
  short = '-m',
  type="character",
  default="",
  help="Enter the meta data file in csv format along with path"
)

parser<-add_argument(
  parser,
  arg="--rna_file",
  short = '-f',
  type="character",
  default="",
  help="Enter the processed sc-RNA seq file in rds format along with path"
)

parser<-add_argument(
  parser,
  arg="--cores",
  short = '-c',
  type="character",
  default='1',
  help= "Number of cores to use when processing samples in parallel. Use 'all' to use all available cores. Windows MUST use 1. Default: 1"
)

parser<-add_argument(
  parser,
  arg="--verbose",
  short = '-v',
  flag=TRUE,
  help="Flag set for verbose output. Warning: May lead to performance drops."
)

argv <- parse_args(parser)

argv

##### Argument Processing #####

#Determine number of cores
if(length(grep('^all$',argv$cores,ignore.case = TRUE)) != 0){
  argv$cores<-detectCores()
}else{
  argv$cores<-as.numeric(argv$cores)
}

#Validate sample names
if(argv$samples==''){
  if(argv$verbose){
    cat("Sample is empty. Terminating Run",sep='\n')
    tic.clear()
    quit(save = 'no')
  }
  quit()
}else
  {samples=unlist(strsplit(argv$samples, split=','))}

# Prepare individual filenames using run ID
run_file_name<-''
if(argv$run_name!=''){
  run_file_name<- gsub(' ','_',argv$run_name)
  run_file_name<-paste0(run_file_name,'-')
}


#Delete Later !!!
add_function_or_delete<-funtion()
{
    
     #Data pre-processing
    DefaultAssay(sample.atac)<-"ACTIVITY"
    sample.atac<-FindVariableFeatures(sample.atac)
    sample.atac<-NormalizeData(sample.atac)
    sample.atac<-ScaleData(sample.atac)

     #Process peak matrix
    DefaultAssay(sample.atac)<-"peaks"
    VariableFeatures(sample.atac)<-names(which(Matrix::rowSums(sample.atac)>100))
    sample.atac<-RunLSI(sample.atac, n=50, scale.max=NULL)
    sample.atac<-RunUMAP(sample.atac,reduction="lsi",dims=1:50)

    #sample.atac$dataset<-sample_name

    return(sample.atac)
sample.atac.filtered<-subset(sample.atac, subset=prediction.score.max>0.5)
sample.atac.filtered$predicted.id<-factor(sample.atac.filtered$predicted.id, levels=levels(sample.rna))

sample<-list(sample.atac,sample.atac.filtered)
return(sample)

#Filter
sample.atac.filtered<-subset(sample.atac, subset=prediction.score.max>0.5)
sample.atac.filtered$predicted.id<-factor(sample.atac.filtered$predicted.id, levels=levels(sample.rna))

sample<-list(sample.atac,sample.atac.filtered)
return(sample)
}

}


##### Main Script #####

## Create atac seq Seurat object ##
n<-length(samples)
atac.obj.list<-list()

atac.obj.list<-lapply(samples[1:n],create_atac_obj,argv$input_dir,argv$annotation_file)

## Sample QC Plots ##
lapply(atac.obj.list[1:n],get_qc_plots,argv$out_dir)
#atac.obj.list<-lapply(samples[1:n],get_tss_enrichment, samples)


## Merge data ##
# Create common peak set
combined.peaks<-UnifyPeaks(object.list=atac.obj.list, mode="reduce")

# Quantify peaks in each set
atac.obj.list<-lapply(atac.obj.list, quantify_peaks, combined.peaks)

# Merge sc.atac object
fragment.path<-dirname(argv$input_dir)
fragment.path<-paste0(fragment.path,'/fragments.tsv.gz')
atac.merged<-merge_atac_data(atac.obj.list,fragment.path,argv$samples)
atac.merged$tech<-'atac'

## Filter cells ##
get_qc_plots(atac.merged,out_dir=argv$out_dir,sample.name='merged')

atac.merged<-filter_cells(atac.merged, peak_frag=thres[1:2],read_pct=thres[3],blk_ratio=thres[4], nuc_signal=thres[5], tss_enrich=thres[6])

## Data pre-process and clustering ##
atac.merged<-data_pre_process_activity(atac.merged)
atac.merged<-data_pre_process_atac(atac.merged)
atac.merged<-get_clusters(atac.merged)

print_umap(atac.merged, plotTitle=argv$run_name,feature=NULL,oDir=argv$output_dir,PrintLabels=TRUE)

## Add RNA seq data ##

#Read RNA seq file
sample.rna<-readRDS(argv$rna_file)
sample.rna$tech<-"rna"

print_umap(sample.rna, feature=NULL,titleGroup='sc-RNA-seq',plotTitle="RNAseq", PrintLabels=TRUE)

#integrate RNA data


atac.merged<-integrate_rna_atac(atac.merged, sample.rna)

#co-embed RNA and atac data
coembed<-run_coembed_process(sample.atac, sample.rna)

print_umap(coembed, feature="tech")
print_umap(coembed, feature="celltype", PrintLabels=TRUE,set_repel=TRUE)
print_umap(coembed, split="tech", feature="celltype",Printlabels=TRUE,set_repel=TRUE)

cat("End of Script", '\n')

