#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-r", "--deglists", 
                    type="character", 
                    help="comma separated DESeq2 DEG lists",
                    required=TRUE)
parser$add_argument("-s", "--samplenames", 
                    type="character", 
                    help="samplenames as comma separated lists",
                    required=TRUE)
parser$add_argument("-i", "--indexcols", 
                    type="character",
                    help="comma separated list columns to together use as index",
                    required=TRUE)
parser$add_argument("-c", "--includecols", 
                    type="character",
                    help="columns to include (will be prefixed with samplenames)",
                    required=TRUE)
parser$add_argument("-o", "--outfile", 
                    type="character",
                    help="outfile",
                    required=TRUE)

args <- parser$parse_args()

suppressPackageStartupMessages(library("tidyverse"))

debug=0

deglists=unlist(strsplit(args$deglists,","))
samplenames=unlist(strsplit(args$samplenames,","))
indexcols=unlist(strsplit(args$indexcols,","))
includecols=unlist(strsplit(args$includecols,","))

if (debug==1){
  deglists="/Volumes/Wolin/mESC_slam_analysis/find_mutation_101521/rsem/DEG.KO_slam_vs_WT_slam.tsv,/Volumes/Wolin/mESC_slam_analysis/find_mutation_101521/rsem/DEG.KO_vs_WT.tsv"
  samplenames="KO_slam_vs_WT_slam,KO_vs_WT"
  indexcols="gene_id,gene_name"
  includecols="log2FoldChange,FDR"
  deglists=unlist(strsplit(deglists,","))
  samplenames=unlist(strsplit(samplenames,","))
  indexcols=unlist(strsplit(indexcols,","))
  includecols=unlist(strsplit(includecols,","))  
}

read_deg<-function(deg,sn,includecols,indexcols){
  degdf=read.csv(deg,header = TRUE,sep = "\t")
  filtercols=c(indexcols,includecols)
  degdf %>% dplyr::select(all_of(filtercols)) %>%
    unite("geneID",all_of(indexcols),sep="##",remove=TRUE) %>%
    column_to_rownames(var="geneID") -> degdf
  oldcols=colnames(degdf)
  newcols=c()
  for (j in 1:length(oldcols)){
    newsn=paste(sn,oldcols[j],sep="_")
    newcols=c(newcols,newsn)
  }
  colnames(degdf)=newcols
  degdf <- rownames_to_column(degdf,var="geneID")
  return(degdf)
}

all_degs=read_deg(deglists[1],samplenames[1],includecols,indexcols)
for (i in 2:length(deglists)){
  newdegs=read_deg(deglists[i],samplenames[i],includecols,indexcols)
  all_degs=merge(all_degs,newdegs,by=c("geneID"))
}

all_degs <- separate(all_degs,col="geneID",into=c("gene_id","gene_name"),sep = "##")

write.table(all_degs,file=args$outfile,row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
