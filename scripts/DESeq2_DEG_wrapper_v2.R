#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
setwd("~/Documents/GitRepos/rbl4/scripts")

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-r", "--rawcountsmatrix", 
                    type="character", 
                    help="raw count matrix",
                    required=TRUE)
parser$add_argument("-s", "--coldata", 
                    type="character", 
                    help="colData or study design",
                    required=TRUE)
parser$add_argument("-i", "--indexcols", 
                    type="character",
                    help="comma separated list columns to together use as index",
                    required=TRUE)
parser$add_argument("-c", "--condition1", 
                    type="character",
                    help="condition1 or group1",
                    required=TRUE)
parser$add_argument("-d", "--condition2", 
                    type="character",
                    help="condition2 or group2.... condition1 vs condition2 contrast will be made",
                    required=TRUE)
parser$add_argument("-p", "--slam", 
                    type="integer",
                    help="slam=1;no slam=0;unmasked=2",
                    required=TRUE)
parser$add_argument("-m", "--masked", 
                    type="integer",
                    help="masked=1;unmasked=0",
                    required=TRUE)
parser$add_argument("-q", "--mutated", 
                    type="integer",
                    help="mutated=1;unmutated=0;neither=2",
                    required=TRUE)
parser$add_argument("-u", "--spliceaware", 
                    type="integer",
                    help="spliceaware=1;spliceunaware=0",
                    required=TRUE)
parser$add_argument("-t", "--cpmcutoff", 
                    type="character",
                    help="cpm threshold ... if cpm less than this value for all samples, then the gene is excluded",
                    required=FALSE,
                    default=1)
parser$add_argument("-o", "--outdir", 
                    type="character",
                    help="outputdir",
                    required=TRUE)

args <- parser$parse_args()

scriptdir=getwd()

suppressPackageStartupMessages(library("rmarkdown"))

debug=0
# debug=1
if (debug==0){
  condition1=args$condition1
  condition2=args$condition2
  rawcountsmatrix=args$rawcountsmatrix
  coldata=args$coldata
  slam=args$slam
  masked=args$masked
  mutated=args$mutated
  spliceaware=args$spliceaware
  indexcols=args$indexcols
  outdir=args$outdir
  cpmcutoff=args$cpmcutoff
}else{
  condition1="KO"
  condition2="WT"
  rawcountsmatrix="/Volumes/Wolin/mESC_slam_analysis/find_mutation_101521/counts/all_w_genename.genecounts.tsv"
  coldata="/Volumes/Wolin/mESC_slam_analysis/find_mutation_101521/counts/coldata.tsv"
  slam=0
  masked=0
  mutated=2
  spliceaware=1
  cpmcutoff=1
  indexcols="gene_id,gene_name"
  outdir="/Users/kopardevn/Documents/Projects/rbl4/find_mutations_101521/DEGs"
}

code=paste0(slam,masked,mutated,spliceaware)

contrast0=paste(condition1,"vs",condition2,sep="_")
contrast=paste(contrast0,code,sep=".")
countsfn=paste(contrast,"counts","tsv",sep=".")
coldatafn=paste(contrast,"coldata","tsv",sep=".")
degfn=paste(contrast,"DEG","tsv",sep=".")

# create counts file
# print(slam)
# print(masked)
# print(mutated)
# print(spliceaware)


rawcounts=read.csv(rawcountsmatrix,header=TRUE,sep="\t",
    check.names = FALSE,comment.char = "#",strip.white = TRUE)
rawcoldata=read.csv(coldata,header = TRUE,sep = "\t",
    check.names = FALSE,comment.char = "#",strip.white = TRUE)

k_slam=rawcoldata$slam==slam
table(k_slam)
k_masked=rawcoldata$masked==masked
table(k_masked)
k_mutated=rawcoldata$mutated==mutated
table(k_mutated)
k_spliceaware=rawcoldata$splice_aware==spliceaware
table(k_spliceaware)

# colnames(rawcoldata)

k = k_slam & k_masked & k_mutated & k_spliceaware
table(k)
cols=unlist(strsplit(indexcols,","))
sample_name=rawcoldata[k,]$sampleName
condition=substr(sample_name,start=1,stop=2)
cols=c(cols,sample_name)
newcounts=rawcounts[,cols]
newcoldata=data.frame(sample_name=sample_name,condition=condition)
write.table(newcounts,file = paste(outdir,countsfn,sep="/"),row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
write.table(newcoldata,file = paste(outdir,coldatafn,sep="/"),row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
# q()


htmlfn=paste(contrast,"html",sep = ".")
htmlfilepath=paste(outdir,htmlfn,sep = "/")
# print(htmlfilepath)

# print(args$cpmcutoff)
rmarkdown::render(
  input  = paste(scriptdir,'DESeq2_DEG.Rmd',sep="/"),
  output_format = 'html_document',
  output_file = htmlfilepath,
  params = list(
    rawcountsmatrix = paste(outdir,countsfn,sep="/"),
    coldata = paste(outdir,coldatafn,sep="/"),
    condition1 = condition1,
    condition2 = condition2,
    indexcols = indexcols,
    cpm_cutoff = cpmcutoff,
    degoutdir = outdir
  )
)
cmd=paste0("mv ",outdir,"/DEG.",contrast,".tsv ",outdir,"/",degfn)
system(cmd)
# q()
# 
# htmlfilepath="/Volumes/Wolin/mESC_slam_analysis/find_mutation_101521/rsem/KO_vs_WT.html"
# rawcountsmatrix="/Volumes/Wolin/mESC_slam_analysis/find_mutation_101521/rsem/rsem.raw_counts_matrix.tsv"
# coldata="/Volumes/Wolin/mESC_slam_analysis/find_mutation_101521/rsem/rsem.raw_counts_matrix.tsv.colData"
# indexcols="gene_id,gene_name"
# condition1="KO"
# condition2="WT"
# cpmcutoff=1
# outdir="/Volumes/Wolin/mESC_slam_analysis/find_mutation_101521/rsem/"
# rmarkdown::render(
#   input  = "/Users/kopardevn/Documents/Projects/rbl4/find_mutations_101521/DESeq2_DEG.Rmd",
#   output_format = 'html_document',
#   output_file = htmlfilepath,
#   params = list(
#     rawcountsmatrix = rawcountsmatrix,
#     coldata = coldata,
#     condition1 = condition1,
#     condition2 = condition2,
#     indexcols = indexcols,
#     cpm_cutoff = cpmcutoff,
#     degoutdir = outdir
#   )
# )
