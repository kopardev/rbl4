#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

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


suppressPackageStartupMessages(library("rmarkdown"))

contrast=paste(args$condition1,"vs",args$condition2,sep="_")
htmlfn=paste(contrast,"html",sep = ".")
htmlfilepath=paste(args$outdir,htmlfn,sep = "/")
print(htmlfilepath)
print(args$cpmcutoff)
rmarkdown::render(
  input  = 'DESeq2_DEG.Rmd',
  output_format = 'html_document',
  output_file = htmlfilepath,
  params = list(
    rawcountsmatrix = args$rawcountsmatrix,
    coldata = args$coldata,
    condition1 = args$condition1,
    condition2 = args$condition2,
    indexcols = args$indexcols,
    cpm_cutoff = args$cpmcutoff,
    degoutdir = args$outdir
  )
)
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
