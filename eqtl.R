ibrary(dplyr)
library(tidyr)
library(data.table)
library(liftOver)
library(RACER)
###perform liftover on gwas snps and prepare input for coloc analysis 
morris_snps_chr6 <- fread("/scratch/aa9gj/wrk_backup/eya4_project/morris_data/chr6_bms_snps.txt")
chr6_lookup <- fread("/scratch/aa9gj/wrk_backup/chr6_lookup_table_use", header = FALSE)
colnames(chr6_lookup) <- c("id", "rsid")
ch <- import.chain("/scratch/aa9gj/wrk_backup/hg19ToHg38.over.chain")
colnames(morris_snps_chr6)[3] <- "seqnames"
cur6 <- makeGRangesFromDataFrame(morris_snps_chr6, keep.extra.columns = TRUE, start.field = "V4", end.field = "V4")
seqlevelsStyle(cur6) = "UCSC"  # necessary
cur38_6 = liftOver(cur6, ch)
class(cur38_6)
cur38_6 = unlist(cur38_6)
genome(cur38_6) = "hg38"
cur38_6 = new("gwaswloc", cur38_6)
cur38_6 <- as.data.table(cur38_6)
cur38_6$maf <- ifelse(cur38_6$V7 > 0.5, 1 - cur38_6$V7, cur38_6$V7)
cur38_6$seqnames <- gsub("chr", "", cur38_6$seqnames)
cur38_6$V13 <- as.numeric(cur38_6$V13)
cur38_6 <-cur38_6[,-3]
cur38_6 <-cur38_6[,-4]
cur38_6[order(cur38_6$V13)]

##reading in list of files from command line to automate
setwd("/scratch/aa9gj/wrk_backup/eqtl_colocalization/Gtex_v8/")
my_files <- read.delim("files", header=FALSE)#

##functions written to automate colocalization analysis
##prep_eqtl takes eqtl files from specific genes (eya4 in my case) amd prep them for coloc analysis  
prep_eqtl <- function(x,y) {
  colnames(x) <- c("gene", "id", "tss_dist", "ma_sample", "ma_count", "maf", "pval_nominal", "slope", "slope_se")
  x <- inner_join(x, y, by = "id")
  x <- filter(x, tss_dist >= -200000 & tss_dist <= 200000)
  x <- separate(x, id, into=c("chr", "position", "ref", "alt", "genome", sep = "_"))
  x$chr <- gsub("chr", "", x$chr)
  x<-x[order(x$pval_nominal),]
  return(x)
}
##coloc performs coloc analysis x= lifover morris snps, y=prepped eqtl file
coloc_tissue_bmd <- function(x,y,n) {
  coloc_res <- coloc.abf(dataset1=list(snp=x$V2,type="quant", MAF=x$maf, N=426824 , pvalues = x$V13), 
                         dataset2=list(snp=y$rsid,type="quant", MAF=y$maf, N=n, pvalues = y$pval_nominal))
  return(coloc_res)
}
##use this method to automate functions and reading files
##for reading files 
my_data <- list()
for (i in seq_along(my_files$V1)) {
  my_data[[i]] <- read.delim(file = my_files$V1[i], header=FALSE)
}

split_df<-function(list){
  for (i in 1:length(list)){
    assign(paste0(my_files$V1[i]), list[[i]], envir = .GlobalEnv)
  }
}
split_df(my_data)

##autamation of prep_eqtl
my_eqtl <- list()
for (i in seq_along(my_files$V1)) {
  my_eqtl[[i]] <- prep_eqtl(my_data[[i]], chr6_lookup)
}

split_df_eqtl<-function(list){
  for (i in 1:length(list)){
    assign(paste0(my_files$V1[i], "_prepped"), list[[i]], envir = .GlobalEnv)
  }
}

##automation of coloc
split_df_coloc<-function(list){
  for (i in 1:length(list)){
    assign(paste0(my_files$V1[i], ".coloc"), list[[i]], envir = .GlobalEnv)
  }
}

split_df_coloc<-function(list){
  for (i in 1:length(list)){
    assign(paste0(my_files$V1[i], ".coloc"), list[[i]], envir = .GlobalEnv)
  }
}

##mirrorplot code completed
mark3_bmd_gwas_f = RACER::formatRACER(assoc_data = cur38_6, chr_col = 1, pos_col = 2, p_col = 14)
mark3_eqtl_f = RACER::formatRACER(assoc_data = ENSG00000112319.17_Muscle_Skeletal.eqtl_prepped, chr_col = 2, pos_col = 3, p_col = 12)
mark3_bmd_gwas_f_ld = RACER::ldRACER(assoc_data = mark3_bmd_gwas_f, rs_col = 5, pops = "EUR", lead_snp = "rs9375951")
mark3_eqtl_f_ld = RACER::ldRACER(assoc_data = mark3_eqtl_f, rs_col = 15, pops = "EUR", lead_snp = "rs9375951")
mirrorPlotRACER(assoc_data1 = mark3_bmd_gwas_f_ld, assoc_data2 = mark3_eqtl_f_ld, chr = 6, plotby = "gene", gene_plot = "EYA4", label_lead = TRUE,build = "hg38")
mirrorPlotRACER(assoc_data1 = mark3_bmd_gwas_f_ld, assoc_data2 = mark3_eqtl_f_ld, chr = 6, plotby = "coord", start_plot = 133000000, end_plot = 133750000, label_lead = TRUE, build = "hg38")
