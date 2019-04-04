#!/apps/software/gcc-6.2.0/R/3.5.0/bin/Rscript
##
##


check.sample.size <- function(tissue){

    OrigDataDir <- paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "sexDE_v8_final/meri/data/")

    count_file <- paste0(OrigDataDir, "Phenotypes/",tissue,
                         ".v8.gene_counts.txt.gz")
    covs_file <- paste0(OrigDataDir, "/Covariates/covs.basalcovs.",tissue,
        ".txt")
    sex.info <- read.table(file = covs_file, as.is=T, header=F, quote="",
                           comment.char = "", sep="\t")
    colnames(sex.info) = c('SUBJID', 'SEX','SMTSISCH','SMRIN','AGE')
    samps <- read.table(file = count_file, nrow=1, header=F, as.is=T)
    samps <- gsub("\\.","-",samps)
    
    males <- samps[samps %in% sex.info$SUBJID[sex.info$SEX=="M"]]
    females <- samps[samps %in% sex.info$SUBJID[sex.info$SEX=="F"]]

    wcs <- list(males = length(males), females=length(females))

    return(wcs)
}


ResultDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/Aracne/"

tissues <- read.table(file = paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "data/support_files/all_v8_tissues_both_sexes.txt"), as.is=T)
tissues <- unlist(tissues)

job.files = c()

for (tissue in tissues[-1]){    

    sampSize <- check.sample.size(tissue)
    if (any(sampSize < 40)){
        print(paste0("Skipping tissue ",tissue, "due to sample size: ",sampSize))
        next
    }
        
    out.dir <- paste0(ResultDir,tissue,"/")
    if (!file.exists(out.dir)){
        dir.create(out.dir)
    }
    d <-paste0(out.dir,"Male/")
    if (!file.exists(d)){
        dir.create(d, recursive=T)
    }
    d <-paste0(out.dir,"Female/")
    if (!file.exists(d)){
        dir.create(d, recursive=T)
    }                       

    system("module load aracne-ap")
    
    cmd <-  paste0("R CMD BATCH  '--args ", tissue,
                   "' /gpfs/data/stranger-lab/askol/GTEx/Coexpression/Code/Expression_Aracne.r ",
                   ResultDir,
                   tissue,".out&")
    system(cmd)
}


 
