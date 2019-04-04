setwd("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Analysis_Summary/")

ScriptDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Analysis_Summary/Scripts/"

create.pbs <- function(ScriptDir, tissue){
    
    file <- paste0(ScriptDir, tissue, "_summary.pbs")
    sh.txt <- rbind(
        "#!/bin/bash",
        "#PBS -l  nodes=1:ppn=1,mem=16gb",
        "#PBS -l walltime=96:00:00",
        paste0("#PBS -o ", ScriptDir),
        "#PBS -j oe",
        paste0("#PBS -N ",tissue),
        "module load R")
    
    write.table(file = file, sh.txt, quote=F, row.names=F, col.names=F)

    return(file=file)
}

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

tissues <- read.table(file = paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "data/support_files/all_v8_tissues_both_sexes.txt"), as.is=T)
tissues <- unlist(tissues)


for (tissue in tissues){
    samp.size <- check.sample.size(tissue)
    if (min(unlist(samp.size)) < 40){ next }
    
    cmd <-  paste0("R CMD BATCH  '--args ",tissue,
                   "' /gpfs/data/stranger-lab/askol/GTEx/Coexpression/Code/",
                   "Compare_WGCNA_ARACNE_BicMix.r ",
                   " /gpfs/data/stranger-lab/askol/GTEx/Coexpression/Analysis_Summary/",
                   tissue,".out")
    
    print(paste0("Working on tissue ",tissue))
    
    file <- create.pbs(ScriptDir, tissue)
    write.table(file = file, cmd, quote=F, row.names=F, col.names=F, append=T)

    cmd <- paste0("qsub ", file)
    system(cmd)
}

