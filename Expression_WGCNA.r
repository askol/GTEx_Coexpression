### DETEMRINE ARGUEMENT FOR ANALYSIS ###
###
source("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Code/Expression_WGCNA_funcs.r")


Expression_WGCNA <- function(tissue){

    CodeDir <- paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                      "sexDE_v8_final/meri/software/")
    DataDir <- paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                  "sexDE_v8_final/meri/data/")
    OutDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Data/WGCNA/"
    ResultDir <- paste0("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/WGCNA/")                        
    
    if (!exists(ResultDir)){
        dir.create(ResultDir, tissue)
    }
    ResultDir <- paste0(ResultDir,tissue,"/")
    
    ## DATA WAS CREATED BY CREATE_PRODUCTION_DATA.R
    data <- get.data(tissue, DataDir = OutDir)

    multiExpr <- data$multiExpr

    ## DETERMINE 'BEST' SOFT THRESHOLD TO CREATE ADJACENCY
    print("Identifying soft thresholds for males and females")
    softThr <- find.soft.threshold(multiExpr, ResultDir)

    ## CALCULATE SCALED TOMS ##
    TOMs <- calculate.TOM(multiExpr, softThr)

    saveRDS(TOMs ,
            file="/gpfs/data/stranger-lab/askol/TCGA/Expression_Analysis/WGCNA/TOMs.RDS")
    
    ## OUTPUT TOM VALUES AS A TABLE WITH ALL TF-GENE PAIRS ##
    TFs <- read.table(file =
                      paste0("/gpfs/data/stranger-lab/askol/TCGA/",
                             "Expression_Analysis/ARACNE/TFs_Ensembl_v_1.01.txt"),
                      as.is=T, header=T)
    TFs <- as.character(unlist(TFs))
    
    output.TOMs(TOMs, TFs, ResultDir)
    
    
    print(paste0("Completed Expression_WGCNA for tissue: ",tissue))
}

## MAIN ##

args <- commandArgs(TRUE)
project = args[1]
Expression_WGCNA(project)
print(date())

finish.file = paste0("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/WGCNA/",
    tissue,".finished")
system(paste0("touch ",finish.file))
