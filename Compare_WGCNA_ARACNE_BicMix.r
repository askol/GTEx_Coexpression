
SummaryDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Analysis_Summary/"
GSEADir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/GSEA/"

## GET NECESSARY ANALYSIS OUTPUT ##
compare.results <- function(tissue){

    SummaryDir <- paste0(SummaryDir, tissue)
    if (!dir.exists(SummaryDir)){
        dir.create(SummaryDir)
    }
    setwd(SummaryDir)

 
    source("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Code/Compare_WGCNA_ARACNE_BicMix_funcs.r")    
    
    ## GENE INFO ##
    gene.annot.file <- paste0("/gpfs/data/stranger-lab/askol/GTEx/",
                              "Coexpression/Data/", tissue,
                              ".gene.annot.txt")
    
    gene.info <- read.table(file = gene.annot.file, as.is=T, header=T)
    
    arac.dir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/Aracne/"
    wgcna.dir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/WGCNA/"
    
    ## DETERMINE NUMBER OF MALE / FEMALE SUBJECTS
    Ns <- check.sample.size(tissue)
    
    ## GET WGCNA DATA ##
    print(paste0("Reading in wgcna data for ",tissue," . . ."))
    wgcna <- get.wgcna.data(tissue, wgcna.dir, gene.info)

    print("Extracting data falling in 97.5 percentile. . .")
    wgcna95 <- get.percentile(wgcna, 0.975)

    rm(wgcna)

    print("Performing regression on wgcna data. . .")
    tmp <- regr.and.errors(wgcna95, Ns)
    wgcna95 <- tmp[[1]]
    names(wgcna95) <- gsub("SRes","SRes.wgcna",names(wgcna95))
    wgcna.R2 <- tmp$R2; wgcna.coefs <- tmp$coefs
    
    ## GET ARACNE DATA ##
    print("Reading in aracne data. . .")
    arac <- get.aracne.data(tissue, arac.dir, gene.info)
    print("Performing regression on aracne data. . .")
    tmp <- regr.and.errors(arac, Ns)
    arac <- tmp[[1]]
    arac.R2 <- tmp$R2; arac.coefs <- tmp$coefs
    print("Performing logistic regression on aracne data. . .")
    tmp <- logreg.and.errors(arac, Ns)
    arac <- tmp[[1]]
    names(arac) <- gsub("SRes","SRes.arac",names(arac))
    arac.R2.miss = tmp$R2
    
    ## GET BicMix DATA ##
    err.thr.a <- calc.max.cutoff(sum(!is.na(arac$MI.1) & !is.na(arac$MI.2)), P=0.25)
    err.thr.w <- calc.max.cutoff(nrow(wgcna95), P=0.25)
    print("Identifying gene pairs that appear as outlier . . .")
    diffs <- find.sex.differences(arac, wgcna95,
                                 err.thr.a = err.thr.a,
                                 err.thr.w = err.thr.w,
                                  err.thr.miss = 2.5)

    file <- paste0(tissue,"_diffs.txt")
    print(paste0("Writing data to file ",file))
    write.table(file = file, diffs, quote=F, row.names=F, col.names=T)

    arac.dist <- calc.connections(arac, Ns)
    names(arac.dist) <- gsub("count","count.arac",names(arac.dist))
    wgcna.dist <- calc.connections(wgcna95, Ns)
    names(wgcna.dist) <- gsub("count", "count.wgcna", names(wgcna.dist))
    dist <- merge(arac.dist, wgcna.dist, by="TF", all=T)

    file <- paste0(tissue,"_base_TF_connections.txt")
    print(paste0("Writing file about distribution of targets to TFs. . ."))
    write.table(file = file, dist, quote=F, row.names=F, col.names=T)
    
    ## GET RANKINGS FOR GSEA ANALYSIS ##
    gsea.scores <- make.gsea.scores(arac, wgcna95)
    ## WRITE OUT ALL RANKS TO GSEA DIRECTOR ##
    rank.names <- names(gsea.scores)[grep("rank", names(gsea.scores))]
    
    RankDir <- paste0(GSEADir, gsub("-","_",tissue), "/Ranks/")    
    GSEAOutDir <- paste0(GSEADir, gsub("-","_",tissue), "/")

    if (!dir.exists(RankDir)){
        dir.create(RankDir, recursive=T)
    }
    if (!dir.exists(GSEAOutDir)){
        dir.create(GSEAOutDir)
    }
    
    gsea.scores <- merge(gsea.scores, gene.info[,c("SYMBOL","ENSEMBL")],
                         by.x="gene", by.y="ENSEMBL", all.x=T, all.y=F)
    for (nm in rank.names){
        
        file <- paste0(RankDir, nm, ".rnk")
        tmp <- gsea.scores[,c("SYMBOL", nm)]       
        rm.ind <- which(is.na(tmp[,2]))
        if (length(rm.ind) > 0){
            tmp <- tmp[-rm.ind,]
        }
        ord <- order(tmp[,2], decreasing=T)
        tmp <- tmp[ord,]
        
        write.table(file = file, tmp, quote=F, col.names=F, sep="\t",
                    row.names=F)
    }
    
    ## RUN GSEA ##
    run.gsea.all(tissue, RankDir, GSEAOutDir)
    run.gsea.custom(tissue, RankDir, GSEAOutDir)
    process.gsea.output(tissue, RankDir, GSEAOutDir)
    
}

## ---- MAIN -----##

source("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Code/Compare_WGCNA_ARACNE_BicMix_funcs.r")
args <- commandArgs(TRUE)
tissue = args[1]
compare.results(tissue)
## run.gsea.all(tissue)
print(date())
## compare.results(tissue)

