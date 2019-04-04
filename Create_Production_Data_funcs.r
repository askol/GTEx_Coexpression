## BASED ON MERIS CODE IN: /gpfs/data/gtex-group/
## sex_biased_regulation_v8/sexDE_v8_final/meri/
## software/calculate.sex.de.genes.no_inv_norm.cmd

system("module load htslib")
system("module load bcftools")
system("module load bedtools")

library(limma)
library(stringr)
library(edgeR)
library(readr)
library(biomaRt)
library(SmartSVA)
library(annotables)
library(dplyr)

CodeDir <- paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                    "sexDE_v8_final/meri/software/")
DataDir <- paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                      "sexDE_v8_final/meri/data/")
OutDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Data/"

Create_Production_Data <- function(tissue){  
    
    make.expr.files(CodeDir, tissue)

}

make.expr.files <- function(CodeDir, tissue){

    covs_file=paste0(DataDir, "/Covariates/covs.basalcovs.",tissue,
        ".txt")

    count_file <- paste0(DataDir, "Phenotypes/",tissue,
                         ".v8.gene_counts.txt.gz")
    
    ## Load Covariates
    covs=read.table(file=covs_file, header=F, row.names=1,
        check.names=F)
    colnames(covs) = c('SEX','SMTSISCH','SMRIN','AGE')
    covs$SEX = as.numeric(covs$SEX == "F") + 1
    covs$SMTSISCH[is.na(covs$SMTSISCH)] = median(covs$SMTSISCH,na.rm=T)
    covs$SMRIN[is.na(covs$SMRIN)] = median(covs$SMRIN , na.rm=T)
    form = "~ SEX+SMTSISCH+SMRIN+AGE"
    if ( sum(as.numeric(!is.na(covs$SMTSISCH))) < 1 ) {
        covs$SMTSISCH <- NULL;
        form = "~ SEX+SMRIN+AGE";
        print("Ischemic time not available. Drop variable")
    }

    ## Load gene filtered counts
    ##    exp_genes_unorm = read.table(pipe('$cmd'),header=T,check.names=F,row.names=1)
    exp_genes_unorm = read_tsv(count_file, col_names=FALSE, skip=1)
    cnames <- read_tsv(count_file, col_names=FALSE, n_max=1)
    cnames <- str_extract(cnames,'\\w+-\\w+')
    names(exp_genes_unorm) <- c("GENE", cnames) ## gsub("-",".",cnames))
    exp_genes_unorm$GENE <- gsub("\\.[0-9].*$", "", exp_genes_unorm$GENE)

    ## GET GENE INFORMATION
    ## gene.info <- get.gene.annotation(tissue, genes, rerun=F)
    gene.info <- grch38
    gene.info <- filter(gene.info, ensgene %in% exp_genes_unorm$GENE) %>%
        mutate(ID = paste(ensgene,symbol,sep=".")) %>%
            filter(duplicated(ID)==FALSE) %>% dplyr::select(-ID,-entrez,-strand)
    
    ## REMOVE GENES NOT IN GENE ANNOTATION ##
    ind <- which(exp_genes_unorm$GENE %in% gene.info$ensgene)
    exp_genes_unorm <- exp_genes_unorm[ind,]
    
    ## Subset samples with exp+metadata
    rownames(covs) = str_extract(rownames(covs),'\\w+-\\w+')
    common = rownames(covs)[rownames(covs) %in%
        colnames(exp_genes_unorm)]
    covs = covs[common,]
    exp_genes_unorm = exp_genes_unorm[,c("GENE", common)]

    ## Return edgeR normalized/rescaled CPM (counts per million)
    design <- model.matrix(eval(parse(text=form)), data=covs)
    colnames(design) <- c('Intercept',
                          gsub('covs','',colnames(design)[-1]))
    rownames(design) <- rownames(covs)
    exp_genes <- DGEList(counts=exp_genes_unorm[, -grep("GENE",
                             names(exp_genes_unorm))],
                         genes=exp_genes_unorm$GENE,
                         group=covs$SEX)

    ## ADD COVARIATES INTO OBJECT
    mtch <- match(rownames(exp_genes$samples), rownames(covs))
    exp_genes$samples <- cbind(exp_genes$samples, covs[mtch,])
    
    ## ADD GENE INFO ##
    mtch <- match(exp_genes$gene$genes, gene.info$ensgene)
    exp_genes$genes <- cbind( exp_genes$gene$genes, gene.info[mtch, ])   
    
    ## REMOVE GENES IS > 20% OF SAMPLES HAVE 
    propGTE6 <-  apply(exp_genes$counts, 1, function(x) mean(x >= 6))
    keep <- which(propGTE6 >= .2)
    exp_genes <- exp_genes[keep,]
    exp_genes <- calcNormFactors(exp_genes, method="TMM")
    
    ## ADD SURAGATE VARIABLES ##
    svs <- generate.svs(exp_genes, covs=names(covs))
    exp_genes$samples <- cbind(exp_genes$samples, svs)
    
    ## Merge covs
    ## covs <- cbind(covs, svs)
	
    ## Return edgeR normalized/rescaled CPM (counts per million)
    covs <- as.data.frame(exp_genes$samples)
    covs <- dplyr::select(covs, -SEX, -group, -lib.size, -norm.factors)
    lcpm <- cpm(exp_genes, log=TRUE)
    lcpmNorm <- removeBatchEffect(lcpm, covariates=covs)        
    rownames(lcpmNorm) <- exp_genes$genes$gene
    genes <- exp_genes$genes$genes

    ## REMOVE X CHROMOSOME GENES ##
    ##xy.genes <- gene.info$ENSEMBL[gene.info$chr %in% c("X","Y")]    
    ##auto.ind <- which(!rownames(lcpmNorm) %in% xy.genes)
    ##lcpmNormAuto <- lcpmNorm[auto.ind,]

    saveRDS(exp_genes, file = paste0(OutDir , tissue, ".RDS"))
            
    ## MAKE ARACNE FILES ##
    make.aracne.files(lcpmNorm, exp_genes, tissue, OutDir)
    make.aracne.files(lcpm, exp_genes, tissue, OutDir, suff="notNmlz")

    ## MAKE WGCNA FILES ##
    make.wgcna.files(lcpmNorm, exp_genes, tissue, OutDir)
    make.wgcna.files(lcpm, exp_genes, tissue, OutDir, suff="notNmlz")
}

get.gene.annotation <- function(tissue, genes, rerun=F){
        
    gene.annot.file <- paste0("/gpfs/data/stranger-lab/askol/GTEx/",
                              "Coexpression/Data/", tissue,
                              ".gene.annot.txt")
    
    if (rerun==F & file.exists(gene.annot.file)){
        
        print("Reading from file. Use rerun=T to regenerate file.")
        gene.info <- read.table(file = gene.annot.file, as.is=T, header=T)
        
    }else{
        
        ## GET ALL POSSIBLE GENE IDS IN GTEX DATA ##
        geneid <- gsub("\\.[0-9].*$","", genes)
        
        ## GET GENE ANNOTATION ##
        ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
        mart <- useMart("ensembl", host = "www.ensembl.org")
        mart <- useDataset(dataset = 'hsapiens_gene_ensembl',
                           mart = mart)        
        attribs = c("ensembl_gene_id", "hgnc_symbol",
            "chromosome_name", "start_position","end_position", "gene_biotype")
        gene.info <- getBM(attributes = attribs,
                           filters = 'ensembl_gene_id',
                           values = geneid, mart = mart)
        names(gene.info) <- c("ENSEMBL", "SYMBOL", "chr", "start", "end",
                              "gene_biotype")
        
        ## REMOVE GENES WITH NO SYMBOL ##
        rm.ind <- which(is.na(gene.info$SYMBOL) | gene.info$SYMBOL == "")
        if (length(rm.ind)>0){
            gene.info <- gene.info[-rm.ind, ]
        }
        
        ## REMOVE DUPLICATES
        rm.ind <- which(duplicated(
            paste(gene.info$ENSEMBL, gene.info$chr, gene.info$SYMBOL,
                  sep=".")))
        if (length(rm.ind) > 0 ){
            gene.info <- gene.info[-rm.ind,]
        }
        
        write.table(file = gene.annot.file, gene.info, quote=F,
                    row.names=F, col.names=T)
    }    
    return(gene.info)
}

make.aracne.files <- function(lcpm, data, tissue, OutDir, suff=""){

    
    file <-  paste0(OutDir,"Aracne/",tissue,"_lcpm_auto_male.txt")
    if (suff != ""){
        file <- gsub("_lcpm",paste0("_",suff,"_lcpm"), file)
    }
    male.ind <- which(data$samples$group==1)
    lcpm.male <- lcpm[,male.ind]
    lcpm.male <- cbind(rownames(lcpm.male), lcpm.male)
    colnames(lcpm.male)[1] <- "gene"
    write.table(file=file, lcpm.male, quote=F, row.names=F, col.names=T, sep="\t")

    file <-  paste0(OutDir,"Aracne/",tissue,"_lcpm_auto_female.txt")
    if (suff != ""){
        file <- gsub("_lcpm",paste0("_",suff,"_lcpm"), file)
    }
    lcpm.female <- lcpm[,-male.ind]
    lcpm.female <- cbind(rownames(lcpm.female), lcpm.female)
    colnames(lcpm.female)[1] <- "gene"
    write.table(file=file, lcpm.female, quote=F, row.names=F, col.names=T, sep="\t")


}

make.wgcna.files <- function(lcpm, data, tissue, OurDir, suff=""){

    file <-  paste0(OutDir,"WGCNA/",tissue,"_lcpm_auto_male.txt")
    if (suff != ""){
        file <- gsub("_lcpm",paste0("_",suff,"_lcpm"), file)
    }
    male.ind <- which(data$samples$group==1)
    lcpm.male <- lcpm[,male.ind]
    write.table(file=file, lcpm.male, quote=F, row.names=T, col.names=T, sep="\t")

    file <-  paste0(OutDir,"WGCNA/",tissue,"_lcpm_auto_female.txt")
    if (suff != ""){
        file <- gsub("_lcpm",paste0("_",suff,"_lcpm"), file)
    }
    lcpm.female <- lcpm[,-male.ind]
    write.table(file=file, lcpm.female, quote=F, row.names=T, col.names=T, sep="\t")

}


## FROM MERI'S CODE IN: ##
## /gpfs/data/gtex-group/sex_biased_regulation_v8/sexDE_v8_final/meri/
## software/generate.svs.R    
generate.svs <- function(data, covs=c()) {
    
    ## Return edgeR normalized/rescaled CPM (counts per million)
    form <- paste("~", paste(covs, collapse="+"))
    design <- model.matrix(eval(parse(text=form)), data=data$sample)
    count.log2.cpm <- voom(counts=data, design=design)
    
    ## Inverse-rank-normalize counts
    norm_counts <- inverse_quantile_normalization(count.log2.cpm$E)
    
    ## Calculate smartSVs blocking for covs
    ## Determine the number of SVs
    cmd <- paste0('Y.r <- t(resid(lm(t(norm_counts) ',form,
                  ' , data=data$samples)))')
    eval(parse(text=cmd))
    
    ## Add one extra dimension to compensate potential loss of 1 degree of freedom
    ## in confounded scenarios (very important)
    n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
    mod <- model.matrix(eval(parse(text=form)),data=data$samples)
    
    ## Modify the default parameters: iteration numbers (B) and learning rate (alpha)
    sv.obj <- smartsva.cpp(norm_counts, mod, mod0=NULL, n.sv=n.sv, B=1000, alpha=1)
    
    allSv <- sv.obj$sv
    colnames(allSv) <- paste0("SV", 1:n.sv)
    rownames(allSv) <- colnames(data)

    return(allSv)
}

inverse_quantile_normalization <- function(gct) {
        gct = t(apply(gct, 1, rank, ties.method = "average"));
        gct = qnorm(gct / (ncol(gct)+1));
        return(gct)
}

