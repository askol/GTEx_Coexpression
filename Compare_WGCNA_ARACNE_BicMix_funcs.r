library(edgeR)
library(limma)
library(dplyr)
library(survival)
library(survminer)
library(WGCNA)
library(pscl)
library(reshape2)

GSEADir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/GSEA/"
source("~/Code/qq_plot.r")

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


get.aracne.data <- function(tissue, arac.dir, gene.info, prob.filt.thresh=c(.01,.01)){

    ## NETWORKS FROM ARACNE ##
    files <- paste0(paste0(arac.dir,tissue,"/",c("Male/","Female/")),
                    "network.txt")
    
    network <- list()
    for (i in 1:2){

        file <- files[i]    
        network[[i]] <- read.table(file = file, header=T, as.is=T)      
        network[[i]]$name <- paste(network[[i]]$Regulator, network[[i]]$Target, sep=".")
        ind <- which(names(network[[i]]) %in% c("MI","pvalue"))
        names(network[[i]])[ind] <- paste(names(network[[i]])[ind], i, sep=".")
    }

    networks <- merge(network[[1]], network[[2]][,-c(1:2)], by="name", all=T)
    networks$Regulator <- gsub("\\..+","",networks$name)
    networks$Target <- gsub(".+\\.","",networks$name)

    rm(network)

    ## APPLY PROBABILITY THRESHOLDS : SET TO MISSING GENE PAIRS WITH VALUES
    ## LESS THAN THRESHOLDS

    for (i in 1:2){
        p.ind <- grep(paste0("pvalue.",i), names(networks))
        mi.ind <- grep(paste0("MI.",i), names(networks))
        ind <- which(networks[,p.ind] > prob.filt.thresh[i])
        networks[ind,mi.ind] <- NA
    }


    ## REMOVE GENE PAIRS THAT BOTH HAVE NAS FOR MI ##
    ind <- which(is.na(networks$MI.1) & is.na(networks$MI.2))
    networks <- networks[-ind,]    

    ## ADD POSITION FOR TFS
    networks <- merge(networks, gene.info[,c("ENSEMBL","chr","start")],
                      by.x="Regulator", by.y="ENSEMBL", all.x=T, all.y=F)
    names(networks)[ncol(networks)-c(1,0)] <- c("Reg.chr","Reg.pos")

    ## ADD POSITIONS FOR TARGETS
    networks <- merge(networks, gene.info[,c("ENSEMBL","chr","start")],
                      by.x="Target", by.y="ENSEMBL", all.x=T, all.y=F)
    names(networks)[ncol(networks)-c(1,0)] <- c("Targ.chr","Targ.pos")        
    
    return(networks)
}

get.wgcna.data <- function(tissue, wgcna.dir, gene.info){

    files <- paste0(paste0(wgcna.dir,tissue,"/",c("Male/","Female/")),"TOM_table.txt")
    network <- list()
    for (i in 1:2){

        file <- files[i]    
        network[[i]] <- read.table(file = file, header=T, as.is=T)
        rm.ind <- which(network[[i]]$TF == network[[i]]$Target)
        network[[i]] <- network[[i]][-rm.ind,]
        
        network[[i]]$name <- paste(network[[i]]$TF, network[[i]]$Target, sep=".")
        ind <- which(names(network[[i]]) == "TOM")
        names(network[[i]])[ind] <- paste(names(network[[i]])[ind], i, sep=".")
    }
    
    networks <- merge(network[[1]], network[[2]][,-c(1:2)], by="name", all=T)
    networks$TF <- gsub("\\..+","",networks$name)
    networks$Target <- gsub(".+\\.","",networks$name)

    rm(network)

    ## ADD POSITION FOR TFS
    networks <- merge(networks, gene.info[,c("ENSEMBL","chr","start")],
                      by.x="TF", by.y="ENSEMBL", all.x=T, all.y=F)
    names(networks)[ncol(networks)-c(1,0)] <- c("Reg.chr","Reg.pos")

    ## ADD POSITIONS FOR TARGETS
    networks <- merge(networks, gene.info[,c("ENSEMBL","chr","start")],
                      by.x="Target", by.y="ENSEMBL", all.x=T, all.y=F)
    names(networks)[ncol(networks)-c(1,0)] <- c("Targ.chr","Targ.pos")        

    return(networks)
}

get.percentile <- function(data, quant.cut = .95){

    ## IDENTIFY GENE PAIRS WITH ASSOCIATION VALUE IN THE TOP
    ## QUANT.CUT PERCENTILE AMONG PAIRS OF GENES ON DIFFERENT
    ## CHROMOSOMES
    ##
    ## GOAL IS TO REDUCE GENE PAIRS TO THOSE MOST LIKELY TO BE
    ## COEXPRESSED FOR REALS
    
    value.ind <- grep("TOM|MI", names(data))
    chr.ind <- data$Reg.chr != data$Targ.chr
    
    ## once for males, once for females
    quant <- c()
    quant.ind <- list()
    for (i in 1:2){
        ## which gene pairs are on different chromosomes

        quant[i] <- quantile(data[chr.ind,value.ind[i]], probs=quant.cut, na.rm=T)
        quant.ind[[i]] <- which(data[,value.ind[i]] >= quant[i])
    }

    keep <- unique(c(quant.ind[[1]], quant.ind[[2]]))

    return(data[keep,])
}


regr.and.errors <- function(data, Ns){

    ## Ns are the number of males and females, respecitively
    sexes <- c("Male","Female")
    max.n.ind <- order(unlist(Ns), decreasing=T)[1]

    ## THE LARGER SAMPLE SIZE VALUE WITH SERVE AS THE X VARIABLE SINCE IN PRINCIPLE
    ## IT IS MEASURED WITHOUT ERROR
    x.ind <- intersect(grep(max.n.ind, names(data)), grep("MI|TOM",names(data)))
    y.ind <- intersect(grep(3-max.n.ind, names(data)), grep("MI|TOM",names(data)))
    
    reg <- lm(data[,y.ind] ~ data[,x.ind], na.action=na.exclude)
    SRes <- rstudent(reg)
    data$SRes <- SRes
    R2 <- summary(reg)$r.squared
    coefs <- summary(reg)$coefficients

    return(list(data, R2 = R2, coefs = coefs))
}


logreg.and.errors <- function(data, Ns){
    
    sexes <- c("Male","Female")
    max.n.ind <- order(unlist(Ns), decreasing=T)[1]

    x.ind <- intersect(grep(max.n.ind, names(data)), grep("MI|TOM",names(data)))
    calc.miss.on.ind <- intersect(grep(3-max.n.ind, names(data)),
                                  grep("MI|TOM",names(data)))
    
    ## CREATE MISSINGNESS VARIABLE FOR THE SEX WITH SMALLER SAMPLE SIZE
    data$miss <- 1*is.na(data[,calc.miss.on.ind])
    reg <- glm(data$miss ~ data[,x.ind], na.action=na.exclude,
               family="binomial")
    SRes <- rstudent(reg)
    E <- predict(reg, type="response")
    R <- E - data$miss
    data$Exp.miss <- E
    data$SRes.miss <- SRes
    
    R2 <- pR2(reg)

    return(list(data, R2=R2))
}
                     
    
find.sex.differences <- function(arac, wgcna,
                                 err.thr.a = 4,
                                 err.thr.w = 4,
                                 err.thr.miss = 2.5){

    wgcna.ind <- which(abs(wgcna$SRes.wgcna) >= err.thr.w)
    arac.ind <- which(abs(arac$SRes.arac) >= err.thr.a )
    arac.miss.ind <- which(abs(arac$SRes.arac.miss) >= err.thr.miss)

    arac.names <- arac$name[arac.ind]
    arac.miss.names <- arac$name[arac.miss.ind]
    wgcna.names <- wgcna$name[wgcna.ind]

    all.names <- unique(c(arac.names, arac.miss.names, wgcna.names))
    
    ind.arac <- which(arac$name %in% all.names)
    ind.wgcna <- which(wgcna$name %in% all.names)

    diffs <- merge(arac[ind.arac, ],
                   wgcna[ind.wgcna,-which(names(wgcna) %in% c("TF","Target"))],
                   by="name", all.x=T, all.y=T)
  
    diffs$arac.chosen <- diffs$arac.miss.chosen <- diffs$wgcna.chosen <- 0

    diffs$arac.chosen[diffs$name %in% arac.names] <- 1
    diffs$arac.miss.chosen[diffs$name %in% arac.miss.names] <- 1
    diffs$wgcna.chosen[diffs$name %in% wgcna.names] <- 1

    ## CLEAN COLUMN NAMES ##
    diffs$Regulator <- gsub("\\..+","",diffs$name)
    diffs$Target <- gsub(".+\\.", "", diffs$name)

    for (gene in c("Reg","Targ")){
        for (Pos in c("chr","pos")){

            col.pre <- paste0(gene,".",Pos)
            col.ind <- grep(col.pre , names(diffs))
            vals <- apply(diffs[,col.ind],1,function(x)unique(x[!is.na(x)]))
            vals <- as.numeric(vals)
            diffs[,col.ind[1]] <- vals
            diffs <- diffs[,-col.ind[-1]]
            names(diffs)[col.ind[1]] <- col.pre
        }
    }
    return(diffs)
    
}

make.gsea.scores <- function(arac, wgcna){

    ## PROCESS WGCNA
    ## CONSOLIDATE POSITIVE SIGNAL ##
    ## ind.pos <- which(wgcna$SRes.wgcna > 0)
    ## ind.neg <- which(wgcna$SRes.wgcna < 0)

    wgcna.stats <- summarize.gsea(wgcna, col="SRes.wgcna")
    names(wgcna.stats)[-1] <- paste0("w.",names(wgcna.stats)[-1]) 
    
    ## PROCESS ARACNE ##
    Arac.stats <- summarize.gsea(arac, col="SRes.arac")
    names(Arac.stats)[-1] <- paste0("a.", names(Arac.stats)[-1])
    
    Arac.miss.stats <- summarize.gsea(arac, col="SRes.arac.miss")
    names(Arac.miss.stats)[-1] <- paste0("am.",names(Arac.miss.stats)[-1])
    
    ## COMBINE RANKS BASED ON MIN RANK ACROSS METHODS ##
    stats <- merge(wgcna.stats, Arac.stats, by="gene", all=T)
    stats <- merge(stats, Arac.miss.stats, by="gene", all=T)

    ## COME UP WITH BEST RANK FROM ALL METHODS ##
    stats$min.rank.best <- best.ranks(stats, "min.rank")
    stats$max.rank.best <- best.ranks(stats, "max.rank")
    stats$any.rank.best <- best.ranks(stats, "any.rank")
    stats$mean.rank.best <- best.ranks(stats, "mean.rank")     
    
    return(stats)
    
}

summarize.gsea <- function(x, col, adj.N=""){

    ## METRICS FOR TFS
    col.ind <- which(names(x) == col)

    ## arac data has TF named Regulator ##
    names(x) <- gsub("Regulator","TF", names(x))    
    
    stats.tf <- calc.stats(x[, col.ind], x$TF)
    stats.targ <- calc.stats(x[, col.ind], x$Target)

    rm.ind <- which(stats.targ$gene %in% stats.tf$gene)
    if (length(rm.ind)>0){
        stats.targ <- stats.targ[-rm.ind,]
    }
    
    stats.tf$type <- "TF"
    stats.targ$type <- "Target"
    
    stats <- rbind(stats.tf, stats.targ)
    
    ## CREATE A SIGNED RANK BASED ON P.ADJ ##
    ## Ranks for p.max ##
    stats$max.rank <- sign.rank(stats$max.p.adj)
    stats$min.rank <- sign.rank(stats$min.p.adj)
    stats$any.rank <- sign.rank(
        (1 - 2*(stats$min.p.adj<stats$max.p.adj)) *
        apply(stats[,c("min.p.adj","max.p.adj")], 1, min))
    stats$mean.rank <- sign.rank(sign(stats$mean) * stats$mean.p.adj)

       return(stats)                        
}

calc.stats <- function(exp, gene, adj.N=""){

    rm.na <- which(is.na(exp))
    if (length(rm.na) > 0){
        exp <- exp[-rm.na]
        gene <- gene[-rm.na]
    }

    tmp.mean <- tapply(exp, gene, mean, na.rm=T)
    tmp.sd <- tapply(exp, gene, sd, na.rm=T)
    tmp.N <- tapply(exp, gene, function(x) sum(!is.na(x)))
    tmp.N.pos <- tapply(exp, gene, function(x) sum(x>0, na.rm=T))
    tmp.prop.pos <- tmp.N.pos / tmp.N
    tmp.prop.pos.p <- 1 - pbinom(tmp.N.pos-1, tmp.N , .5)
    tmp.mean.p <- tapply(exp, gene, function(x)
                         if (sum(!is.na(x)) < 2)
                         { NA } else {
                             if (var(x,na.rm=T) > 0){
                             t.test(x)$p.value } else {NA}})
    tmp.min <- tapply(exp, gene, min, na.rm=T)
    tmp.max <- tapply(exp, gene, max, na.rm=T)
    tmp.min.p <- 1-(1 - pnorm(tmp.min))^tmp.N
    tmp.max.p <- 1 - (pnorm(tmp.max))^tmp.N

    tmp.min.p.adj = NA
    tmp.max.p.adj = NA
    
    ## SHOULD BE TRUE WHEN USING BOTH
    if (length(adj.N) == 1){
        
        adj.N <- apply(cbind(tmp.sd^2*tmp.N,  tmp.N), 1,  min)
        adj.N <- apply(cbind(adj.N, 1), 1, max)
    }
    
    ## ADJUST FOR VARIANCE DEFLATION ##
    ## adjustment based on that E(sigma.hat^2) = N
    tmp.min.p.adj <-  1-(1 - pnorm(tmp.min))^adj.N
    tmp.max.p.adj <- 1 - (pnorm(tmp.max))^adj.N
    
    tmp.mean.p.adj <- pt(qt(tmp.mean.p/2 , tmp.N), adj.N)*2           

    tmp.prop.pos.p.adj <- 1 - pbinom(ceiling(tmp.prop.pos*adj.N)-1,
                                     ceiling(adj.N), .5)
    out <- data.frame(gene = names(tmp.mean), mean=tmp.mean,
                      sd = tmp.sd, N = tmp.N, Prop.pos = tmp.prop.pos,
                      mean.p = tmp.mean.p, Prop.pos.p = tmp.prop.pos.p,
                      mean.p.adj = tmp.mean.p.adj, min = tmp.min,
                      max = tmp.max, min.p = tmp.min.p,
                      max.p = tmp.max.p, min.p.adj = tmp.min.p.adj,
                      max.p.adj = tmp.max.p.adj, prop.pos.p.adj = tmp.prop.pos.p.adj,
                      adj.N = adj.N)
    return(out)
    
}

make.ranks <- function(x){
   
    x$max.rank <- sign.rank(x$max.p.adj)
    x$min.rank <- sign.rank(x$min.p.adj)
    x$best.min.max.rank <- sign.rank(apply(x[,c("min.p.adj","max.p.adj")],
                                      1, min))
    x$mean.rank <- sign.rank(sign(x$mean)*x$mean.p.adj)

    return(x)
}

sign.rank <- function(x){

    ## expecting signed p-values ## 
    pos.ind <- which(x>0)
    ## multiply by -1 so the smallest p-values will have the largest rank
    ## (rank ranks smallest first)
    rank.pos <- rank(x[pos.ind])
    neg.ind <- which(x<0)
    rank.neg <- rank(-1*x[neg.ind])
    z.ind <- which(x == 0)

    ranks <- rep(NA,length(x))
    if (length(pos.ind)>0){
        ranks[pos.ind] <- rank.pos
    }
    if (length(neg.ind)>0){
        ranks[neg.ind] <- -1*rank.neg
    }

    if (length(z.ind) == 1){
        ranks[z.ind] <- 0
    }
    if (length(z.ind) > 1){
        ranks[z.ind] <- seq(-.9,.9, length.out=length(z.ind))
    }

    ## ADD A SMALL AMOUNT OF RANDOMNESS TO BREAK TIES ##
    ranks[!is.na(ranks)] <- ranks[!is.na(ranks)] + .001*runif(sum(!is.na(ranks)),0,1)
    
    return(ranks)
}


    
best.ranks <- function(x, col.name){

    col.ind <- grep(col.name, names(x))
    tmp <- x[,col.ind]
    
    norm <- c()
    tmp.c <- rep(NA, nrow(tmp))
    for (i in 1:3){

        pos.ind <- which(tmp[,i] > 0)
        neg.ind <- which(tmp[,i] < 0)

        if (length(pos.ind)>0){
            tmp.c[pos.ind] <- (rank(tmp[pos.ind,i])+.5)/(length(pos.ind)+.5)
        }
        if (length(neg.ind)>0){
            tmp.c[neg.ind] <- -1*(rank(tmp[neg.ind,i])+.5)/(length(neg.ind)+.5)
        }

        norm <- cbind(norm, tmp.c)
    }
    
    best.rank <- apply(norm,  1, function(x) {o <- order(abs(x), decreasing=T)[1];
                                              return(x[o])})
        
    return(best.rank)
}
        

                           
calc.max.cutoff <- function(N, P=0.25){

    ## N= number of tests ##
    ## P = probability of observing a max value great than inv(P)
    ## E.g. P = 0.5 would give the expected maximum value
    ## ASSUMING THAT SAMPLE IS INDEPENDENT STANDARD NORMAL DISTRIBUTED 

    solv.func <- function(X,N,P){ abs(log(pnorm(X))-log(1-P)/N)}
    tmp <- optimize(solv.func, c(0,10), N=N, P=P)
    return(tmp$minimum)
}

calc.connections <- function(data, Ns){

    TF.ind <- grep("TF|Regulator", names(data))
    val.ind <- grep("MI|TOM", names(data))
    miss.ind <- grep("miss", names(data))

    tbl.miss <- c()
    if (length(miss.ind) > 0){
        
        max.n.ind <- order(unlist(Ns), decreasing=T)[1]
        x.ind <- intersect(grep(max.n.ind, names(data)), grep("MI|TOM",names(data)))
        TF.miss <- data[,TF.ind]
        TF.miss <- TF.miss[!is.na(data[,x.ind])]
        tbl.miss <- as.data.frame(table(TF.miss))
        names(tbl.miss) <- c("TF","miss01.count")
    }

    TF.both <- data[,TF.ind]
    ind.both <- which(!is.na(data[,val.ind[1]]) & !is.na(data[,val.ind[2]]))
    tbl.obs.both <- as.data.frame(table(TF.both[ind.both]))
    names(tbl.obs.both) <- c("TF","both.count")

    TF <- c()
    if (length(tbl.miss)){    
        TF <- merge(tbl.obs.both, tbl.miss, by="TF")
    }else{
        TF <- tbl.obs.both
    }
    
    return(TF)
}
    
                                

summarize.diffs <- function(diffs){

    TFs <- unique(diffs$Regulator)
    Targs <- unique(diffs$Target)
    Targs <- Targs[!Targs %in% TFs]

    out <-  data.frame(gene = as.character(c(TFs, Targs)),
                       Regulator = 0)
    out$Regulator[out$gene %in% TFs] = 1

    for (value in c("arac","arac.missing","wgcna")){

        ## DETERMINE WHICH GENE ARE SIGN IN THE CURRENT ANALYSIS(VALUE) ##
        ind <- diffs[,names(diffs)==paste0(value,".chosen")]==1
        N.found <- sum(ind)
        print(paste0("Number of interesting ",value," found = ",N.found))        
        genes <- sort(table(diffs$Regulator[ind]), decreasing=T)               
        N.TFs <- length(genes)      
        
        print(paste0("Number of TFs identified = ",length(genes)))

        genes <- data.frame(genes)
        names(genes) <- c("gene", paste0("Reg.Count.",value))
        eval(parse(text=paste0("genes$",value,".sig <- 1")))
        
        out <- merge(out, genes, all=T, by = "gene")     

        genes <- data.frame(sort(table(diffs$Target[ind & (diffs$Target %in% TFs==F)]),
                                 decreasing=T))
        names(genes) <- c("gene", paste0("Targ.Count.",value))

        
        out <- merge(out, genes, all=T, by="gene")                   
    }
    
    out[is.na(out)] = 0
    ord <- order(out$Reg.Count.wgcna, decreasing=T)
    out <- out[ord,]    

    return(out)
}



get.bicmix.results <- function(){

    BMdir <- "/gpfs/data/stranger-lab/askol/TCGA/Expression_Analysis/BicMix/"

    f.edges <- read.table(file = paste0(BMdir,"Sex_F1_M0_spec_ensemble_nDup2_edges.csv"),
                          as.is=T, header=T, sep=";")
    m.edges <- read.table(file = paste0(BMdir,"Sex_F0_M1_spec_ensemble_nDup2_edges.csv"),
                          as.is=T, header=T, sep=";")
    
    f.nodes <- read.table(file = paste0(BMdir,"Sex_F1_M0_spec_ensemble_nDup2_nodes.csv"),
                          as.is=T, header=T, sep=";")
    m.nodes <- read.table(file = paste0(BMdir,"Sex_F0_M1_spec_ensemble_nDup2_nodes.csv"),
                          as.is=T, header=T, sep=";")
    
    f.edges <- merge(f.edges, f.nodes, by.x="Source", by.y="Id", all.x=T, all.y=F)
    f.edges$Source <- f.edges$Label
    f.edges <- f.edges[,-ncol(f.edges)]
    f.edges <- merge(f.edges, f.nodes, by.x="Target", by.y="Id", all.x=T, all.y=F)
    f.edges$Target <- f.edges$Label
    f.edges <- f.edges[,-ncol(f.edges)]
    
    m.edges <- merge(m.edges, m.nodes, by.x="Source", by.y="Id", all.x=T, all.y=F)
    m.edges$Source <- m.edges$Label
    m.edges <- m.edges[,-ncol(m.edges)]
    m.edges <- merge(m.edges, m.nodes, by.x="Target", by.y="Id", all.x=T, all.y=F)
    m.edges$Target <- m.edges$Label
    m.edges <- m.edges[,-ncol(m.edges)]
    
    return(list(f.bicmix = f.edges, m.bicmix=m.edges))
    
}
    
    
test.clinical.significance <- function(genes,
                                       tissue="TCGA-LIHC",
                                       print.it = T,
                                       plot.name="Surv_qq_plots.pdf"){

    ## genes <- unique(summ$gene);  tissue="TCGA-LIHC"; plot.name="Surv_qq_plots.pdf"
    
    DIR <- "/gpfs/data/stranger-lab/askol/TCGA/"
    META.DIR <- "/gpfs/data/stranger-lab/askol/TCGA/TCGA_MetaData/"

    keep.col <- c("hits.cases.case_id",
                  "hits.cases.demographic.gender",
                  "hits.cases.diagnoses.days_to_death",
                  "hits.cases.diagnoses.days_to_last_follow_up",
                  "hits.cases.diagnoses.age_at_diagnosis",
                  "hits.cases.demographic.race",
                  "hits.cases.demographic.ethnicity")
    meta <- read.table(file <- paste0(META.DIR, tissue,".METAdata.txt"),
                       header=T, as.is=T, sep="\t", quote="")
    meta <- meta[,keep.col]
    names(meta) <- gsub("hits\\..+\\.","",names(meta))

    meta$case_id <- paste0("A",meta$case_id)
    meta <- subset(meta, case_id %in% colnames(count.auto$counts))

    ## REMOVE DUPES ##
    ind <- which(duplicated(meta$case_id))
    meta <- meta[-ind,]

    ## REMOVE DATA POINTS WITH NO TIME OF EVENT ##
    ind <- which(is.na(meta$days_to_last_follow_up) &
                 is.na(meta$days_to_death))
    meta <- meta[-ind,]
    ## FOR THE SUBSET OF SAMPLES THAT HAVE BOTH DAYS TO DEATH AND TO LAST FOLLOW-UP
    ## ALL DEATH DAYS > FOLLOW-UP, SO WILL USE DEATH DAYS ##
    meta$event.days <- apply(meta[,c("days_to_last_follow_up", "days_to_death")], 1,
                             max, na.rm=T)
    cens.ind <- which(is.na(meta$days_to_death)==T & is.na(meta$days_to_last_follow_up)==F)
    meta$censored = 2
    meta$censored[cens.ind] = 1    
    
    ## GET LIMMA/VOOM RESULTS ##
    ## source("../Code/process_expression_data_funcs.r")
    ## DE_by_covariate(tissue) 
    
    count.auto <- readRDS(
        paste0(DIR,"/TCGA_Expression/Limma_Voom_Use/TCGA-LIHC_mrna_autosome_DGEList.RDS"))
    genes.all <- unique(count.auto$genes$ENSEMBL)
    lcpm <- cpm(count.auto, log=TRUE, prior.count=3)
    keep <- which(rownames(lcpm) %in% genes)
    lcpm <- as.data.frame(t(lcpm[keep,]))
    lcpm$id <- rownames(lcpm)
       
    data <- merge(meta[,-c(3,4)], lcpm, by.x="case_id", by.y="id", all.y=T)
    
    file = paste0("Plots/",plot.name)
    if (print.it){
        pdf(file = file)        
        sfit <- survfit(Surv(event.days, censored)~gender, data=data)
        print(ggsurvplot(sfit, data=data))
    }
    
    out <- c()
    for (gene in genes){

        data$exp <- data[,which(names(data)==gene)]

        sfit.all <- coxph(Surv(event.days, censored)~exp, data=data)
        sfit.male <- coxph(Surv(event.days, censored)~exp,
                           data=subset(data, gender=="male"))
        sfit.female <- coxph(Surv(event.days, censored)~exp,
                             data=subset(data, gender=="female"))

        out <- rbind(out, c(gene, sfit.all$coef, sfit.male$coef, sfit.female$coef,
                            summary(sfit.all)$logtest[c(1,3)],
                            summary(sfit.male)$logtest[c(1,3)],
                            summary(sfit.female)$logtest[c(1,3)]))
        
    }
    out <- as.data.frame(out)
    out[,-1] <- apply(out[,-1], 2, function(x) as.numeric(as.character(x)))
    out[,2:4] <- exp(out[,2:4])
    names(out) <- c("gene", paste0("HR.",c("all","male","female")),
                    paste0(c("LR.test.","LR.p."), "all"),
                    paste0(c("LR.test.","LR.p."), "male"),
                    paste0(c("LR.test.","LR.p."), "female"))

    ##      
    if (print.it){
        print(qq_plot(out$LR.p.all, title = "QQ plot males and females combined"))
        print(qq_plot(out$LR.p.male, title = "QQ plot males"))
        print(qq_plot(out$LR.p.female, title = "QQ plot females"))
        
        p <- ggplot(out, aes(-log10(LR.p.male), -log10(LR.p.female))) + geom_point()    
        print(p)
        
        p <- ggplot(out, aes(HR.male, HR.female)) + geom_point()
        print(p)
        
        dev.off()    
        print(paste("Created Plot ",file))
        
        ## PLOT TOP HITS ##
        plot.surv(out, data)
    }
    
    ## RETURN A SET OF THE SAME NUMBER OF RANDOME GENES ##
    genes.rand <- sample(genes.all, nrow(out))
    return(list(out=out, genes.rand=genes.rand))
}

plot.surv <- function(anly.info, data, p.cutoff = 10^-4, plot.name = "survPlots.pdf"){

    file <- paste0("Plots/survPlots.pdf")
    pdf(file=file)
    genes <- anly.info$gene[anly.info$LR.p.male <= p.cutoff]
    
    for (gene in genes){

        data$exp <- data[,which(names(data)==gene)]
        data.male <- subset(data, gender=="male")
        q <- quantile(data.male$exp, probs=c(0,.2,.4,.6,.8,1), names=T)
        data.male$q <- cut(data.male$exp, q)        
        levels(data.male$q) <- paste0("Q",1:5)
        data.male$q <- as.character(data.male$q)
        sfit.male <- survfit(Surv(event.days, censored)~q, data=data.male)

        print(ggsurvplot(sfit.male, data=data.male, title = gene))
        
    }
    dev.off()
    print(paste0("Wrote file ",file))
    
}


run.gsea.all <- function(tissue, RankDir, GSEAOutDir){

    geneSet <- "msigdb.v6.2.symbols.gmt"
    rankFiles <- dir(pattern=".rnk", RankDir)
                     
    rm.ind <- grep("a\\.|w\\.|am\\.", rankFiles)
    rankFiles <- rankFiles[-rm.ind]
    
    for (rankFile in rankFiles){                  
        
        print(paste("Working on ",geneSet," and ",rankFile))
        
        dir <- gsub("\\.rank.+","",rankFile)
        geneSetName <- gsub("\\.v6.+","",geneSet)
        rpt.lab <- paste0(dir,"_",geneSetName)
        params <- rbind(c("rnk",
                          paste0(RankDir, rankFile)),
                        c("out",
                          paste0(GSEAOutDir, rpt.lab)))
        param.file <- paste0("~/GSEA_parameters_", gsub("-","_",tissue), ".txt")
        write.table(file = param.file, params, quote=F, col.names=F,
                    row.names=F, sep="\t")
        cmd <- paste0("java -cp /home/askol/bin/gsea-3.0.jar -Xmx5000m ",
                      "xtools.gsea.GseaPreranked ",
                      "-param_file ", param.file,
                      " -gmx /home/askol/bin/GSEA_genesets/",
                      "msigdb_v6.2_GMTs/",geneSet," -norm meandiv -nperm 1000 ",
                      "-scoring_scheme weighted -rpt_label ", rpt.lab,
                      " -create_svgs false -make_sets true -plot_top_x 100 ",
                      "-rnd_seed timestamp ",
                      "-set_max 500 -set_min 15 -zip_report false  -gui false")
        print(cmd)
        system(cmd)
    }
}

run.gsea.custom <- function(tissue, RankDir, GSEAOutDir){

    geneSet <- "Hormone_Immune_Custom.gmt"
    rankFiles <- dir(pattern=".rnk", RankDir)
    
    rm.ind <- grep("a\\.|w\\.|am\\.", rankFiles)
    rankFiles <- rankFiles[-rm.ind]
    
    for (rankFile in rankFiles){                  
        
        print(paste("Working on ",geneSet," and ",rankFile))
        
        dir <- gsub("\\.rank.+","",rankFile)
        geneSetName <- gsub("\\.gmt","",geneSet)
        rpt.lab <- paste0(dir,"_",geneSetName)
        params <- rbind(c("rnk",
                          paste0(RankDir, rankFile)),
                        c("out",
                          paste0(GSEAOutDir, rpt.lab)))
        param.file <- paste0("~/GSEA_parameters_", gsub("-","_",tissue), ".txt")
        write.table(file = param.file, params, quote=F, col.names=F,
                    row.names=F, sep="\t")
        cmd <- paste0("java -cp /home/askol/bin/gsea-3.0.jar -Xmx5000m ",
                      "xtools.gsea.GseaPreranked ",
                      "-param_file ", param.file,
                      " -gmx /home/askol/bin/GSEA_genesets/Custom/",geneSet,
                      " -norm meandiv -nperm 1000 ",
                      "-scoring_scheme weighted -rpt_label ", rpt.lab,
                      " -create_svgs false -make_sets true -plot_top_x 100 ",
                      "-rnd_seed timestamp ",
                      "-set_max 500 -set_min 15 -zip_report false  -gui false")
        print(cmd)
        system(cmd)
    }
}
    
process.gsea.output <- function(tissue, RankDir, GSEAOutDir, qThresh=0.1){

    ## READ IN RESULTS AND SORT BY NORMALIZED ES ##

    geneSet <- "msigdb.v6.2.symbols.gmt"

    rankFiles <- dir(pattern=".rnk", RankDir)
    rm.ind <- grep("a\\.|w\\.|am\\.", rankFiles)
    rankFiles <- rankFiles[-rm.ind]

    gsea <- c()
    genesets <- list()
    gseaTbl <- c()

    print(paste0("Working on tissue: ",tissue))
    
    for (rankFile in rankFiles){
        
        for (sign in c("pos","neg")){
            
            if (sign=="neg" &
                length(grep("mean",rankFile) == 0)){
                next
            }
            
            print(paste0("Working on rankFile: ",rankFile))
            
            dir <- gsub("\\.rank.+","",rankFile)
            geneSetName <- gsub("\\.v6.+","",geneSet)
            prefix <- paste0(dir,"_",geneSetName)
            
            out.dir <- paste0(GSEAOutDir, prefix,"/")
            
            dirs <- list.dirs(out.dir)
            files <- dir(dirs[2], pattern="gsea")
            files <- files[grep("xls",files)]
            file <- paste0(dirs[2],"/",files[grep(sign,files)])
            if (length(files) == 0){
                print(paste0("Didn't find file: ",file))
                next
            }
            out <- read.table(file = file, as.is=T, header=T, sep="\t")
            if (nrow(out) == 0){
                print(paste0("No data found in file ",file))
                next
            }
            nr <- nrow(out)
            qpass.ind <- which(out$FDR.q.val < .25)
            if (length(qpass.ind)>0){
                gsea <- rbind(gsea, cbind(dir, geneSetName, out[qpass.ind,]))
            }
            
        }
    }    

    rm.cols <- which(names(gsea) %in% c("GS.DETAILS","X"))
    gsea <- gsea[, -rm.cols]

    ## REMOVE SETS WITH FDR.P.VAL == NA ##
    rm.ind <- which(is.na(gsea$FDR.q.val))
    if (length(rm.ind) > 0){
        gsea <- gsea[-rm.ind,]
    }

    ## CREATE A KEEP INDEX
    gsea$keep = 1

    ## KEEP ONLY GENESETS WITH FDR < 10% ##
    lose.ind <- which(gsea$FDR.q.val > qThresh)
    if (length(lose.ind) > 0){
        gsea$keep[lose.ind] <- 0
    }

    gsea <- gsea[order(gsea$FDR.q.val),]
    out.file <- paste0(GSEAOutDir, tissue, "_GSEA_summary.txt")
    
    write.table(file = out.file, gsea, quote=F, row.names=F, col.names=T, sep="\t")
    
}


  
process.gsea.output.custom <- function(tissue, RankDir, GSEAOutDir, qThresh=0.25){

    ## READ IN RESULTS AND SORT BY NORMALIZED ES ##

    geneSet <- "Hormone_Immune_Custom.gmt"

    rankFiles <- dir(pattern=".rnk", RankDir)
    rm.ind <- grep("a\\.|w\\.|am\\.", rankFiles)
    rankFiles <- rankFiles[-rm.ind]

    gsea <- c()
    genesets <- list()
    gseaTbl <- c()

    print(paste0("Working on tissue: ",tissue))
    
    for (rankFile in rankFiles){
        
        for (sign in c("pos","neg")){
            
            if (sign=="neg" &
                length(grep("mean",rankFile) == 0)){
                next
            }
            
            print(paste0("Working on rankFile: ",rankFile))
            
            dir <- gsub("\\.rank.+","",rankFile)
            geneSetName <- gsub("\\.gmt+","",geneSet)
            prefix <- paste0(dir,"_",geneSetName)
            
            out.dir <- paste0(GSEAOutDir, prefix,"/")
            
            dirs <- list.dirs(out.dir)
            dirs <- dirs[-grep("error", dirs)]
            files <- dir(dirs[2], pattern="gsea")
            files <- files[grep("xls",files)]
            file <- paste0(dirs[2],"/",files[grep(sign,files)])
            if (length(files) == 0){
                print(paste0("Didn't find file: ",file))
                next
            }
            out <- read.table(file = file, as.is=T, header=T, sep="\t")
            if (nrow(out) == 0){
                print(paste0("No data found in file ",file))
                next
            }
            nr <- nrow(out)
            qpass.ind <- which(out$FDR.q.val < .25)
            if (length(qpass.ind)>0){
                gsea <- rbind(gsea, cbind(dir, geneSetName, out[qpass.ind,]))
            }
            
        }
    }    

    rm.cols <- which(names(gsea) %in% c("GS.DETAILS","X"))
    gsea <- gsea[, -rm.cols]

    ## REMOVE SETS WITH FDR.P.VAL == NA ##
    rm.ind <- which(is.na(gsea$FDR.q.val))
    if (length(rm.ind) > 0){
        gsea <- gsea[-rm.ind,]
    }

    ## CREATE A KEEP INDEX
    gsea$keep = 1

    ## KEEP ONLY GENESETS WITH FDR < 25% ##
    lose.ind <- which(gsea$FDR.q.val > qThresh)
    if (length(lose.ind) > 0){
        gsea$keep[lose.ind] <- 0
    }

    gsea <- gsea[order(gsea$FDR.q.val),]
    out.file <- paste0(GSEAOutDir, tissue, "_GSEA_custom_summary.txt")
    
    write.table(file = out.file, gsea, quote=F, row.names=F, col.names=T, sep="\t")
    
}


