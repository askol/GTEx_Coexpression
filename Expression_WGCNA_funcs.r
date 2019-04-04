library(edgeR)
library(limma)
library(WGCNA)
library(reshape2)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()
DataDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Data/WGCNA/"

get.data <- function(tissue, DataDir){
    
    data.file <- paste0(DataDir, tissue,"_lcpm_auto_male.txt")
    lcpm.m <- read.table(data.file, as.is=T, header=T, row.names=1)

    data.file <- paste0(DataDir, tissue,"_lcpm_auto_female.txt")
    lcpm.f <- read.table(data.file, as.is=T, header=T, row.names=1)
    
    setLabels = shortLabels = c("Male", "Female")
    
    ## Form multi-set expression data:
    ## columns starting from 9 contain actual expression data.
    
    multiExpr = vector(mode = "list", length = 2)
    
    multiExpr[[1]] = list(data = as.data.frame(t(lcpm.m)))
    names(multiExpr[[1]]$data) = rownames(lcpm.m)
    rownames(multiExpr[[1]]$data) = colnames(lcpm.m)
    
    multiExpr[[2]] = list(data = as.data.frame(t(lcpm.f)))
    names(multiExpr[[2]]$data) = rownames(lcpm.f)
    rownames(multiExpr[[2]]$data) = colnames(lcpm.f)

    ## Check that the data has the correct format for many functions
    ## operating on multiple sets:
    exprSize = checkSets(multiExpr)

    ## Check that all genes and samples have sufficiently low numbers of missing values.
    gsg = goodSamplesGenesMS(multiExpr, verbose = 3);

    ## REMOVE GENES IF LOW VARIANCE OR TOO MUCH MISSINGNESS ##
    
    if (!gsg$allOK){
        ## Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
            printFlush(paste("Removing genes:",
                             paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                   collapse = ", ")))
        for (set in 1:exprSize$nSets){
            if (sum(!gsg$goodSamples[[set]]))
                printFlush(paste("In set", setLabels[set], "removing samples",
                                 paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]],
                                       collapse = ", ")))
            ## Remove the offending genes and samples
            multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
        }
        ## Update exprSize
        exprSize = checkSets(multiExpr)
    }
    
    return(list(multiExpr = multiExpr, lcpm.f = lcpm.f, lcpm.m = lcpm.m))
}


find.soft.threshold <- function(multiExpr, ResultDir, powers){
    
    ## Choose a set of soft-thresholding powers
    powers = c(seq(4,10,by=1), seq(12,20, by=2));
    
    ## Initialize a list to hold the results of scale-free analysis
    powerTables = list()

    ## Call the network topology analysis function for each set in turn
    for (set in 1:2)
        powerTables[[set]] <- list(data = pickSoftThreshold(multiExpr[[set]]$data,
                                       powerVector=powers,
                                       verbose = 2)[[2]])
    
    collectGarbage()

    ## ############## ##
    ##  PLOT RESULTS  ## 

    plot.softthresh(powerTables, ResultDir, powers)
     

    ## IDENTIFY BEST SOFT.THRESHOLD ##
    best.power <- select.best.power(powerTables)

    return(best.power)
    

}

select.best.power <- function(powerTables){

    best.ind <- c(0,0)    
    
    for (i in 1:2){

        x = powerTables[[i]]$data$Power
        y = powerTables[[i]]$data$SFT.R.sq

        deltas <- y[-1] - y[-length(y)]
        perc.max <- y/max(y)
        r.gt.8 <- y> .8
        pm.9 <- perc.max > .9

        best.ind[i] <- min(which(r.gt.8 & pm.9))

    }

    print(paste("Best powers for males and females found: ",x[best.ind]))
    print(paste("Using power of ",x[min(best.ind)]," for both sexes"))
    return(best.power = x[min(best.ind)])
                    
}


plot.softthresh <- function(powerTables, out.dir, powers){

    file <- paste0(out.dir,"/SoftThreshold.pdf")
    pdf(file = file)
    
    colors = c("black", "red")
    setLabels = shortLabels = c("Male", "Female")
    
    ## Will plot these columns of the returned scale free analysis tables    
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity",
        "Median connectivity", "Max connectivity")
    
    ## Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:2){
        for (col in 1:length(plotCols)){
            
            ylim[1, col] = min(ylim[1, col],
                    powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
            ylim[2, col] = max(ylim[2, col],
                    powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
        }
    }

    ## Plot the quantities in the chosen columns vs. the soft thresholding power
    sizeGrWindow(8, 6)
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:2){
        
        if (set==1){
            
            plot(powerTables[[set]]$data[,1],
                 -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                 xlab="Soft Threshold (power)",ylab=colNames[col],type="n",
                 ylim = ylim[, col],
                 main = colNames[col])
            addGrid();
        }
        if (col==1){
            
            text(powerTables[[set]]$data[,1],
                 -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                 labels=powers,cex=cex1,col=colors[set])
        } else {
            text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
                 labels=powers,cex=cex1,col=colors[set]);
            if (col==1){
                
                legend("bottomright", legend = setLabels, col = colors, pch = 20)
            } else {
                legend("topright", legend = setLabels, col = colors, pch = 20) 
            }
        }                
    }

    dev.off()
    ## DONE PLOTTING ##

    print(paste0("Data used for determined soft threshold in ",file))

}
        



calculate.TOM <- function(multiExpr, softPower){

    
    exprSize <- checkSets(multiExpr)
    nSets <- exprSize$nSets; nGenes <- exprSize$nGenes
##    softPower = 5;
    ## Initialize an appropriate array to hold the adjacencies
    adjacencies = list()

    print("Calculating adjacencies. . .")

    ## Calculate adjacencies in each individual data set
    for (set in 1:nSets)
        adjacencies[[set]] = abs(cor(multiExpr[[set]]$data,
                       use = "p"))^softPower
    

    print("Calculating TOMs. . .")
    TOM = list()

    ## Calculate TOMs in each individual data set
    for (set in 1:nSets)
        TOM[[set]] = TOMsimilarity(adjacencies[[set]])

    
    ## Define the reference percentile
    scaleP = 0.95
    
    ## Set RNG seed for reproducibility of sampling
    set.seed(12345)
    
    ## Sample sufficiently large number of TOM entries
    nSamples = as.integer(1/(1-scaleP) * 1000)
    
    ## Choose the sampled TOM entries
    scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
    TOMScalingSamples = list()
    
    ## These are TOM values at reference percentile
    scaleQuant = rep(1, nSets)
    
    ## Scaling powers to equalize reference TOM values
    scalePowers = rep(1, nSets)

    print("Scaling TOMs. . .")
    ## Loop over sets
    for (set in 1:nSets){
        
        ## Select the sampled TOM entries
        TOMScalingSamples[[set]] = as.dist(TOM[[set]])[scaleSample]

        ## Calculate the 95th percentile
        scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8)

        ## Scale the male TOM
        if (set>1){
  
            scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
            TOM[[set]] = TOM[[set]]^scalePowers[set];
        }
    }

    ## ADD GENE NAMES ##
    for (set in 1:nSets){
        rownames(TOM[[set]]) <- names(multiExpr[[set]]$data)
        colnames(TOM[[set]]) <- names(multiExpr[[set]]$data)
    }
    
    return(TOM)
}



output.TOMs <- function(TOMs, TFs, out.dir){

    Sexes <- c("Male", "Female")
    for (sex.ind in 1:2){

        out <- TOMs[[sex.ind]]
        tf.ind <- which(rownames(out) %in% TFs)

        out <- out[tf.ind,]
        out <- as.data.frame(out)
        out$TF <- rownames(out)
        
        out <- melt(out, id = "TF", variable.name = "Target", value.name="TOM")
        
        file.name <- paste0(out.dir, Sexes[sex.ind],"/TOM_table.txt")

        write.table(file = file.name, out, quote=F, row.names=F, col.names=T)

        print(paste0("Wrote output to ",file.name))

    }
}
