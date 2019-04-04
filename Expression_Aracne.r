AracDir <- paste0("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/",
                  "Results/Aracne/")
DataDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Data/Aracne/"
ResultsDir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/Aracne/"

### DETEMRINE ARGUEMENT FOR ANALYSIS ###
###
Expression_Aracne <- function(tissue, nseeds=50){

    seeds.p.job = 5

    arac.dir <- AracDir
    ## CHECK IF DATA ALREADY EXISTS FOR TISSUE ##
    tmp <- file.count(tissue, arac.dir)
    number.seed.files <- tmp$file.count
    if (sum(number.seed.files) > 0){
        print(paste("Previous bootstrapNetwork files found for ", tissue))
        print("Removing files . . .")
        for (i in 1:2){
            cmd <- paste0("rm ", tmp$dirs[i],"bootstrapNetwork*")
            system(cmd)
        }
    }
    
    ## CREATE FILES and DETERMINE THRESHOLDS ##
    data.dir <- DataDir
    out.dir <- paste0(ResultsDir ,tissue,"/")

    ## MAKE ARACNE FILES
    out.dirs <- make.aracne.files(tissue, nseeds, data.dir, out.dir)

    sex.inds <- c(1,2)
    n.files = 0
    while (n.files < 2){

        Sys.sleep(30) ## wait 60 seconds between seeing if threshold files exist ##
        m.file <- length(dir(out.dirs[1], pattern="miThreshold"))>0
        f.file <- length(dir(out.dirs[2], pattern="miThreshold"))>0

        n.files <- m.file + f.file

        print(n.files)
    }
    
    for (sex.ind in 1:2){
        
        ## NO DPI VERSION ##
        dir.out <- out.dirs[sex.ind]        
        runAracne(nseeds = nseeds, seeds.p.job=seeds.p.job, dir.out=dir.out)
        
    }
    
    ## CONSOLIDATE BOOTSTRAP OUTPUT ##
    
    sex.inds <- c(1,2)
    completed.seeds <- c(0,0)
    while (length(sex.inds) > 0){

        Sys.sleep(60)
        for (sex.ind in sex.inds){
            
            dir.out <- out.dirs[sex.ind]
            lns <- length(dir(dir.out, pattern="bootstrapNetwork"))

            if (lns > completed.seeds[sex.ind]){
                completed.seeds[sex.ind] <- lns
                print(paste0("Number of completed seeds for sex.ind ",sex.ind," = ", lns))
            }
            
            ## WAIT UNTIL ALL BOOTSTRAPNETWORK FILES EXIST BEFORE CONSOLIDATING ##
            if (lns == nseeds){                

                print(paste0("Consolidating for sex.ind ",sex.ind))
                runAracneConsolidate(dir.out)
                sex.inds <- sex.inds[-which(sex.inds == sex.ind)]
                print(paste0("Consolidation complete for sex.ind ",sex.ind))
            }
            
        }

    }
    
    print(paste0("Completed Expression_Aracne for tissue: ", tissue))
}

## MAIN ##

source("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Code/Expression_Aracne_funcs.r")
args <- commandArgs(TRUE)
tissue = args[1]
Expression_Aracne(tissue)
print(date())

finish.file = paste0(tissue,".finished")
system(paste0("touch ",finish.file))
