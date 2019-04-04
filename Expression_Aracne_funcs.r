library(readr)


make.aracne.files <- function(tissue, nseeds, data.dir, out.dir){

    sexes <- c("Male", "Female")
    
    print(paste0("Making scripts for project ", tissue))

    dir.out.male <-  paste0(out.dir, "Male/")
    dir.out.female <-  paste0(out.dir, "Female/")
    out.dirs <- c(dir.out.male, dir.out.female)
    
    file.expression.male <- paste0(data.dir, tissue, "_lcpm_auto_male.txt")
    file.expression.female <- paste0(data.dir, tissue, "_lcpm_auto_female.txt")    
    expression.files <- c(file.expression.male, file.expression.female)

    n.male <- as.numeric(system(paste0("head -n1 ",file.expression.male,
                                       " | wc -w"), intern=TRUE))
    n.female <- as.numeric(system(paste0("head -n1 ", file.expression.female,
                                         " | wc -w"), intern=TRUE))
    Ns <- c(n.male, n.female)
    
    ps <- rep(10^-8, 2)
    min.pos <- order(Ns)[1]
    max.pos <- 3-min.pos
    
    p.alt <- 0.25
    ps[min.pos] <- p.alt
    ps <- as.character(ps)
    ps <- gsub("e","E",ps)
    ps <- gsub("-0","-", ps)
    
    for (sex.ind in 1:2){

        file.expression <- expression.files[sex.ind]
        dir.out <- out.dirs[sex.ind]
        
        make.all.files(tissue, sex.ind, file.expression, dir.out, ps, nseeds)

        print(paste0("Determining threshold value for ", sexes[sex.ind]))
        runAracneThresh(dir.out)
    }

    return(out.dirs)
}

make.all.files <- function(tissue, sex.ind, file.expression, dir.out, ps, nseeds=50) {

    ## MAKE DIRECTORY IF IT DOES NOT ALREADY EXIST ##
    if (!dir.exists(dir.out)){
        dir.create(dir.out, recursive=T)
    }

    ## FEMALE IS CODED AS 2, MALE AS 1 ##
    sexes=c("Male","Female")     
    
    seeds <- round(runif(nseeds)/10^-8)
    file <- paste0(dir.out,"seeds.txt")
    write.table(file = file, seeds, quote=F, row.names=F, col.names=F)
    
    min.pos <- order(ps)[2]
    pboots <- rep(0.05, 2)
    pboots[min.pos] <- 0.50

    for (seed in 0:nseeds){
        
        if (seed == 0){

            ## DETERMINES THE APPOPRIATE MI THRESHOLD ##
            file = paste0(dir.out,"aracne_job_",seed,".sh")
            pre = prelude(tissue, seed, sex=sexes[sex.ind], dir.out)
            write.table(file = file, pre, quote=F, row.names=F, col.names=F)
            cmd <- paste0("java -Xmx5G -jar $ARACNEAP -e ", file.expression,
                          "  -o ", dir.out, " --pvalue ", ps[sex.ind],
                          " --consolidatepvalue ", pboots[sex.ind],

                          "/ARACNE/TFs_Ensembl_v_1.01.txt --seed 1 --calculateThreshold")
            write.table(file = file, cmd, append=T, quote=F, row.names=F, col.names=F)
                        
            system(paste0("chmod +x ",file))

            ## # WRITE CONSOLIDATION SCRIPT # ##

            file <- paste0(dir.out,"aracne_job_consolidate.sh")
            pre = prelude(tissue, "_CON", sex=sexes[sex.ind], dir.out)
            write.table(file = file, pre, quote=F, row.names=F, col.names=F)
            cmd <- paste("java -Xmx5G -jar $ARACNEAP --nodpi -o ",
                         dir.out, "--consolidate --consolidatepvalue ",pboots[sex.ind])
            write.table(file = file, cmd, append=TRUE, quote=F, row.names=F, col.names=F)
            system(paste0("chmod +x ",file))      
            
            next
        }
    
        if (seed%%5 == 1){

            ## NO DPI VERSION ##
            file <- paste0(dir.out,"aracne_jobs_",seed,"-",seed+4,".sh")
            sh.txt <- prelude(tissue, seed, sex=sexes[sex.ind], dir.out)
            write.table(file =file, sh.txt, quote=F, row.names=F, col.names=F)
        }

        cmd <-  paste0("java -Xmx5G -jar $ARACNEAP -e ",
                       file.expression, "  -o ", dir.out, "  --pvalue " ,ps[sex.ind],
                       " --consolidatepvalue ",pboots[sex.ind], 
                       " -t /gpfs/data/stranger-lab/askol/TCGA/Expression_Analysis/ARACNE/",
                       "TFs_Ensembl_v_1.01.txt --seed ", seeds[seed])
        write.table(file = file, cmd, append=T, quote=F, row.names=F, col.names=F)
        touch.file <- paste0(dir.out,"finished.",seed,".txt")
        cmd <- paste0("touch ",touch.file)
        write.table(file = file, cmd, append=T, quote=F, row.names=F, col.names=F)
        system(paste0("chmod +x ",file))
          
    }
}
        

runAracneThresh <- function(dir.out){
    file <- paste0(dir.out,"aracne_job_0.sh")
    cmd <- paste0("qsub ",file)
    print(paste0("Submitting ",file))
    system(cmd)
}


runAracne <- function(nseeds = 100, seeds.p.job = 5, dir.out){
    
    for (seed in 1:nseeds){

        if (seed%%seeds.p.job == 1){
            file <- paste0(dir.out,"aracne_jobs_",seed,"-",seed+4,".sh")
            cmd <- paste0("qsub ",file)
            print(paste0("Submitting ",file))
            system(cmd)
        }
    }
}

runAracneConsolidate <- function(dir.out){
    file <- paste0(dir.out,"aracne_job_consolidate.sh")    
    cmd <- paste0("qsub ",file)
    print(paste0("Submitting ",file))
    system(cmd)
}

        
prelude <- function(tissue, seed, sex="", dir){

    sh.txt <- rbind(
        "#!/bin/bash",
        "#PBS -l  nodes=1:ppn=1,mem=8gb",
        "#PBS -l walltime=96:00:00",
        paste0("#PBS -o ",dir),
        "#PBS -j oe",
        paste0("#PBS -N AR", tissue,"-", seed,"-",sex),
        paste("module load aracne-ap"))
    
    return(sh.txt)
}


file.count <- function(tissue, arac.dir){

    sexes <- c("Male","Female")
    arac.dir <- paste0(arac.dir, tissue,"/")
    arac.dirs <- paste0(arac.dir, sexes,"/")

    file.count <- c(0,0)
    for (sex.ind in c(1,2)){        
        
        file.count[sex.ind] <- length(dir(pattern="bootstrapNetwork", arac.dirs[sex.ind]))

    }
    return(list(file.count=file.count, dirs = arac.dirs))
}
  
