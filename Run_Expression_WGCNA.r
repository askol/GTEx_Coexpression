setwd("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/WGCNA/")

create.pbs <- function(tissue, dir){

    file <- paste0(dir,tissue,"/",tissue,".pbs")
    sh.txt <- rbind(
        "#!/bin/bash",
        "#PBS -l  nodes=1:ppn=2,mem=16gb",
        "#PBS -l walltime=96:00:00",
        paste0("#PBS -o ",dir),
        "#PBS -j oe",
        paste0("#PBS -N WGC.",tissue),
        "module load gcc/5.4.0",
        "module load R",
        paste0("R CMD BATCH  '--args ",tissue,
               "' /gpfs/data/stranger-lab/askol/GTEx/Coexpression/Code/Expression_WGCNA.r ",
               dir, "WGCNA_",tissue,"_pbs.out")
        )

    write.table(file = file, sh.txt, quote=F, row.names=F, col.names=F)
    system(paste0("chmod +x ",file))
    
    return(file)
}

tissues <- read.table(file = paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "data/support_files/all_v8_tissues_both_sexes.txt"), as.is=T)
tissues <- unlist(tissues)

for (tissue in tissues){

    out.dir <- paste0("/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/WGCNA/",tissue,"/")
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
     d <-paste0(out.dir,"Plots/")
    if (!file.exists(d)){
        dir.create(d, recursive=T)
    } 

    out.dir <- "/gpfs/data/stranger-lab/askol/GTEx/Coexpression/Results/WGCNA/"
    
    file <- create.pbs(tissue, out.dir)

    system(paste0("qsub ",file))
}
