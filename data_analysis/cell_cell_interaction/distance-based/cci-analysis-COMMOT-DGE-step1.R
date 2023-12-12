library(stringr)
library(data.table)
library(tradeSeq)
library(clusterExperiment)

args <- commandArgs(trailingOnly=TRUE)

sample_dir <- args[1]
thr_type <- args[2]
# sample_dir <- "/projects/b1131/SpatialT/10x/PID1/DS1D/DS1D.1"
# thr_type <- "short"

out_dir <- paste0(sample_dir, "/analysis/Distance/COMMOT_dec/", thr_type)
pathways <- fread(paste0(out_dir, "/pathways.txt"), header = FALSE)$V1
pathways_success <- character(0)

for (pathway in pathways) {
	# pathway = "MIF"
	pathway_dir <- paste0(out_dir, "/", pathway)
	
	possibleError <- tryCatch( {
		X <- fread(paste0(pathway_dir, "/step1_X.csv"), header = TRUE)
		pseudoTime <- fread(paste0(pathway_dir, "/step1_pseudoTime.csv"), header = FALSE)$V1
		cellWeight <- fread(paste0(pathway_dir, "/step1_cellWeight.csv"), header = FALSE)$V1
		spot_ids <- X$V1
		X$V1 <- NULL
		X <- t(as.matrix(X))
		colnames(X) <- spot_ids
		
		source(paste0(pathway_dir, "/step1.R"))
		fwrite(assoRes, paste0(pathway_dir, "/assoRes.txt"), sep = "\t", row.names = TRUE)
		saveRDS(sce, paste0(pathway_dir, "/step1_sce.RDS"))
	},
		error=function(e) e
	)
	if(!inherits(possibleError, "error")) {
		pathways_success <- c(pathways_success, pathway)
	}
}

pathways_failed <- pathways[!pathways %in% pathways_success]
fwrite(data.frame(pathways_success), paste0(out_dir, "/pathways_step1_success.txt"), col.names = FALSE)
fwrite(data.frame(pathways_failed), paste0(out_dir, "/pathways_step1_failed.txt"), col.names = FALSE)
