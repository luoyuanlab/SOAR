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
pathways <- fread(paste0(out_dir, "/pathways_step1_success.txt"), header = FALSE)$V1
pathways_success <- character(0)

for (pathway in pathways) {
	# pathway = "CXCL"
	pathway_dir <- paste0(out_dir, "/", pathway)
	
	possibleError <- tryCatch( {
		sce <- readRDS(paste0(pathway_dir, "/step1_sce.RDS"))
		assoRes <- fread(paste0(pathway_dir, "/assoRes.txt"), sep = "\t")
		genes <- assoRes$V1
		assoRes$V1 <- NULL
		assoRes <- data.frame(assoRes)
		rownames(assoRes) <- genes
		
		source(paste0(pathway_dir, "/step3.R"))
		fwrite(yhatScaled, paste0(pathway_dir, "/yhatScaled.txt"), sep = "\t", row.names = TRUE)
	},
		error=function(e) e
	)
	if(!inherits(possibleError, "error")) {
		pathways_success <- c(pathways_success, pathway)
	}
}

pathways_failed <- pathways[!pathways %in% pathways_success]
fwrite(data.frame(pathways_success), paste0(out_dir, "/pathways_step3_success.txt"), col.names = FALSE)
fwrite(data.frame(pathways_failed), paste0(out_dir, "/pathways_step3_failed.txt"), col.names = FALSE)
