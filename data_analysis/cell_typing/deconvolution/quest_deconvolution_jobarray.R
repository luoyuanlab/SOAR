### Author: Yiming Li
### Example usage:
### quest_deconvolution_jobarray.R $sample_dir $ref_dir

library(stringr)
library(Seurat)
library(BayesPrism)
library(data.table)

# cd /share/fsmresfiles/SpatialT/10x/PID27/DS27A/DS27A_1160920F/processed
# conda activate R4
args <- commandArgs(trailingOnly=TRUE)

### Change the below parameters if needed
if (length(args) > 2) {
        n_cores <- as.integer(args[3])
} else {
        n_cores <- 20
}
chain.length <- 1000
burn.in <- 500
maxit <- 10000

st_dir <- "/projects/b1131/SpatialT"
sample_dir <- args[1]
ds_name <- str_split(sample_dir, '/')[[1]][7]
tech <- str_split(sample_dir, '/')[[1]][5]
p_name <- str_split(sample_dir, '/')[[1]][6]
ds_dir <- paste(c(st_dir, tech, p_name, ds_name), collapse = "/")
sample_name <- str_split(sample_dir, '/')[[1]][8]

### Assumes two files are present under this path:
### * sc.dat.filtered.pc.sig.RDS
### * cell_types.txt
ref_dir <- args[2]

##### Validate if reference files exist

if (!file.exists(paste0(ref_dir, "/cell_types.txt"))) {
	stop(paste0("[", ref_dir, "/cell_types.txt] does not exist\n"))
} else if (!file.exists(paste0(ref_dir, "/sc.dat.filtered.pc.sig.RDS"))) {
	stop(paste0("[", ref_dir, "/sc.dat.filtered.pc.sig.RDS] does not exist\n"))
}
final_ref <- readRDS(paste0(ref_dir, "/sc.dat.filtered.pc.sig.RDS"))
cell_types <- fread(paste0(ref_dir, "/cell_types.txt"))
cell.type.labels <- cell_types$label
cell.state.labels <- cell_types$label

##### Validate if transposed count matrix exist

if (!file.exists(paste0(sample_dir, "/processed/bk.dat.RDS"))) {
	stop(paste0("[", sample_dir, "/processed/bk.dat.RDS] does not exist\n"))
} else {
	bk.dat <- readRDS(paste0(sample_dir, "/processed/bk.dat.RDS"))
	cat(paste0("\n[", sample_name, "] - ST count matrix read from file\n"))
}



#####



### Create output dir
output_dir <- paste0(sample_dir, "/analysis")
if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}
output_dir <- paste0(sample_dir, "/analysis/deconvolution")
if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}

### Run BayesPrism
myPrism <- new.prism(
	reference = final_ref,
	mixture = bk.dat,
	input.type = "count.matrix", 
	cell.type.labels = cell.type.labels, 
	cell.state.labels = cell.state.labels,
	# key="tumor",
	key = NULL,
	outlier.cut = 0.01,
	outlier.fraction = 0.1
)
cat(paste0("\n[", sample_name, "] - deconvolution started [", Sys.time(), "]\n"))
bp.res <- run.prism(prism = myPrism, n.cores = n_cores, gibbs.control = list(burn.in = burn.in, chain.length = chain.length), opt.control = list(maxit = maxit))
saveRDS(bp.res, paste0(output_dir, "/BayesPrism_results.RDS"))
gc()
cat(paste0("\n[", sample_name, "] - deconvolution done [", Sys.time(), "]\n"))



#####



### Save thetas
theta.cv <- bp.res@posterior.theta_f@theta.cv
theta <- get.fraction(bp = bp.res, which.theta = "final", state.or.type = "type")
# BayesPrism advises to mask theta with CV above 0.5 (Visium)
theta[theta.cv > 0.5] <- 0
theta <- t(apply(theta, 1, function(x) x / sum(x)))
theta <- theta[,sort(colnames(theta))]
saveRDS(theta, paste0(output_dir, "/BayesPrism_theta.RDS"))

### Save deconvoluted cell-type-specific expressions
seurat_object <- readRDS(paste0(sample_dir, "/processed/Seurat.RDS"))
location <- fread(paste0(sample_dir, "/processed/coordinates.csv"))
meta <- seurat_object@meta.data
spot_id_mapping <- meta$new_spot_id
names(spot_id_mapping) <- rownames(meta)

deconv_genes <- colnames(bp.res@prism@mixture)
not_deconv_genes <- colnames(bk.dat)
not_deconv_genes <- not_deconv_genes[!not_deconv_genes %in% deconv_genes]
not_deconv_exp <- bk.dat[,not_deconv_genes]
not_deconv_exp <- not_deconv_exp / length(colnames(theta))

all_cell_types <- character(0)
ct_i <- 1
for (cell_type in colnames(theta)) {
	# cell_type <- "CAFs"
	cell_type_s <- gsub("/", ".", cell_type)
	cell_type_s <- gsub(" ", ".", cell_type_s)
	cell_type_s <- gsub("-", ".", cell_type_s)
	cell_type_s <- gsub("\\*", ".", cell_type_s)
	cell_type_s <- gsub("\\+", ".", cell_type_s)
	
	ct_exp <- get.exp(bp = bp.res, state.or.type = "type", cell.name = cell_type)
	ct_exp <- cbind(ct_exp, not_deconv_exp)
	
	counts_df <- data.table(ct_exp)
	counts_df$spot <- spot_id_mapping[rownames(ct_exp)]
	counts_df <- transpose(counts_df, keep.names = "gene", make.names = "spot")
	keep_spots <- intersect(colnames(counts_df), location$barcode)
	keep_spots <- c("gene", keep_spots)
	counts_df <- counts_df[,..keep_spots]
	
	#### Write counts, coordinates, and meta_spots to file
	fwrite(counts_df, paste0(output_dir, "/counts_", cell_type_s, ".csv"), sep = ",")
	all_cell_types[ct_i] <- cell_type_s
	ct_i <- ct_i + 1
}
write.table(all_cell_types, paste0(output_dir, "/all_cell_types.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
cat(paste0("\n[", sample_name, "] - cell-type-specific expression matrices saved [", Sys.time(), "]\n"))
