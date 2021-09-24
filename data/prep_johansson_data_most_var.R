# Rscript prep_johansson_data_most_var.R

nMostVar <- 1000

inFile <- "johansson_data_relative_ratios_to_pool.csv"


### prep count
all_dt <- read.delim(inFile, sep=",", header=T)

stopifnot(!duplicated(all_dt$gene_symbol))

tmp_dt <- all_dt
rownames(tmp_dt) <- tmp_dt$gene_symbol
tmp_dt$gene_symbol <- tmp_dt$ensembl_id <- NULL

all_vars <- apply(tmp_dt, 1, var)
stopifnot(!is.na(all_vars))

all_vars_s <- sort(all_vars, decreasing=TRUE)


to_keep <- names(all_vars_s)[1:nMostVar]

stopifnot(to_keep %in% all_dt$gene_symbol)

out_dt <- all_dt[all_dt$gene_symbol %in% to_keep,]

stopifnot(ncol(out_dt) == ncol(all_dt))
stopifnot(setequal(out_dt$gene_symbol, to_keep))
stopifnot(nrow(out_dt) == nMostVar)

outFile <-  paste0("johansson_data_relative_ratios_to_pool_", nMostVar, "mostVar.txt")
write.table(out_dt, file = outFile, sep=",",
            col.names=TRUE, row.names=FALSE, quote=F)

cat(paste0("... written: ", outFile, "\n"))