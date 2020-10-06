library(dplyr)
library(tidyr)
library(parallel)
library(pheatmap)
library(RColorBrewer)

max_legend <- function(a) {
    cond <- TRUE
    n <- 0
    while (cond) {
        n <- n + 1
        b <- a %/% 2
        if (b == 1) {
            cond <- FALSE
            return(n + 1)
        } else {
            a <- b
        }
    }
}

fill_mat <- function(chrom) {
    input_fil <- filter(input, X1 == chrom)
    input_fil$X2 <- input_fil$X2 %/% bin_size + 1
    input_fil <- unite(input_fil, X1, c(X1, X2), sep = "-")
    mat_temp <- mat
    for (i in seq_len(nrow(input_fil))) {
        mat_temp[input_fil$X3[i], input_fil$X1[i]] <- mat_temp[input_fil$X3[i], input_fil$X1[i]] + 1
    }
    return(mat_temp)
}

chrom_size_list <- c(23513712, 25286936, 28111227, 32079331, 1348131, 23648458)
chrom_list <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
bin_size <- 100000

input_path <- "/data/agl_data/ningjun/RNA-DNA/GSE97131/GSE97131_mRNA.txt"
output_tiff <- "/data/agl_data/ningjun/RNA-DNA/GSE97131/all-mRNA.tiff"
output_pdf <- "/data/agl_data/ningjun/RNA-DNA/GSE97131/all-mRNA.pdf"
output_csv <- "/data/agl_data/ningjun/RNA-DNA/GSE97131/all-mRNA.csv"


cat("-------------- loading file -------------", "\n")
input <- readr::read_delim(input_path, delim = "\t", col_names = FALSE)
input <- arrange(input, desc(X4), desc(X5))

cat("------------ creating matrix ------------", "\n")
input <- unite(input, X3, c(X4, X3),sep = "-")
mat_row <- unique(input$X3)
mat_col <- vector()
for (i in seq_len(length(chrom_list))) {
    mat_col <- append(mat_col, paste0(chrom_list[i], "-", seq(1, chrom_size_list[i] %/% bin_size + 1)))
}
mat <- matrix(0, length(mat_row), length(mat_col))
rownames(mat) <- mat_row
colnames(mat) <- mat_col
cat("----------- size (", dim(mat), ") ------------", "\n")

cat("---------- fullfilling matrix -----------", "\n")
core <- detectCores() %/% 3
cl <- makeForkCluster(core)
mat_list <- parLapply(cl, chrom_list, fill_mat)
mat_out <- Reduce("+", mat_list)
stopCluster(cl)

cat("--------------- making csv ---------------", "\n")
write.csv(mat_out, output_csv)

cat("------------- making heatmap ------------", "\n")
col <- colorRampPalette(brewer.pal(9, "YlOrRd"))
breaks <- c(0, 2^(0:max_legend(max(mat_out))))
pic <- pheatmap(mat_out, cluster_rows = FALSE, cluster_cols = FALSE,
			 	show_rownames = TRUE, show_colnames = TRUE,
                cellheight = 1, cellwidth = 5, fontsize = 1,
				col = col(length(breaks)), breaks = breaks,
                legend = FALSE, border_color = NA,
                filename = output_pdf)

cat("----------------- Done! -----------------", "\n")