library(parallel)
library(dplyr)
library(tidyr)
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

fill_gap <- function(list_in) {
    list_out <- vector()
    list_out <- append(list_out, list_in[1])
    for (i in 2:length(list_in)) {
        a <- list_in[i]
        b <- list_in[i - 1]
        if (a == b + 1000) {
            list_out <- append(list_out, a)
        } else {
            list_out <- append(list_out, seq(b, a, 1000)[-1])
        }
    }
    return(list_out)
}

make_matrix <- function(chrom) {
    df_temp <- filter(input, X1 == chrom)
    mat_temp <- mat_blank
    for (i in seq_len(nrow(df_temp))) {
        mat_temp[df_temp$X4[i], df_temp$X7[i]] <- mat_temp[df_temp$X4[i], df_temp$X7[i]] + df_temp$X9[i]
    }
    return(mat_temp)
}


file <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/GSM2188866_MDA231_merged.ghits.pkbin.net.txt"
output_pdf <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/all-all.pdf"
output_csv <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/all-all.csv"

chrom_list <- paste0("chr", c(1:2))
col <- colorRampPalette(brewer.pal(9, "YlOrRd"))

cat("-------------- loading file -------------", "\n")
input <- readr::read_delim(file, delim = "\t", col_names = FALSE)
input <- data.frame(input)

cat("------------ creating matrix ------------", "\n")
# col names
for (i in 1:length(chrom_list)) {
    cat("-----------------",chrom_list[i], "------------------", "\r")
    position <- sort(unique((filter(input, X7 == chrom_list[i]))$X8))
    position <- fill_gap(position)
    position <- paste0(chrom_list[i], "-", position)
    if (i == 1) {
        mat_col <- position
    } else {
        mat_col <- append(mat_col, position)
    }
}
# row names
for (i in 1:length(chrom_list)) {
    position <- filter(input, X1 == chrom_list[i])
    position <- arrange(position, X2)
    position <- unique(position$X4)
    if (i == 1) {
        mat_row <- position
    } else {
        mat_row <- append(mat_row, position)
    }
}
# create matrix
mat_blank <- matrix(0, length(mat_row), length(mat_col))
colnames(mat_blank) <- mat_col
rownames(mat_blank) <- mat_row

cat("------------- arrange input -------------", "\n")
input <- unite(input, X7, c("X7", "X8"), sep = "-")

cat("----------- fullfiling matrix -----------", "\n")
cl <- makeForkCluster(8)
# clusterExport(cl, c("input", "mat_blank"))
mat_list <- parLapply(cl, chrom_list, make_matrix)
data_mat <- Reduce("+", mat_list)
stopCLuster(cl)

# cat("--------------- making csv ---------------", "\n")
# write.csv(data_mat, output_csv)

cat("------------- making heatmap -------------", "\n")
breaks <- c(0, 2^(0:max_legend(max(data_mat))))
pic <- pheatmap(data_mat, cluster_rows = FALSE, cluster_cols = FALSE,
			 	show_rownames = FALSE, show_colnames = FALSE,
				col = col(length(breaks)), breaks = breaks,
                legend = FALSE, border_color = NA,
                filename = output_pdf)
