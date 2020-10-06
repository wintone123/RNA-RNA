library(tidyr)
library(dplyr)
library(parallel)
library(pheatmap)
library(RColorBrewer)

peak_map <- function(chrom) { # chromosome of target position
    anno <- filter(anno, chr == chrom)
    input_1 <- filter(input, X1 == chrom & X7 == chrom)
    rna_id_list <- unique(input_1$X4)
    rna_list_length <- length(rna_id_list)
    out_mat <- matrix(0, rna_list_length, 2 * bin_num + 1)
    rownames(out_mat) <- rna_id_list
    colnames(out_mat) <- seq(2 * bin_num + 1)

    for (i in seq_len(rna_list_length)) {

        rna_id <- rna_id_list[i]
        if (rna_id %in% anno$id) {
            gene <- filter(anno, id == rna_id)
            start <- gene$start[1]
            up <- start - range
            end <- gene$end[1]
            down <- end + range
        } else {
            next
        }

        rna_list <- filter(input_1, X4 == rna_id)

        for (j in seq_len(nrow(rna_list))) {
            position <- rna_list$X8[j]
            # print(c(rna_id, position))
            if (position >= up & position <= down) {
                if (position >= start & position <= end) {
                    out_mat[rna_id, center] <- out_mat[rna_id, center] + input_1$X9[j]
                } else if (position < start) {
                    position <- start - position
                    if (position %% binsize == 0) {
                        out_mat[rna_id, center - position %/% binsize] <- out_mat[rna_id, center - position %/% binsize] + input_1$X9[j]
                    } else {
                        out_mat[rna_id, center - position %/% binsize - 1] <- out_mat[rna_id, center - position %/% binsize - 1] + input_1$X9[j]
                    }
                } else if (position > end) {
                    position <- position - end
                    if (position %% binsize == 0) {
                        out_mat[rna_id, center + position %/% binsize] <- out_mat[rna_id, center + position %/% binsize] + input_1$X9[j]
                    } else {
                        out_mat[rna_id, center + position %/% binsize + 1] <- out_mat[rna_id, center + position %/% binsize + 1] + input_1$X9[j]
                    }
                }
            }
        }
    }

    return(out_mat)    
}

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

file_path <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/GSM2188866_MDA231_merged.ghits.pkbin.net.txt"
anno_path <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/hs_gene.csv"
output_csv <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/heatmap.csv"
output_tiff <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/heatmap.tiff"

range <- 2000000
binsize <- 20000
bin_num <- range %/% binsize
center <- bin_num + 1
chrom_list <- paste0("chr", c(1:22, "X"))


cat("-------------- loading file -------------", "\n")
input <- readr::read_delim(file_path, delim = "\t", col_names = FALSE)
input <- data.frame(input)
anno <- readr::read_delim(anno_path, delim = ",", col_names = TRUE)
anno <- data.frame(anno)

cat("---------- fullfilling matrix -----------", "\n")
# anno <- filter(anno, chr == chrom)
# input_1 <- filter(input, X1 == chrom & X7 == chrom)
# rna_id_list <- unique(input_1$X4)
# rna_list_length <- length(rna_id_list)
# out_mat <- matrix(0, rna_list_length, 2 * bin_num + 1)
# rownames(out_mat) <- rna_id_list
# colnames(out_mat) <- seq(2 * bin_num + 1)

# n <- 1
# for (i in seq_len(rna_list_length)) {
#     if (i >= (rna_list_length * n) %/% 100) {
#         cat("----------------- ", n, "% ------------------", "\r", sep = "")
#         n <- n + 1
#     }

#     rna_id <- rna_id_list[i]
#     if (rna_id %in% anno$id) {
#         gene <- filter(anno, id == rna_id)
#         start <- gene$start[1]
#         up <- start - range
#         end <- gene$end[1]
#         down <- end + range
#     } else {
#         next
#     }

#     rna_list <- filter(input_1, X4 == rna_id)

#     for (j in seq_len(nrow(rna_list))) {
#         position <- rna_list$X8[j]
#         # print(c(rna_id, position))
#         if (position >= up & position <= down) {
#             if (position >= start & position <= end) {
#                 out_mat[rna_id, center] <- out_mat[rna_id, center] + input_1$X9[j]
#             } else if (position < start) {
#                 position <- start - position
#                 if (position %% binsize == 0) {
#                     out_mat[rna_id, center - position %/% binsize] <- out_mat[rna_id, center - position %/% binsize] + input_1$X9[j]
#                 } else {
#                     out_mat[rna_id, center - position %/% binsize - 1] <- out_mat[rna_id, center - position %/% binsize - 1] + input_1$X9[j]
#                 }
#             } else if (position > end) {
#                 position <- position - end
#                 if (position %% binsize == 0) {
#                     out_mat[rna_id, center + position %/% binsize] <- out_mat[rna_id, center + position %/% binsize] + input_1$X9[j]
#                 } else {
#                     out_mat[rna_id, center + position %/% binsize + 1] <- out_mat[rna_id, center + position %/% binsize + 1] + input_1$X9[j]
#                 }
#             }
#         }
#     }
# }
# out_mat <- out_mat[order(rowSums(out_mat), decreasing = TRUE),]

core <- detectCores() %/% 3
cl <- makeForkCluster(core)
mat_list <- parLapply(cl, chrom_list, peak_map)
data_mat <- Reduce("rbind", mat_list)
stopCluster(cl)
data_mat <- data_mat[order(rowSums(data_mat), decreasing = TRUE),]

cat("--------------- making csv ---------------", "\n")
# write.csv(out_mat, output_csv)
write.csv(data_mat, output_csv)

cat("------------- making heatmap ------------", "\n")
col <- colorRampPalette(brewer.pal(9, "YlOrRd"))
breaks <- c(0, 2^(0:max_legend(max(data_mat))))
pic <- pheatmap(data_mat, cluster_rows = FALSE, cluster_cols = FALSE,
                show_rownames = FALSE, show_colnames = FALSE,
				col = col(length(breaks)), breaks = breaks,
                legend = FALSE, border_color = NA,
                filename = output_tiff)

cat("----------------- Done! -----------------", "\n")