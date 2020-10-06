library(tidyr)
library(dplyr)


count_distribute <- function(input, anno) {
    gene_list <- unique(input$X4)
    genelist_length <- length(gene_list)
    output <- data.frame()
    n <- 1
    for (i in seq_len(genelist_length)) {
        if (i >= (genelist_length * n) %/% 100) {
            cat("----------------- ", n, "% ------------------", "\r", sep = "")
            n <- n + 1
        }

        if (gene_list[i] %in% anno$id == FALSE) {
            next
        }

        # anno data
        gene_anno <- filter(anno, id == gene_list[i])
        chrom_anno <- gene_anno$chr
        start_anno <- gene_anno$start - 10000
        end_anno <- gene_anno$end + 10000

        # input
        gene_input <- filter(input, X4 == gene_list[i])
        gene_output <- data.frame(id = gene_anno$id,
                                  name = gene_anno$name,
                                  local = 0,
                                  cis = 0,
                                  trans = 0,
                                  stringsAsFactors = FALSE)
        for (j in seq_len(nrow(gene_input))) {
            chrom_gene <- gene_input$X7[j]
            local_gene <- gene_input$X8[j]
            if (chrom_gene != chrom_anno) {
                gene_output$trans <- gene_output$trans + gene_input$X9[j] # trans
            } else {
                if (local_gene >= start_anno & local_gene <= end_anno) {
                    gene_output$local <- gene_output$local + gene_input$X9[j] # local
                } else {
                    gene_output$cis <- gene_output$cis + gene_input$X9[j] # cis
                }
            }
        }
        if (i == 1) {
            output <- gene_output
        } else {
            output <- rbind(output, gene_output)
        }
    }
    return(output)
}

file_path <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/GSM2188866_MDA231_merged.ghits.pkbin.net.txt"
anno_path <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/hs_gene.csv"
output_csv <- "/data/agl_data/ningjun/RNA-DNA/GSE82312/distribute.csv"

cat("-------------- loading file -------------", "\n")
input <- readr::read_delim(file_path, delim = "\t", col_names = FALSE)
input <- data.frame(input)
anno <- readr::read_delim(anno_path, delim = "\t", col_names = TRUE)
anno <- data.frame(anno)

cat("-------------- calculating -------------", "\n")
output <- count_distribute(input, anno)

output$total <- output$local + output$cis + output$trans
output$local_per <- output$local / output$total
output$cis_per <- output$cis / output$total
output$trans_per <- 1 - output$local_per - output$cis_per

cat("-------------- making csv --------------", "\n")
write.csv(output, output_csv, row.names = FALSE)

cat("---------------- Done! -----------------", "\n")