library(tidyr)
library(pheatmap)

csv_input <- "/Volumes/Data/RNA-DNA/GSE79230/THSS_chr13.csv"

input <- readr::read_csv(csv_input)
input <- as.data.frame(input)
input <- spread(input, position, count)

mat <- as.matrix(input[,2:ncol(input)])
rownames(mat) <- input$gene
mat <- mat[rev(order(rowSums(mat))),]

pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, show_rownames = FALSE)


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

library(RColorBrewer)
col <- colorRampPalette(brewer.pal(9, "YlOrRd"))
breaks <- c(0, 2^(0:max_legend(max(mat))))
pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
                show_rownames = FALSE, show_colnames = FALSE,
                col = col(length(breaks)), breaks = breaks,
                legend = FALSE, border_color = NA, cellheight = 2, cellwidth = 2)
