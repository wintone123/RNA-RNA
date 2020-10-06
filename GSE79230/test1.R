# load libraries
library(BSgenome)
library(tidyr)
library(dplyr)

# load functions
TileSequence_po <- function(seqname, gene, start, end, tilewidth){
    start_list <- seq(start, end, by = tilewidth)
    end_list <- start_list + tilewidth -1
    if (start_list[length(start_list)] == end) {
        end_list <- end_list[c(1:length(end_list)-1)]
        end_list[length(end_list)] <- start_list[length(start_list)]
        start_list <- start_list[c(1:length(start_list)-1)]
    } else {
        end_list[length(end_list)] <- end
    }
    GRanges_temp <- GRanges(seqnames = Rle(seqname),
                            ranges = IRanges(start = start_list, end = end_list),
                            strand = Rle("*"),
                            gene = rep(gene, length(start_list)),
                            position = seq(length(start_list)),
                            count = rep(0, length(start_list)))
    return(GRanges_temp)
}

def_pro_po <- function(anno_df) {
    df_len <- nrow(anno_df)
    n <- 1
    for (i in seq_len(df_len)) {
        if (i >= (df_len %/% 100) * n) {
            cat("-----------------", n, "%----------------", "\r", sep = "")
            n <- n + 1
        }
        seqname <- paste0("chr", anno_df$chr[i])
        gene <- anno_df$gene_id[i]
        start <- anno_df$start[i]
        tss_up <- start - promoter_range
        tss_down <- start + promoter_range
        if (i == 1) {
            output <- TileSequence_po(seqname, gene, tss_up, tss_down, bin_size)
        } else {
            output <- append(output, TileSequence_po(seqname, gene, tss_up, tss_down, bin_size))
        }
    }
    return(output)
}

TileSequence_ne <- function(seqname, gene, start, end, tilewidth){
    start_list <- seq(start, end, by = tilewidth)
    end_list <- start_list + tilewidth -1
    if (start_list[length(start_list)] == end) {
        end_list <- end_list[c(1:length(end_list)-1)]
        end_list[length(end_list)] <- start_list[length(start_list)]
        start_list <- start_list[c(1:length(start_list)-1)]
    } else {
        end_list[length(end_list)] <- end
    }
    GRanges_temp <- GRanges(seqnames = Rle(seqname),
                            ranges = IRanges(start = start_list, end = end_list),
                            strand = Rle("*"),
                            gene = rep(gene, length(start_list)),
                            position = rev(seq(length(start_list))),
                            count = rep(0, length(start_list)))
    return(GRanges_temp)
}

def_pro_ne <- function(anno_df) {
    df_len <- nrow(anno_df)
    n <- 1
    for (i in seq_len(df_len)) {
        if (i >= (df_len %/% 100) * n) {
            cat("-----------------", n, "%----------------", "\r", sep = "")
            n <- n + 1
        }
        seqname <- paste0("chr", anno_df$chr[i])
        gene <- anno_df$gene_id[i]
        start <- anno_df$start[i]
        tss_up <- start - promoter_range
        tss_down <- start + promoter_range
        if (i == 1) {
            output <- TileSequence_ne(seqname, gene, tss_up, tss_down, bin_size)
        } else {
            output <- append(output, TileSequence_ne(seqname, gene, tss_up, tss_down, bin_size))
        }
    }
    return(output)
}

sum_overlap <- function(query, subject, overlap_list) {
    unique_overlap <- unique(queryHits(overlap_list))
    list_len <- length(unique_overlap)
    n <- 1
    for (i in seq_len(list_len)) {
        if (i >= (list_len %/% 100) * n) {
            cat("-----------------", n, "%----------------", "\r", sep = "")
            n <- n + 1
        }
        q <- unique_overlap[i]
        s_list <- subjectHits(overlap_list[queryHits(overlap_list) == q])
        query[q]$count <- query[q]$count + sum(subject[s_list]$score)
    }
    return(query)
}

# parameter settings
chrom_list <- paste0("chr", c(1:19, "X", "Y"))
promoter_range <- 1500
bin_size <- 100

# input and output
anno_path <- "/data/agl_data/ningjun/anno/anno_mouse/"
wig_input <- "/data/agl_data/ningjun/RNA-DNA/GSE79230/SPERM_ATACseq_THSS.wig"
csv_output <- "/data/agl_data/ningjun/RNA-DNA/GSE79230/THSS.csv"

# load annotations
cat("==========load annotations==========", "\n")
for (n in seq_len(length(chrom_list))) {
    if (n == 1) {
        anno_mm <- readr::read_csv(paste0(anno_path, "anno_", chrom_list[n], ".csv"))
    } else {
        anno_mm <- rbind(anno_mm, readr::read_csv(paste0(anno_path, "anno_", chrom_list[n], ".csv")))
    }
}

# promoter definition
cat("========promoter definition========", "\n")
# positive strand
anno_po <- filter(anno_mm, strand == "+")
anno_po <- def_pro_po(anno_po)
# negative strand
anno_ne <- filter(anno_mm, strand == "-")
anno_ne <- def_pro_ne(anno_ne)

# load THSS data
cat("===========load THSS data===========", "\n")
thss_mm <- import.wig(wig_input) 

# data procession
cat("============data process============", "\n")
# count overlap
po_list <- findOverlaps(anno_po, thss_mm)
anno_po <- sum_overlap(anno_po, thss_mm, po_list)
ne_list <- findOverlaps(anno_ne, thss_mm)
anno_ne <- sum_overlap(anno_ne, thss_mm, ne_list)
# extract data from Grange
output_po <- data.frame(gene = anno_po$gene,
                        position = anno_po$position,
                        count = anno_po$count,
                        stringsAsFactors = FALSE)
output_ne <- data.frame(gene = anno_ne$gene,
                        position = anno_ne$position,
                        count = anno_ne$count,
                        stringsAsFactors = FALSE)
output_merge <- rbind(output_po, output_ne)
# output_merge <- spread(output_merge, position, count)

# output
# csv
cat("============writing csv============", "\n")
write.csv(output_merge, csv_output, col.names = TRUE, sep = "\t")

