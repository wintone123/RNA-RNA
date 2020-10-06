library(BSgenome)
library(tidyr)
library(dplyr)
library(Gviz)

# function
sum_overlap <- function(query, subject, overlap_list) {
    query$count <- 0
    overlap_list <- as.data.frame(overlap_list)
    overlap_list$count <- subject[overlap_list$subjectHits]$score
    overlap_list <- group_by(overlap_list, queryHits) %>% summarise(count = sum(count))
    query[overlap_list$queryHits]$count <-  overlap_list$count
    return(query)
}

# import file path
anno_path <- "/Volumes/Data/RNA-DNA/Drosophila_melanogaster.BDGP6.28.47.chr.gff3"
roX2_chrX_path <- "/Volumes/Data/RNA-DNA/GSE97456/roX2_ChIRP.bedGraph"

# anno load
anno <- import.gff3("/Volumes/Data/RNA-DNA/Drosophila_melanogaster.BDGP6.28.47.chr.gff3",
                    sequenceRegionsAsSeqinfo = FALSE)

# anno process
# chromosome X
anno_chrX <- anno[anno@seqnames == "X"]
anno_chrX@seqnames <- Rle(rep("chrX", length(anno_chrX)))
# promoter area (TSS +- 2kb)
anno_chrX_promoter <- anno_chrX[anno_chrX$type == "gene"]
anno_chrX_promoter <- anno_chrX_promoter[,0]
start(anno_chrX_promoter) <- start(anno_chrX_promoter) - 2000
end(anno_chrX_promoter) <- start(anno_chrX_promoter) + 4000
anno_chrX_promoter <- GRanges(seqnames = anno_chrX_promoter@seqnames,
                              ranges = anno_chrX_promoter@ranges,
                              strand = anno_chrX_promoter@strand)
# exon area
anno_chrX_exon <- anno_chrX[anno_chrX$type == "exon"]
anno_chrX_exon <- anno_chrX_exon[,0]
anno_chrX_exon <- unique(anno_chrX_exon)
anno_chrX_exon <- GRanges(seqnames = anno_chrX_exon@seqnames,
                              ranges = anno_chrX_exon@ranges,
                              strand = anno_chrX_exon@strand)

# read load
roX2_chrX <- import.bedGraph(roX2_chrX_path)
roX2_chrX@seqnames <- Rle(rep("X", length(roX2_chrX)))

# alignment
roX2_chrX_promoter <- sum_overlap(anno_chrX_promoter, roX2_chrX, findOverlaps(anno_chrX_promoter, roX2_chrX))
roX2_chrX_promoter$count <- ifelse(roX2_chrX_promoter@strand == "+", 
                                   roX2_chrX_promoter$count, -roX2_chrX_promoter$count)
roX2_chrX_promoter@strand <- Rle(rep("*", length(roX2_chrX_promoter)))
roX2_chrX_exon <- sum_overlap(anno_chrX_exon, roX2_chrX, findOverlaps(anno_chrX_exon, roX2_chrX))
roX2_chrX_exon$count <- ifelse(roX2_chrX_exon@strand == "+", 
                                   roX2_chrX_exon$count, -roX2_chrX_exon$count)
roX2_chrX_exon@strand <- Rle(rep("*", length(roX2_chrX_exon)))

# image
GA_track <- GenomeAxisTrack()
Id_track <- IdeogramTrack(chromosome = "chrX", genome = "dm6")
promoter_track <- DataTrack(roX2_chrX_promoter, type = "h", name = "promoter")
exon_track <- DataTrack(roX2_chrX_exon, type = "h", name = "exon")
plotTracks(c(Id_track, GA_track, promoter_track, exon_track))
