# binning sequence as certain length
TileSequence <- function(seqname, start, end, tilewidth){
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
                          count = rep(0, length(start_list)))
    return(GRanges_temp)
}

chr3L_mRNA_1K <- TileSequence("chr3L", 1, 28110227, 1000)
chrX_1k <- TileSequence("chrX", 1, 23648458, 1000)

DpnII <- read_delim("/Users/njduan/OneDrive/RNA-DNA/GSE97131/DpnII.txt",
                    delim = "\t", col_names = FALSE)
DpnII <- GRanges(seqnames = DpnII$X1,
                 ranges = IRanges(start = DpnII$X2, end = DpnII$X2),
                 strand = Rle("*"))
DpnII <- DpnII[DpnII@seqnames == "chrX"]

mrna <- read_delim("/Users/njduan/OneDrive/RNA-DNA/GSE97131/chr3L_mRNA_filted.txt",
                   delim = "\t", col_names = FALSE)
mrna <- GRanges(seqnames = Rle("chr3L"),
                 ranges = IRanges(start = mrna$X2, end = mrna$X2),
                 strand = Rle("*"))
roX2 <- read_delim("/Users/njduan/OneDrive/RNA-DNA/GSE97131/chrX_roX2.txt",
                   delim = "\t", col_names = FALSE)
roX2 <- GRanges(seqnames = roX2$X1,
                ranges = IRanges(start = roX2$X2, end = roX2$X2),
                strand = Rle("*"))
roX2 <- roX2[roX2@seqnames == "chrX"]

ten_m <- GRanges(seqnames = Rle("chr3L"),
                 ranges = IRanges(start = ten_m$X2, end = ten_m$X2),
                 strand = Rle("*"))

chr3L_mRNA_1k$count <- countOverlaps(chr3L_mRNA_1K, mrna)
chr3L_mRNA_1K$DPn_count <- countOverlaps(chr3L_mRNA_1K, DpnII)
chr3L_mRNA_1K$nor_count <- ifelse(chr3L_mRNA_1K$DPn_count == 0, 0, chr3L_mRNA_1K$count/chr3L_mRNA_1K$DPn_count)
chr3L_1k_ten_m$count <- countOverlaps(chr3L_1k_ten_m, ten_m)
chrX_1k$count <- countOverlaps(chrX_1k, roX2)
chrX_1k$Dpn_count <- countOverlaps(chrX_1k, DpnII)
chrX_1k$nor_count <- ifelse(chrX_1k$Dpn_count == 0, 0, chrX_1k$count/chrX_1k$Dpn_count)

GA_track <- GenomeAxisTrack()
I_track <- IdeogramTrack(chromosome = "chrX", genome = "dm6")
chrX_1k_track <- DataTrack(chrX_1k, type = "h")
chr3L_1k_ten_m_track <- DataTrack(chr3L_1k_ten_m, type = "h")
plotTracks(c(I_track, GA_track, chrX_1k_track),
           from = 1, to = 23648458)
