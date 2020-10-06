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

for (i in seq_len(length(chr4_100k))) {
    start <- start(chr4_100k[i])
    end <- end(chr4_100k[i])
    count <- 0
    filted <- filter(input, X8 >= start & X8 <= end)
    for (j in seq_len(nrow(filted))) {
        count <- count + filted$X9[j]
    }
    chr4_100k[i]$count <- count
}

anno <- GRanges(seqnames = Rle("chr4"),
                ranges = IRanges(start = gene_pc$start, end = gene_pc$end),
                strand = Rle("*"),
                name = gene_pc$name)

GA_track <- GenomeAxisTrack()
I_track <- IdeogramTrack(chromosome = "chr4", genome = "hg38")
GR_track <- GeneRegionTrack(anno)
data_track <- DataTrack(chr4_100k, type = "h")
plotTracks(c(I_track, GA_track, data_track, GR_track), 
           from = 100000000, to = 175000000)
