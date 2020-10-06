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

sum_overlap <- function(query, subject, overlap_list) {
    overlap_list <- as.data.frame(overlap_list)
    overlap_list$count <- subject[overlap_list$subjectHits]$score
    overlap_list <- group_by(overlap_list, queryHits) %>% summarise(count = sum(count))
    query[overlap_list$queryHits]$count <-  overlap_list$count
    return(query)
}

roX2_ChIRP <- import.bedGraph("/Volumes/Data/RNA-DNA/GSE97456/GSE69208_MEL-roX2.merge.bedGraph")

roX2_ChIRP_X <- roX2_ChIRP[roX2_ChIRP@seqnames == "chrX"]

chrX <- TileSequence("chrX", 1, 23542271, 200)
overlap_list <- findOverlaps(chrX, roX2_ChIRP_X)
chrX <- sum_overlap(chrX, roX2_ChIRP_X, overlap_list)
chrX$count <- ifelse(chrX$count / 100 <= 1, 0, log2(chrX$count / 100))

ga_track <- GenomeAxisTrack()
ig_track <- IdeogramTrack(chromosome = "chrX", 'dm6')
a_track <- DataTrack(chrX, type = "h", )
plotTracks(c(ig_track, ga_track, a_track), from = 10000000, to = 23540000)
