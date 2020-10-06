library(dplyr)
library(tidyr)
library(bsgenome)
library(circlize)

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

chrom_size_list <- c(23513712, 25286936, 28111227, 32079331, 1348131, 23648458)
chrom_list <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
bin_size <- 20000

input <- readr::read_delim("/Users/njduan/OneDrive/RNA-DNA/GSE97131/roX2.txt",
                           delim = "\t", col_names = FALSE)
input_g <- GRanges(seqnames = input$X1,
                   ranges = IRanges(start = input$X2, end = input$X2),
                   strand = Rle("*"))
DpnII <- readr::read_delim("/Users/njduan/OneDrive/RNA-DNA/GSE97131/DpnII.txt",
                                 delim = "\t", col_names = FALSE)
DpnII_g <- GRanges(seqnames = DpnII$X1,
                         ranges = IRanges(start = DpnII$X2, end = DpnII$X2),
                         strand = Rle("*"))

for (i in seq_len(length(chrom_list))) {
    chrom <- chrom_list[i]
    end <- chrom_size_list[i]
    bin_temp <- TileSequence(chrom, 1, end, bin_size)
    if (i == 1) {
      bin <- bin_temp
    } else {
      bin <- append(bin, bin_temp)
    }
}

bin$count <- countOverlaps(bin, input_g)
bin$DpnII_count <-countOverlaps(bin, DpnII_g)

bin_df <- data.frame(chrom = bin@seqnames,
                     start = start(bin),
                     end = end(bin),
                     count = bin$count / bin$DpnII_count)

bin_filte <- bin[bin$count != 0]
links <- data.frame(chrom1 = bin_filte@seqnames,
                    start1 = start(bin_filte),
                    end1 = end(bin_filte))
links$chrom2 <- rep("chrX", nrow(links))
links$start2 <- rep(11474129, nrow(links))
links$end2 <- rep(11474129, nrow(links))


circos.clear()
circos.genomicInitialize(bin_df)
circos.genomicTrackPlotRegion(bin_df, ylim = c(0,1), track.height = 0.1, bg.border = '#ffffffff',
                              panel.fun = function(region, value, ...) {
                                for (i in seq_len(nrow(region))) {
                                  reg <- region[i, ]
                                  val <- value[i, ] / 28.67
                                  circos.genomicRect(reg, val, ybottom = 0, ytop = val, border = '#00000090')
                                }
                              })
circos.track(ylim = c(0,1), 
             bg.co = c("#0000FFFF", "#0000FFFF", "#00FF00FF","#00FF00FF", "#FF00FFFF", "#FF0000FF"),
             bg.border = NA, track.height = 0.05)
circos.genomicLink(links[1:3], links[4:6], col = "#00000001", lwd = 0.7)
