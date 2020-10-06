resize <- function(input, chrom_list, size) {
    for (i in seq_len(length(chrom_list))) {
        temp1 <- filter(input, X7 == chrom_list[i])
        rna_list <- unique(input$X4)
        for (j in seq_len(length(rna_list))) {
            temp2 <- filter(temp1, X4 == rna_list[j])
            local_list <- vector()
            for (k in seq_len(nrow(temp2))) {
                if (temp2$X8[k] %% size == 0) {
                    local_list <- append(local_list, temp2$X8[k] %/% size)
                } else {
                    local_list <- append(local_list, temp2$X8[k] %/% size + 1)
                }
            }
            temp2$X8 <- local_list
            temp3 <- temp2 %>% group_by(X1, X2, X3, X4, X5, X6, X7, X8) %>% summarize(X9 = sum(X9))
            temp3 <- data.frame(temp3)
            if (j == 1) {
                output1 <- temp3
            } else {
                output1 <- rbind(output1, temp3)
            }
        }
        if (i == 1) {
            output2 <- output1
        } else {
            output2 <- rbind(output2, output1)
        }
    }
    return(output2)
}

hgenome <- data.frame (chrom = paste0("chr", c(1:22, "X", "Y")), 
                       size = c(248965422, 242193529, 198295559, 190214555, 181538259, 
                                170805979, 159345973, 145138636, 138394717, 133797422, 
                                135086622, 133275309, 114364328, 107043718, 101991189,
                                90338345, 83257441, 80373285, 58617616, 64444167,
                                46709983, 50818468, 156040895, 57227415))
for (i in seq_len(length(chrom_list))) {
    cat("-----------------", chrom_list[i], "------------------", "\r")
    chrom_size <- (filter(hgenome, chrom == chrom_list[i]))$size
    position <- seq(1, chrom_size %/% size + 1, 1)
    position <- paste0(chrom_list[i], "-", position)
    if (i == 1) {
        mat_col <- position
    } else {
        mat_col <- append(mat_col, position)
    }
}


for (i in seq_len(length(chrom_list))) {
    cat("-----------------", chrom_list[i], "------------------", "\r")
    position <- sort(unique((filter(re_input, X7 == chrom_list[i]))$X8))
    position <- fill_gap(position)
    position <- paste0(chrom_list[i], "-", position)
    if (i == 1) {
        mat_col <- position
    } else {
        mat_col <- append(mat_col, position)
    }
}