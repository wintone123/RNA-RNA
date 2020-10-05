rm(list = ls())

data_500k <- cooler2bedpe("/Volumes/RaEVO/RNA-RNA/GSE127188/temp_500000.cool")
data_cis <- data_500k$cis
data_trans <- data_500k$trans

# cis interaction
data_1v1 <- data_cis$`1`
df_1v1 <- data.frame(chr1 = data_1v1$chr1,
                     start1 = data_1v1$start1 / 500000,
                     chr2 = data_1v1$chr2,
                     start2 = data_1v1$start2 / 500000,
                     value = log2(data_1v1$IF))
df_1v1$start1 <- str_pad(df_1v1$start1, width = 3, pad = 0)
df_1v1$start2 <- str_pad(df_1v1$start2, width = 3, pad = 0)
df_1v1 <- unite(df_1v1, "from", c("chr1", "start1"))
df_1v1 <- unite(df_1v1, "to", c("chr2", "start2"))

ggplot(df_1v1, aes(to, from)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_gradient(low = "#ffffff", high = "#ff0000") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

# autosome interaction
data_full <- data_trans

for (i in 1:length(data_cis)){
    data_full <- rbind(data_full, data_cis[[i]])
}

data_full <- filter(data_full, (chr1 %in% 1:22) & (chr2 %in% 1:22) )

df_full <- data.frame(chr1 = data_full$chr1,
                      start1 = data_full$start1 / 500000,
                      chr2 = data_full$chr2,
                      start2 = data_full$start2 / 500000,
                      value = log2(data_full$IF))
df_full$chr1 <- str_pad(df_full$chr1, width = 2, pad = 0)
df_full$start1 <- str_pad(df_full$start1, width = 3, pad = 0)
df_full$chr2 <- str_pad(df_full$chr2, width = 2, pad = 0)
df_full$start2 <- str_pad(df_full$start2, width = 3, pad = 0)
df_full <- unite(df_full, "from", c("chr1","start1"))
df_full <- unite(df_full, "to", c("chr2", "start2"))

ggplot(df_full, aes(to, from)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_gradient(low = "#ffffff", high = "#ff0000") +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
