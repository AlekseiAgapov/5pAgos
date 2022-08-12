library(ggplot2)

df <- read.table("coverage.tsv", header = T)
total_reads_number <- read.table("total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]
df$RPKM <- df$coverage/df$interval_length*1000/total_reads_number*1000000
df$interval_name = df$interval_name * 10000

df_plus <- read.table("plus_coverage.tsv", header = T)
df_minus <- read.table("minus_coverage.tsv", header = T)
df_comb <- data.frame(df$interval_name, df$interval_length)
df_comb$RPKM_plus <- df_plus$coverage/df_plus$interval_length*1000/total_reads_number*1000000
df_comb$RPKM_minus <- df_minus$coverage/df_minus$interval_length*1000/total_reads_number*1000000*-1

df_comb$ratio <- df_comb$RPKM_plus/df_comb$RPKM_minus*-1

av_mean <- c()
for(i in 1:length(df_comb$ratio)){
  x <- (df_comb$ratio[(i-2)] + df_comb$ratio[(i-1)] + df_comb$ratio[i] + df_comb$ratio[(i+1)] + df_comb$ratio[(i+2)])/5
  av_mean[i] <- x
}

av_mean[1] <- (df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[2] + df_comb$ratio[3])/5
av_mean[2] <- (df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[2] + df_comb$ratio[3] + df_comb$ratio[4])/5
av_mean[(length(av_mean) - 1)] <- (df_comb$ratio[(length(df_comb$ratio)-3)] + df_comb$ratio[(length(df_comb$ratio)-2)] + df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1])/5
av_mean[length(av_mean)] <- (df_comb$ratio[(length(df_comb$ratio)-2)] + df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[1])/5

df_comb$av_mean_ratio <- av_mean

ggplot(df_comb, aes(x = df.interval_name, y = av_mean_ratio)) +
  geom_line(col = "black", size = 2) +
  ylab("RPKM + strand / RPKM - strand") +
  xlab("chromosome coordinate, Mb") +
  scale_x_continuous(breaks = c(0, 1000000, 2000000, 3000000, 4000000),
                     labels = c(0, 1.0, 2.0, 3.0, 4.0)) +
  scale_y_continuous(limits = c(0.2, 5), breaks = c(0.25, 0.5, 1, 2, 4), trans = "log2") +
#  geom_vline(xintercept = c(73000), linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = 749000, linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = c(1329000, 1558000), linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = 3815000, linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
  geom_hline(yintercept = 1, colour = "black", size = 1) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("strand_ratio.png", width = 7, height = 7, dpi = 400)
