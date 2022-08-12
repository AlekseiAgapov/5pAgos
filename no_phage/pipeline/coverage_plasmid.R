library(ggplot2)

df <- read.table("coverage.tsv", header = T)
total_reads_number <- read.table("total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]
df$plus_RPKM <- df$plus_coverage/(df$end - df$start)*1000/total_reads_number*1000000
df$minus_RPKM <- df$minus_coverage/(df$end - df$start)*1000/total_reads_number*1000000*-1


ggplot(df, aes(x = coordinate)) +
  geom_col(aes(y = plus_RPKM), fill = "#005BBB", width = 10) +
  geom_col(aes(y = minus_RPKM), fill = "#CCAA00", width = 10) +
  ylab("RPKM × 1000") +
  xlab("plasmid coordinate, kb") +
  scale_x_continuous(breaks = c(0, 2000, 4000, 6000),
                     labels = c(0, 2, 4, 6)) +
  scale_y_continuous(limits = c(-4200, 4200),
                     breaks = c(-4000, -2000, 0, 2000, 4000),
                     labels = c("-4", "-2", "0", "2", "4")) +
#  geom_vline(xintercept = c(1200, 1800), linetype = "dashed", colour = "grey30", size = 1) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("plasmid.png", width = 7, height = 7, dpi = 400)
