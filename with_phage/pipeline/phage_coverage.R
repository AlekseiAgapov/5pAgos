library(ggplot2)

df <- read.table("coverage.tsv", header = T)
total_reads_number <- read.table("total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]
df$RPKM <- df$coverage/df$interval_length*1000/total_reads_number*1000000
df$interval_name <- df$interval_name * 1000

ggplot(df, aes(x = interval_name, y = RPKM)) +
  geom_area(fill = "grey66", col = 'black', size=1) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, kb") +
  ggtitle("sDNA coverage") +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000),
                     labels = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_y_continuous(limits = c(0, 10000),
                     breaks = c(0, 2500, 5000, 7500, 10000),
                     labels = c("0", "2.5", "5", "7.5", "10")) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'),
        title = element_text(size = 40),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "white",
                                          linetype = "longdash", 
                                          size = 0.3))

ggsave("both.png", width = 15, height = 10, dpi = 400)

df_plus <- read.table("plus_coverage.tsv", header = T)
df_minus <- read.table("minus_coverage.tsv", header = T)
df_comb <- data.frame(df$interval_name, df$interval_length)
df_comb$RPKM_plus <- df_plus$coverage/df_plus$interval_length*1000/total_reads_number*1000000
df_comb$RPKM_minus <- df_minus$coverage/df_minus$interval_length*1000/total_reads_number*1000000*-1

ggplot(df_comb, aes(x = df.interval_name)) +
  geom_area(aes(y = RPKM_plus), fill = "#005BBB") +
  geom_area(aes(y = RPKM_minus), fill = "#FFD500") +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, kb") +
  ggtitle("sDNA coverage") +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000),
                     labels = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_y_continuous(limits = c(-5000, 5000),
                     breaks = c(-5000, -2500, 0, 2500, 5000),
                     labels = c("-5", "-2.5", "0", "2.5", "5")) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'),
        title = element_text(size = 40),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "white",
                                          linetype = "longdash", 
                                          size = 0.3))

ggsave("plus_minus.png", width = 15, height = 10, dpi = 400)


chi_plus_df <- read.table("plus_chi.txt", header = F)
chi_plus_v <- chi_plus_df[,1]
chi_minus_df <- read.table("minus_chi.txt", header = F)
chi_minus_v <- chi_minus_df[,1]

ggplot(df_comb, aes(x = df.interval_name)) +
  geom_area(aes(y = RPKM_plus), fill = "#005BBB") +
  geom_area(aes(y = RPKM_minus), fill = "#FFD500") +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, kb") +
  ggtitle("sDNA coverage") +
  scale_x_continuous(breaks = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000),
                     labels = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_y_continuous(limits = c(-5500, 5000),
                     breaks = c(-5000, -2500, 0, 2500, 5000),
                     labels = c("-5", "-2.5", "0", "2.5", "5")) +
  annotate("segment", x = chi_plus_v, xend = chi_plus_v, y = -5500, yend = -4500, col = "#005BBB", size = 1) +
  annotate("segment", x = chi_minus_v, xend = chi_minus_v, y = -5500, yend = -4500, col = "#FFD500", size = 1) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'),
        title = element_text(size = 40),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "white",
                                          linetype = "longdash", 
                                          size = 0.3))
ggsave("plus_minus_chi.png", width = 15, height = 10, dpi = 400)
