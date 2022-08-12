library(ggplot2)

df <- read.table("coverage.tsv", header = T)
total_reads_number <- read.table("total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]
df$RPKM <- df$coverage/df$interval_length*1000/total_reads_number*1000000
df$interval_name <- df$interval_name * 1000

ggplot(df, aes(x = interval_name, y = RPKM)) +
  geom_col(fill = "black", width = 1000) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, Mb") +
  scale_x_continuous(breaks = c(0, 1000000, 2000000, 3000000, 4000000),
                     labels = c(0, 1.0, 2.0, 3.0, 4.0)) +
  scale_y_continuous(limits = c(-100, 4500),
                     breaks = c(0, 1000, 2000, 3000, 4000, 5000),
                     labels = c("0.0", "1", "2", "3", "4", "5")) +
#  geom_vline(xintercept = c(73000), linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = 749000, linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = c(1329000, 1558000), linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = 3815000, linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
  geom_hline(yintercept = 0, colour="black", size=1) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 50, colour = 'black'))

ggsave("both.png", width = 10, height = 7, dpi = 400)

df_plus <- read.table("plus_coverage.tsv", header = T)
df_minus <- read.table("minus_coverage.tsv", header = T)
df_comb <- data.frame(df$interval_name, df$interval_length)
df_comb$RPKM_plus <- df_plus$coverage/df_plus$interval_length*1000/total_reads_number*1000000
df_comb$RPKM_minus <- df_minus$coverage/df_minus$interval_length*1000/total_reads_number*1000000*-1

ggplot(df_comb, aes(x = df.interval_name)) +
  geom_col(aes(y = RPKM_plus), fill = "#005BBB", width = 1000) +
  geom_col(aes(y = RPKM_minus), fill = "#CCAA00", width = 1000) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, Mb") +
#  ggtitle("sDNA coverage") +
  scale_x_continuous(breaks = c(0, 1000000, 2000000, 3000000, 4000000),
                     labels = c(0, 1.0, 2.0, 3.0, 4.0)) +
  scale_y_continuous(limits = c(-2400, 3500),
                     breaks = c(-1500, 0, 1500, 3000),
                     labels = c("-1.5", "0", "1.5", "3.0")) +
#  geom_vline(xintercept = c(73000), linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = 749000, linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = c(1329000, 1558000), linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
#  geom_vline(xintercept = 3815000, linetype = "dashed", colour = "grey30", size = 1, alpha=0.5) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("plus_minus.png", width = 7, height = 7, dpi = 400)


chi_plus_df <- read.table("plus_chi.txt", header = F)
chi_plus_v <- chi_plus_df[,1]
length(chi_plus_v)
nrow(subset(chi_plus_df, V1 < 1200000 | V1 > 1700000)) 
chi_minus_df <- read.table("minus_chi.txt", header = F)
chi_minus_v <- chi_minus_df[,1]
length(chi_minus_v)
nrow(subset(chi_minus_df, V1 < 1200000 | V1 > 1700000)) 

ggplot(df_comb, aes(x = df.interval_name)) +
  geom_col(aes(y = RPKM_plus), fill = "#005BBB", width = 1000) +
  geom_col(aes(y = RPKM_minus), fill = "#CCAA00", width = 1000) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, Mb") +
  scale_x_continuous(limits = c(1200000, 1700000),
                     breaks = c(1300000, 1600000),
                     labels = c(1.3, 1.6)) +
  scale_y_continuous(limits = c(-2400, 3500),
                     breaks = c(-1500, 0, 1500, 3000),
                     labels = c("-1.5", "0", "1.5", "3.0")) +
#  geom_vline(xintercept = c(1267000, 1328000, 1557000, 1629000), linetype = "dashed", colour = "grey30", size = 1) +
#  annotate("segment", x = chi_plus_v, xend = chi_plus_v, y = -2500, yend = -2100, col = "#005BBB", size = 1) +
#  annotate("segment", x = chi_minus_v, xend = chi_minus_v, y = -2500, yend = -2100, col = "#CCAA00", size = 1) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("Ter_region.png", width = 7, height = 7, dpi = 400)


ggplot(df_comb, aes(x = df.interval_name)) +
  geom_col(aes(y = RPKM_plus), fill = "#005BBB", width = 1000) +
  geom_col(aes(y = RPKM_minus), fill = "#CCAA00", width = 1000) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, Mb") +
  #  ggtitle("sDNA coverage") +
  scale_x_continuous(limits = c(0, 210000),
                     breaks = c(0, 100000, 200000),
                     labels = c(0, 0.1, 0.2)) +
  scale_y_continuous(limits = c(-450, 400),
                     breaks = c(-400, -200, 0, 200, 400),
                     labels = c(-0.4, -0.2, 0, 0.2, 0.4)) +
#  geom_vline(xintercept = c(73000), linetype = "dashed", colour = "grey30", size = 1) +
#  annotate("segment", x = chi_plus_v, xend = chi_plus_v, y = -450, yend = -350, col = "#005BBB", size = 1) +
#  annotate("segment", x = chi_minus_v, xend = chi_minus_v, y = -450, yend = -350, col = "#CCAA00", size = 1) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("araC_region.png", width = 7, height = 7, dpi = 400)

ggplot(df_comb, aes(x = df.interval_name)) +
  geom_col(aes(y = RPKM_plus), fill = "#005BBB", width = 1000) +
  geom_col(aes(y = RPKM_minus), fill = "#CCAA00", width = 1000) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, Mb") +
  #  ggtitle("sDNA coverage") +
  scale_x_continuous(limits = c(650000, 850000),
                     breaks = c(700000, 800000),
                     labels = c(0.7, 0.8)) +
  scale_y_continuous(limits = c(-250, 260),
                     breaks = c(-200, -100, 0, 100, 200),
                     labels = c(-0.2, -0.1, 0, 0.1, 0.2)) +
#  geom_vline(xintercept = c(748000, 791000), linetype = "dashed", colour = "grey30", size = 1) +
#  annotate("segment", x = chi_plus_v, xend = chi_plus_v, y = -250, yend = -200, col = "#005BBB", size = 1) +
#  annotate("segment", x = chi_minus_v, xend = chi_minus_v, y = -250, yend = -200, col = "#CCAA00", size = 1) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("DE3_region.png", width = 7, height = 7, dpi = 400)
