library(ggplot2)

total_reads_number <- read.table("total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]

df_plus_plus <- read.table("plus_plus.tsv", header = T)
df_plus_plus$RPKM <- df_plus_plus$coverage/df_plus_plus$interval_length*1000/total_reads_number*1000000

df_minus_plus <- read.table("minus_plus.tsv", header = T)
df_minus_plus$RPKM <- df_minus_plus$coverage/df_minus_plus$interval_length*1000/total_reads_number*1000000

df_minus_minus <- read.table("minus_minus.tsv", header = T)
df_minus_minus$RPKM <- df_minus_minus$coverage/df_minus_minus$interval_length*1000/total_reads_number*1000000

df_plus_minus <- read.table("plus_minus.tsv", header = T)
df_plus_minus$RPKM <- df_plus_minus$coverage/df_plus_minus$interval_length*1000/total_reads_number*1000000

df_sum_right <- rbind(df_plus_plus, df_minus_minus)
df_sum_wrong <- rbind(df_plus_minus, df_minus_plus)

agr_plus_plus <- aggregate(x = df_plus_plus, by = list(df_plus_plus$interval_name), FUN = mean)
agr_minus_plus <- aggregate(x = df_minus_plus, by = list(df_minus_plus$interval_name), FUN = mean)
agr_minus_minus <- aggregate(x = df_minus_minus, by = list(df_minus_minus$interval_name), FUN = mean)
agr_plus_minus <- aggregate(x = df_plus_minus, by = list(df_plus_minus$interval_name), FUN = mean)
agr_right <- aggregate(x = df_sum_right, by = list(df_sum_right$interval_name), FUN = mean)
agr_wrong <- aggregate(x = df_sum_wrong, by = list(df_sum_wrong$interval_name), FUN = mean)


agr_df <- data.frame(interval = agr_plus_plus$interval_name,
                     agr_plus_plus$RPKM, agr_minus_plus$RPKM, agr_minus_minus$RPKM, agr_plus_minus$RPKM,
                     right_direction = agr_right$RPKM, wrong_direction = agr_wrong$RPKM, ratio = agr_right$RPKM / agr_wrong$RPKM)


ggplot(agr_df, aes(x = interval)) +
  geom_vline(xintercept = 20.5, linetype = "dashed", colour = "black", size = 1) +
  geom_point(aes(y = right_direction), col = '#FFD500', size = 5) +
  geom_point(aes(y = wrong_direction), col = 'grey50', size = 5) +
  geom_line(aes(y = right_direction), col = '#FFD500', size = 2) +
  geom_line(aes(y = wrong_direction), col = 'grey50', size = 2) +
  ylab("RPKM") +
  xlab("coordinate relative to Chi-site, kb") +
  ggtitle("") +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20.5, 24, 28, 32, 36, 40),
                     labels = c("", -8, "", -4, "", "Chi-site", "", 4, "", 8, "")) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text = element_text(size = 40, colour = 'black'),
        title = element_text(size = 40),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "white",
                                          linetype = "longdash", 
                                          size = 0.3))

ggsave("Chi-site_metaplot_summary.png", width = 10, height = 7, dpi = 400)
