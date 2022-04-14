library(ggplot2)

table <- read.table("lengths.tsv", header = T)

ggplot(table, aes(x = length, y = percent)) +
  geom_col(fill = "#FFD500", alpha=0.5, width=0.65) +
  scale_x_continuous(breaks = c(14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)) +
  ggtitle("smDNA length distribution") +
  ylab("% of all aligned reads") +
  xlab("length, nt") +
  theme(text = element_text(size = 30),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(size = 30, colour = 'black'),
        title = element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "gray33",
                                          linetype = "longdash", 
                                          size = 0.5))

ggsave("lengths_distribution.png", width = 10, height = 10, dpi = 300)
