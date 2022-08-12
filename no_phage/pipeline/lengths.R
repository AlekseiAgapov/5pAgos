library(ggplot2)

table <- read.table("lengths.tsv", header = T)

ggplot(table, aes(x = length, y = percent)) +
  geom_col(fill = "#B1602A", alpha=0.6, width=0.65) +
  scale_x_continuous(breaks = c(14, 16, 18, 20, 22, 24)) +
  #  ggtitle("smDNA length distribution") +
  ylab("% of all aligned reads") +
  scale_y_continuous(limits = c(0, 30), 
                     breaks = c(0, 10, 20, 30)) +
  xlab("length, nt") +
  theme(text = element_text(size = 30),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text = element_text(size = 40, colour = 'black'),
        title = element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey33",
                                          linetype = "dashed", 
                                          size = 0.5))

ggsave("lengths_distribution.png", width = 7, height = 7, dpi = 300)
