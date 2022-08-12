library(ggplot2)
library(ggseqlogo)

logo <- read.table("logo.txt")
logo_v <- as.character(logo[,1])

ggseqlogo(logo_v) +
  ylab("bits") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4),
                     limits = c(0, 0.47)) +
  xlab("nucleotide position") +
  scale_x_continuous(breaks = c(5, 10, 15)) +
  theme(text = element_text(size = 30),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text = element_text(size = 40, colour = 'black'),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())
  
ggsave("logo.png", width = 7, height = 7, dpi = 300)
