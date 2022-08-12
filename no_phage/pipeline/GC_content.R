library(ggplot2)

table <- read.table("GC_content.tsv", header = T)
total_genome_GC <- read.table("total_genome_GC.txt", header = F)
total_genome_GC <- total_genome_GC[,1]
total_plasmid_GC <- read.table("total_plasmid_GC.txt", header = F)
total_plasmid_GC <- total_plasmid_GC[,1]

ggplot(table, aes(x = position)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "black", size = 1) +
  geom_hline(yintercept = total_genome_GC, linetype = "dotted", colour = "black", size = 1) +
  geom_point(aes(y = genome), col = 'grey30', size = 4) +
  geom_line(aes(y = genome), col = 'grey30', size = 2) +
  scale_x_continuous(breaks = c(-10, 1, 10, 20, 30)) +
  expand_limits(y = c(30, 70)) +
  scale_y_continuous(breaks = c(30, 40, 50, 60, 70)) +
  ylab("GC, %") +
  xlab("nucleotide position") +
#  ggtitle("genome") +
#  annotate("text", x = c(22, 28, 30), y = 32,
#           label = c('mean GC = ', round(total_genome_GC), '%'), size = 8,
#           col = 'black') +
  theme(text = element_text(size = 30),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("GC_genome.png", width = 7, height = 5, dpi = 300)
