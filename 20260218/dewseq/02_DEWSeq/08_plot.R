setwd(getwd())

library(ggplot2)
library(dplyr)

data <- read.csv("DEWSeq_results.csv")

head(data)


# Volcano plot
ggplot(data, aes(x = log2FoldChange, y = -log10(p_adj_IHW))) +
  geom_point(aes(color = significant), alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() +
  xlab("log2 Fold Change (IP / IN)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("DEWSeq Volcano Plot (Protein-Coding Only)") +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color="blue") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="blue")

library(ggrepel)
# install.packages("ggrepel")

top_hits <- data %>%
  filter(significant) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(1000)

ggplot(data, aes(x = log2FoldChange, y = -log10(p_adj_IHW))) +
  geom_point(aes(color = significant), alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() +
  geom_text_repel(data = top_hits, aes(label = gene_name), max.overlaps = 25) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color="blue") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="blue")

top_hits$gene_name
