library(ggplot2)

set.seed(42)
df <- data.frame(
  group = factor(rep(c("Q1", "Q5"), each = 80), levels = c("Q1", "Q5")),
  mad   = c(rnorm(80, mean = 0.3, sd = 0.08),
            rnorm(80, mean = 0.6, sd = 0.10))
)

df2 <- data.frame(
  group = factor(rep(c("Q1", "Q5"), each = 80), levels = c("Q1", "Q5")),
  mad   = c(rnorm(80, mean = 0.6, sd = 0.18),
            rnorm(80, mean = 0.3, sd = 0.09))
)

p_box <- ggplot(df, aes(x = group, y = mad, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.7,
               color = "black", linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("Q1" = "#74ADD1", "Q5" = "#F46D43")) +
  labs(x = "Expression level of Gene A", y = "Transcriptomic MAD") +
  theme_classic(base_size = 22) +
  theme(
    legend.position  = "none",
    axis.text        = element_text(size = 24, color = "black"),
    axis.title.y     = element_text(size = 24, margin = margin(r = 10)),
    axis.title.x     = element_text(size = 16),
    axis.line        = element_line(linewidth = 0.8),
    axis.ticks       = element_line(linewidth = 0.8)
  )

dir.create("results/schematic", showWarnings = FALSE)
ggsave("results/schematic/boxplot_transcriptomic_mad.pdf",
      p_box, width = 4, height = 5)
ggsave("results/schematic/boxplot_transcriptomic_mad.png",
      p_box, width = 4, height = 5, dpi = 300)

p_box2 <- p_box %+% df2 +
  labs(x = "Expression level of Gene B", y = "Transcriptomic MAD")

ggsave("results/schematic/boxplot_transcriptomic_mad_geneB.pdf",
      p_box2, width = 4, height = 5)
ggsave("results/schematic/boxplot_transcriptomic_mad_geneB.png",
      p_box2, width = 4, height = 5, dpi = 300)