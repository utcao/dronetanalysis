#!/usr/bin/env Rscript
# Demo dot plot illustrating an eQTL association (simulated data)
# Gene expression (y) vs SNP genotype (x), one dot per sample

library(ggplot2)
library(dplyr)

set.seed(42)

# --- Simulate DGRP-style inbred line data ---
# 31 lines, each with ~15-25 individual samples
# SNP has two alleles: REF/REF or ALT/ALT (inbred, so no heterozygotes)
n_lines <- 10
samples_per_line <- round(runif(n_lines, 15, 25))

# Assign genotypes to lines: 60% REF/REF, 40% ALT/ALT
genotype_per_line <- sample(
  c("REF/REF", "ALT/ALT"),
  n_lines,
  replace = TRUE,
  prob = c(0.6, 0.4)
)

# eQTL effect: ALT/ALT lines have higher mean expression
expr_mean <- ifelse(genotype_per_line == "REF/REF", 5.2, 6.7)
expr_sd_between <- 2.4   # between-line variance
expr_sd_within  <- 2.3  # within-line variance (individual noise)

# Line-level mean (random effect for line)
line_mean <- rnorm(n_lines, mean = expr_mean, sd = expr_sd_between)

# Expand to individual samples
dat <- lapply(seq_len(n_lines), function(i) {
  n <- samples_per_line[i]
  data.frame(
    line      = paste0("line_", sprintf("%03d", i)),
    genotype  = genotype_per_line[i],
    expr      = pmax(0.1, rnorm(n, mean = line_mean[i], sd = expr_sd_within))
  )
}) |> bind_rows()

dat$genotype <- factor(dat$genotype, levels = c("REF/REF", "ALT/ALT"))

# Additive numeric encoding (0 = REF/REF, 2 = ALT/ALT): standard eQTL model
dat$geno_num <- ifelse(dat$genotype == "REF/REF", 0, 2)

# --- Linear regression (eQTL model) ---
fit   <- lm(expr ~ geno_num, data = dat)
beta  <- coef(fit)["geno_num"]           # effect size per allele copy
r2    <- summary(fit)$r.squared
pval  <- summary(fit)$coefficients["geno_num", "Pr(>|t|)"]
plab       <- ifelse(pval < 1e-300, "p < 1e-300",
               formatC(pval, format = "e", digits = 2))
# plotmath strings (parse = TRUE renders beta and superscript properly)
anno_beta  <- sprintf("beta == %.2f", beta)
anno_r2    <- sprintf("R^2 == %.3f", r2)
anno_p     <- sprintf("p = %s", plab)   # plain text; no parse needed

# Regression line endpoints at factor positions x = 1 and x = 2
pred_ref <- predict(fit, newdata = data.frame(geno_num = 0))
pred_alt <- predict(fit, newdata = data.frame(geno_num = 2))

# --- Summary stats per genotype ---
summary_dat <- dat |>
  group_by(genotype) |>
  summarise(
    mean_expr = mean(expr),
    se        = sd(expr) / sqrt(n()),
    .groups   = "drop"
  )

# --- Colours ---
col_ref <- "#4393C3"
col_alt <- "#D6604D"
pal     <- c("REF/REF" = col_ref, "ALT/ALT" = col_alt)

# --- Plot ---
p <- ggplot(dat, aes(x = genotype, y = expr, colour = genotype)) +

  # Individual dots (jittered)
  geom_jitter(
    width = 0.18, size = 1.8, alpha = 0.55, stroke = 0
  ) +

  # Mean ± SE bar
  geom_errorbar(
    data    = summary_dat,
    aes(x = genotype, ymin = mean_expr - se, ymax = mean_expr + se),
    colour  = "grey20",
    width   = 0.12,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_point(
    data        = summary_dat,
    aes(x = genotype, y = mean_expr),
    size        = 3.5,
    colour      = "grey10",
    inherit.aes = FALSE
  ) +

  # Regression line spanning both genotype positions
  annotate(
    "segment",
    x = 1, xend = 2,
    y = pred_ref, yend = pred_alt,
    colour = "grey20", linewidth = 0.9, linetype = "solid"
  ) +

  # Statistics annotation: three separate lines so parse = TRUE works per line
  annotate("text", x = 1.5, y = Inf,
    label = anno_beta, parse = TRUE,
    hjust = 0.5, vjust = 1.8, size = 4.8, colour = "grey20") +
  annotate("text", x = 1.5, y = Inf,
    label = anno_r2,   parse = TRUE,
    hjust = 0.5, vjust = 3.4, size = 4.8, colour = "grey20") +
  annotate("text", x = 1.5, y = Inf,
    label = anno_p,    parse = FALSE,
    hjust = 0.5, vjust = 7.0, size = 4.8, colour = "grey20") +

  scale_colour_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +

  labs(
    title    = "eQTL example: SNP 2L_5390 × gene A",
    subtitle = paste0(
      "n = ", nrow(dat), " samples |  ",
      "mean ± SE"
    ),
    x        = "SNP genotype",
    y        = "Gene expression"
  ) +

  theme_classic(base_size = 14) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 12),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4)
  )

# --- Save ---
outfile <- "eqtl_demo_dotplot.pdf"
ggsave(outfile, plot = p, width = 5, height = 5.5)
message("Saved: ", outfile)
