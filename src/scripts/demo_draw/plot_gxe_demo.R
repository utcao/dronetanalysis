#!/usr/bin/env Rscript
# Demo reaction-norm plot illustrating a genotype-by-environment (GxE)
# interaction (simulated data).  Two environments: CT (control) and HS
# (heat shock).  Two genotypes: REF/REF and ALT/ALT.
# GxE is visible when the two reaction-norm lines are non-parallel.

library(ggplot2)
library(dplyr)

set.seed(42)

# --- Simulation parameters ---
n_lines         <- 10
samples_per_line <- round(runif(n_lines, 15, 25))
expr_sd_between  <- 0.5   # between-line variance (shared across environments)
expr_sd_within   <- 1.2   # within-group individual noise

# Assign genotypes: 60 % REF/REF, 40 % ALT/ALT (inbred, no heterozygotes)
genotype_per_line <- sample(
  c("REF/REF", "ALT/ALT"), n_lines,
  replace = TRUE, prob = c(0.6, 0.4)
)

# Group means: REF/REF shows a modest HS response;
# ALT/ALT shows a strong HS response → crossing / diverging reaction norms
group_means <- list(
  "REF/REF_CT" = 5.0,
  "REF/REF_HS" = 5.8,   # small increase
  "ALT/ALT_CT" = 4.8,
  "ALT/ALT_HS" = 8.2    # large increase → GxE
)

# Line-level random intercepts (same offset applied in both environments)
line_offset <- rnorm(n_lines, 0, expr_sd_between)

# Expand to individual samples (both environments per line)
dat <- lapply(seq_len(n_lines), function(i) {
  n    <- samples_per_line[i]
  geno <- genotype_per_line[i]
  bind_rows(
    data.frame(
      line      = paste0("line_", sprintf("%03d", i)),
      genotype  = geno,
      env       = "CT",
      expr      = pmax(0, rnorm(n,
                    mean = group_means[[paste0(geno, "_CT")]] + line_offset[i],
                    sd   = expr_sd_within))
    ),
    data.frame(
      line      = paste0("line_", sprintf("%03d", i)),
      genotype  = geno,
      env       = "HS",
      expr      = pmax(0, rnorm(n,
                    mean = group_means[[paste0(geno, "_HS")]] + line_offset[i],
                    sd   = expr_sd_within))
    )
  )
}) |> bind_rows()

dat$genotype <- factor(dat$genotype, levels = c("REF/REF", "ALT/ALT"))
dat$env      <- factor(dat$env,      levels = c("CT", "HS"))

# Numeric encodings for the linear model
dat$geno_num <- ifelse(dat$genotype == "REF/REF", 0, 2)
dat$env_num  <- ifelse(dat$env == "CT", 0, 1)

# --- Linear model with interaction (GxE) ---
fit <- lm(expr ~ geno_num * env_num, data = dat)
cf  <- summary(fit)$coefficients

extract_stat <- function(term) {
  beta <- cf[term, "Estimate"]
  pval <- cf[term, "Pr(>|t|)"]
  plab <- ifelse(pval < 1e-300, "p < 1e-300",
            formatC(pval, format = "e", digits = 2))
  list(beta = beta, plab = plab)
}

s_geno <- extract_stat("geno_num")
s_env  <- extract_stat("env_num")
s_gxe  <- extract_stat("geno_num:env_num")

# Plotmath annotation strings
anno_geno <- sprintf("beta[geno] == %.2f~~(p == '%s')", s_geno$beta, s_geno$plab)
anno_env  <- sprintf("beta[env]  == %.2f~~(p == '%s')", s_env$beta,  s_env$plab)
anno_gxe  <- sprintf("beta[GxE]  == %.2f~~(p == '%s')", s_gxe$beta,  s_gxe$plab)

# --- Per-group summary (mean ± SE) for reaction-norm lines ---
summary_dat <- dat |>
  group_by(genotype, env) |>
  summarise(
    mean_expr = mean(expr),
    se        = sd(expr) / sqrt(n()),
    .groups   = "drop"
  )

# --- Colours and shapes ---
col_ref <- "#4393C3"
col_alt <- "#D6604D"
pal     <- c("REF/REF" = col_ref, "ALT/ALT" = col_alt)

# --- Plot ---
p <- ggplot(dat, aes(x = env, y = expr, colour = genotype)) +

  # Individual dots (jittered within each env × genotype group)
  geom_jitter(
    aes(group = interaction(env, genotype)),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5),
    size = 1.6, alpha = 0.45, stroke = 0
  ) +

  # Error bars on group means
  geom_errorbar(
    data = summary_dat,
    aes(x = env, ymin = mean_expr - se, ymax = mean_expr + se,
        group = genotype),
    position  = position_dodge(width = 0.5),
    width     = 0.1,
    linewidth = 0.8,
    colour    = "grey25",
    inherit.aes = FALSE
  ) +

  # Mean points
  geom_point(
    data = summary_dat,
    aes(x = env, y = mean_expr, group = genotype, colour = genotype),
    position = position_dodge(width = 0.5),
    shape = 18, size = 4.5,
    inherit.aes = FALSE
  ) +

  # Reaction-norm lines connecting CT → HS means per genotype
  geom_line(
    data = summary_dat,
    aes(x = env, y = mean_expr, group = genotype, colour = genotype),
    position  = position_dodge(width = 0.5),
    linewidth = 1.0,
    inherit.aes = FALSE
  ) +

  # GxE statistics annotation (three parse = TRUE lines)
  annotate("text", x = Inf, y = Inf,
    label = anno_geno, parse = TRUE,
    hjust = 1.05, vjust = 1.8,  size = 3.6, colour = "grey20") +
  annotate("text", x = Inf, y = Inf,
    label = anno_env,  parse = TRUE,
    hjust = 1.05, vjust = 3.5,  size = 3.6, colour = "grey20") +
  annotate("text", x = Inf, y = Inf,
    label = anno_gxe,  parse = TRUE,
    hjust = 1.05, vjust = 5.2,  size = 3.6, colour = col_alt,
    fontface = "bold") +

  scale_colour_manual(values = pal, name = "Genotype") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +

  labs(
    title    = "GxE example: SNP 2L_5390 × gene A × environment",
    subtitle = paste0("n = ", nrow(dat), " samples  |  ",
                      n_lines, " DGRP lines  |  diamond = mean ± SE"),
    x        = "Environment",
    y        = "Gene expression"
  ) +

  theme_classic(base_size = 14) +
  theme(
    legend.position    = "inside",
    legend.position.inside = c(0.12, 0.90),
    legend.background  = element_rect(fill = "white", colour = "grey80",
                                      linewidth = 0.3),
    legend.text        = element_text(size = 10),
    legend.title       = element_text(size = 10, face = "bold"),
    legend.key.size    = unit(0.4, "cm"),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, colour = "grey40"),
    axis.title         = element_text(size = 12),
    axis.text          = element_text(size = 12),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4)
  )

# --- Save ---
outfile <- "gxe_demo_dotplot.pdf"
ggsave(outfile, plot = p, width = 5.5, height = 5.5)
message("Saved: ", outfile)
