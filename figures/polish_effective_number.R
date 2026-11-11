# Polish the effective-number plot, one request at a time --------------------
#
# The plot: the "effective number of features" per sample (exp of Shannon
# diversity -- how many equally-common features a sample behaves like it has),
# for the species table and the gene table of the instructional dataset.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

source(here::here("R/count_tables.R"))
d <- read_counts("instructional")

dir.create(here::here("figures"), showWarnings = FALSE)

# Same per-sample stats as Step 1 of the notebook. read_counts() returns a
# MATRIX (features in rows, samples in columns, sample IDs as column names),
# so there is no ID column to drop -- use it directly.
per_sample_stats <- function(counts) {
  mat <- as.matrix(counts)
  lib_size <- colSums(mat)
  p <- sweep(mat, 2, lib_size, "/")
  shannon <- -colSums(p * log(p), na.rm = TRUE)
  tibble::tibble(
    sample_id          = colnames(mat),
    library_size       = lib_size,
    effective_features = exp(shannon)
  )
}

per_sample <- bind_rows(
  per_sample_stats(d$species) |> mutate(table = "species"),
  per_sample_stats(d$genes)   |> mutate(table = "genes")
)

# Numbers you can use to check the plot later --------------------------------
cat("---- check these against the plot ----\n")
check <- per_sample |>
  group_by(table) |>
  summarise(n = n(), median = median(effective_features),
            min = min(effective_features), max = max(effective_features))
print(check)
small <- per_sample |> filter(table == "species") |>
  slice_min(library_size, n = 3) |> select(sample_id, library_size, effective_features)
cat("\nthree smallest species libraries:\n"); print(small)

# Axis labels that read 1k, 10k instead of 1e+03 -----------------------------
# (one small function: 1047 -> "1k", 1500 -> "1.5k", 3.23 -> "3.23")
k_labels <- function(x)
  ifelse(x >= 1000,
         paste0(round(x / 1000, 1), "k"),
         format(signif(x, 3), scientific = FALSE, trim = TRUE))

# The two colours I chose: colourblind-safe Okabe-Ito blue and vermillion
species_colour <- "#0072B2"
genes_colour   <- "#D55E00"

# v0: the starting point -- plain defaults ------------------------------------
v0 <- ggplot(per_sample, aes(x = effective_features, y = table, colour = table)) +
  geom_jitter(width = 0, height = 0.07, alpha = 0.5, show.legend = FALSE) +
  scale_x_log10() +
  labs(x = "effective_features", y = "table")
ggsave(here::here("figures/effective-number-v0.png"), v0,
       width = 6, height = 4, dpi = 300)

# v1: request 1 -- axis titles in words, not column names ---------------------
v1 <- v0 + labs(x = "Effective number of features", y = NULL)
ggsave(here::here("figures/effective-number-v1.png"), v1,
       width = 6, height = 4, dpi = 300)

# v2: request 2 -- log-axis labels that read 1k and 10k -----------------------
v2 <- v1 + scale_x_log10(labels = k_labels)
ggsave(here::here("figures/effective-number-v2.png"), v2,
       width = 6, height = 4, dpi = 300)

# v3: request 3 -- two colours of my choice -----------------------------------
v3 <- v2 +
  scale_colour_manual(values = c(species = species_colour, genes = genes_colour))
ggsave(here::here("figures/effective-number-v3.png"), v3,
       width = 6, height = 4, dpi = 300)

# v4: request 4 -- text big enough to read from the back of the room ----------
v4 <- v3 + theme_minimal(base_size = 18)
ggsave(here::here("figures/effective-number-v4.png"), v4,
       width = 6, height = 4, dpi = 300)

# v5: request 5 -- already saved as a 6 by 4 inch PNG -------------------------
# (every ggsave above already used width = 6, height = 4, dpi = 300, so v5
# is simply the re-save under the final name once we are done polishing.)

# v6: the "make it look nicer" request -- deliberately vague ------------------
# What a vague request gets you: a theme swap and recoloured points. Nothing
# is wrong afterwards, but nothing got *better* either -- the plot still says
# nothing about which sample is which, or what the typical value is.
v6 <- v4 +
  theme_bw(base_size = 18) +
  scale_colour_manual(values = c(species = "#1B9E77", genes = "#7570B3"))
ggsave(here::here("figures/effective-number-v6-vague-nicer.png"), v6,
       width = 6, height = 4, dpi = 300)

# v7: the request I would not know how to code myself ------------------------
# "Label the 438-read sample with its ID."
# How it works: ggrepel draws text labels and nudges them apart so they do
# not sit on top of each other or on the points. Docs:
#   https://cran.r-project.org/web/packages/ggrepel/
sample_438 <- per_sample |>
  filter(table == "species", library_size == min(library_size))
cat("\nthe sample being labelled:\n"); print(sample_438)

v7 <- v4 +
  geom_text_repel(
    data = sample_438,
    aes(label = paste0(sample_id, "\n(", library_size, " reads)")),
    size = 4.2, colour = "black", fontface = "bold",
    nudge_y = 0.33,            # lift the label clear of the species strip
    segment.colour = "grey35",
    min.segment.length = 0, seed = 1
  )
ggsave(here::here("figures/effective-number-v7-labelled.png"), v7,
       width = 6, height = 4, dpi = 300)

# v8: "draw the median as a line with its value printed beside it" ------------
# Each label sits just inside its panel edge: right of the line when there is
# room (species), left of it when the line is near the right edge (genes).
medians <- per_sample |>
  group_by(table) |>
  summarise(median = median(effective_features), .groups = "drop") |>
  mutate(
    row       = as.numeric(factor(table, levels = c("genes", "species"))),
    lab_x     = ifelse(median >= 1000, median * 0.93, median * 1.12),
    lab_hjust = ifelse(median >= 1000, 1, 0)
  )

v8 <- v7 +
  # a short dashed segment, confined to its own row, not a full-height line
  geom_segment(
    data = medians,
    aes(x = median, xend = median, y = row - 0.3, yend = row + 0.3,
        colour = table),
    linewidth = 0.9, linetype = "dashed"
  ) +
  geom_label(
    data = medians,
    aes(x = lab_x, y = row, label = paste0("median ", k_labels(median)),
        colour = table, hjust = lab_hjust),
    vjust = 1.7, size = 4.6, fill = "white", linewidth = 0,
    show.legend = FALSE
  ) +
  guides(colour = "none")   # the y axis already says which row is which

ggsave(here::here("figures/effective-number-final.png"), v8,
       width = 6, height = 4, dpi = 300)

cat("\nall versions written to figures/\n")
