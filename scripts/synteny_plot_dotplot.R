#!/usr/bin/env Rscript
# Plot a pairwise NUCmer dotplot on cumulative genome coordinates.
#
# Each retained alignment is drawn as a line segment. Reference and query
# sequences occupy axis space proportional to their physical length within
# each axis. The cross-axis aspect is set by sequence counts rather than
# cumulative length, so cells stay roughly square and the figure fits a
# standing A4 page. Diagonal patterns, inversions, and translocations are
# preserved regardless of axis stretching.

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

table_path  <- snakemake@input[["table"]]
layout_path <- snakemake@input[["layout"]]
png_out     <- snakemake@output[["png"]]
pdf_out     <- snakemake@output[["pdf"]]

ref_name   <- snakemake@wildcards[["ref"]]
query_name <- snakemake@wildcards[["query"]]

dir.create(dirname(png_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(pdf_out), recursive = TRUE, showWarnings = FALSE)

df     <- read_tsv(table_path,  show_col_types = FALSE)
layout <- read_tsv(layout_path, show_col_types = FALSE)

ref_layout   <- layout %>% filter(side == "reference") %>% arrange(order)
query_layout <- layout %>% filter(side == "query")     %>% arrange(order)

ref_total   <- sum(ref_layout$length)
query_total <- sum(query_layout$length)

pretty_species <- function(x) {
  x <- gsub("F_", "F. ", x)
  x <- gsub("_", " ", x)
  x
}

short_label <- function(x) {
  x <- gsub("^>", "", x)
  ifelse(nchar(x) > 16, paste0(substr(x, 1, 13), "..."), x)
}

ref_pretty   <- pretty_species(ref_name)
query_pretty <- pretty_species(query_name)

ref_breaks   <- ref_layout$center
ref_labels   <- short_label(ref_layout$sequence_id)
query_breaks <- query_layout$center
query_labels <- short_label(query_layout$sequence_id)

ref_boundaries   <- c(0, cumsum(ref_layout$length))
query_boundaries <- c(0, cumsum(query_layout$length))

n_x <- length(ref_breaks)
n_y <- length(query_breaks)

# Adaptive X-axis label angle: vertical when crowded, slanted otherwise.
if (n_x > 15) {
  x_angle <- 90
  x_hjust <- 1
  x_vjust <- 0.5
} else {
  x_angle <- 60
  x_hjust <- 1
  x_vjust <- 1
}

size_for_count <- function(n) {
  if (n > 22) 7 else if (n > 14) 8 else 9
}
x_label_size <- size_for_count(n_x)
y_label_size <- size_for_count(n_y)

title_expr <- bquote(
  bold("Pairwise genome alignment: ") *
  bolditalic(.(ref_pretty)) *
  bold(" vs ") *
  bolditalic(.(query_pretty))
)

x_expr <- bquote(italic(.(ref_pretty)) ~ "sequences")
y_expr <- bquote(italic(.(query_pretty)) ~ "sequences")

subtitle_text <- paste0(
  "One-to-one NUCmer alignments (delta-filter -1) at ",
  "\u2265 500 bp and \u2265 80% identity.\n",
  "Forward alignments shown in blue, reverse in purple."
)

p <- ggplot(df) +
  geom_segment(
    aes(x = x_start, xend = x_end,
        y = y_start, yend = y_end,
        colour = orientation),
    linewidth = 1.0,
    alpha = 0.85,
    lineend = "round"
  ) +
  geom_vline(
    xintercept = ref_boundaries,
    colour = "grey75",
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = query_boundaries,
    colour = "grey75",
    linewidth = 0.35
  ) +
  scale_x_continuous(
    breaks = ref_breaks,
    labels = ref_labels,
    limits = c(0, ref_total),
    expand = c(0, 0)
  ) +
  scale_y_reverse(
    breaks = query_breaks,
    labels = query_labels,
    limits = c(query_total, 0),
    expand = c(0, 0)
  ) +
  scale_colour_manual(
    values = c("Forward" = "#56B4E9", "Reverse" = "#8A2BE2"),
    guide = "none"
  ) +
  labs(
    title = title_expr,
    subtitle = subtitle_text,
    x = x_expr,
    y = y_expr
  ) +
  theme_grey(base_size = 13) +
  theme(
    plot.title.position = "plot",
    plot.title          = element_text(size = 14, hjust = 0.5,
                                       margin = margin(b = 4)),
    plot.subtitle       = element_text(size = 10, hjust = 0.5,
                                       margin = margin(b = 8)),
    panel.background    = element_rect(fill = "grey95", colour = NA),
    axis.title.x        = element_text(size = 13, margin = margin(t = 8)),
    axis.title.y        = element_text(size = 13, margin = margin(r = 8)),
    axis.text.x         = element_text(angle = x_angle,
                                       hjust = x_hjust,
                                       vjust = x_vjust,
                                       size = x_label_size),
    axis.text.y         = element_text(size = y_label_size),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank()
  )

# Canvas size driven by sequence counts. Cells stay roughly square; clamp
# to a standing-A4 friendly range.
target_cell_w <- 0.30
target_cell_h <- 0.36
margin_w <- 2.2
margin_h <- 2.2

fig_w <- min(max(n_x * target_cell_w + margin_w, 6.5), 7.5)
fig_h <- min(max(n_y * target_cell_h + margin_h, 5.5), 9.0)

ggsave(png_out, p, width = fig_w, height = fig_h, dpi = 300)
ggsave(pdf_out, p, width = fig_w, height = fig_h)
