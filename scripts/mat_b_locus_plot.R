#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(patchwork)

# =============================================================================
# FILES & CONFIG
# =============================================================================
REFERENCE <- "fpindikaryon"

GFF_FILE <- file.path("results/annotation", REFERENCE, "annotate_results",
                       paste0("Fomitopsis_pinicola_", REFERENCE, ".gff3"))
HET_FILE <- file.path("external_inputs/assembly_pipeline/results/qc/het_density_hifiasm",
                       REFERENCE, "het_windows.tsv")
OUT_FIG  <- "mat_b_locus.png"

CONTIG       <- "ptg000007l"
REGION_START <- 820000
REGION_END   <- 920000

STE3_GENES <- c("FUN_005089", "FUN_005092", "FUN_005101", "FUN_005102", "FUN_005103")
# Small gene candidate for pheromone precursor (identified manually)
PRECURSOR_CANDIDATES <- c("FUN_005099")

# =============================================================================
# READ GFF
# =============================================================================
gff_cmd <- sprintf(
  "awk '$1==\"%s\" && $3==\"gene\" && $4<%d && $5>%d' %s",
  CONTIG, REGION_END, REGION_START, GFF_FILE
)
gff <- read.delim(pipe(gff_cmd), header = FALSE,
                  col.names = c("chr","src","type","start","end",
                                "score","strand","phase","attr"))
gff$gene_id <- gsub(".*ID=([^;]+);.*", "\\1", gff$attr)
gff$length_bp <- gff$end - gff$start

gff$category <- ifelse(gff$gene_id %in% STE3_GENES, "STE3 pheromone receptor",
                ifelse(gff$gene_id %in% PRECURSOR_CANDIDATES, "Pheromone precursor candidate",
                       "Other gene"))

cat("Genes in region:", nrow(gff), "\n")
cat("STE3 receptors:", sum(gff$category == "STE3 pheromone receptor"), "\n")
cat("Precursor candidates:", sum(gff$category == "Pheromone precursor candidate"), "\n\n")
print(gff[, c("gene_id","start","end","length_bp","strand","category")])

# =============================================================================
# READ HETEROZYGOSITY
# =============================================================================
het <- read.delim(HET_FILE, header = TRUE)
het_region <- het[het$contig == CONTIG &
                  het$window_end > REGION_START &
                  het$window_start < REGION_END, ]

# =============================================================================
# COLOURS
# =============================================================================
cat_cols <- c(
  "STE3 pheromone receptor"       = "#D7301F",
  "Pheromone precursor candidate" = "#FC8D59",
  "Other gene"                    = "#BDBDBD"
)


# =============================================================================
# PANEL A: Heterozygosity comparison across genome regions
# =============================================================================
het_stats <- data.frame(
  region = c("Genome\naverage",
             "MAT-B locus\n(ptg7, 820-920 kb)",
             "ptg10 hot-spot\n(3.4-3.5 Mb)"),
  variants_per_100kb = c(1.21, 0, 225)
)
het_stats$region <- factor(het_stats$region, levels = het_stats$region)

pA <- ggplot(het_stats, aes(x = region, y = variants_per_100kb, fill = region)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.55) +
  geom_text(aes(label = sprintf("%.1f", variants_per_100kb)),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("#999999", "#4DAF4A", "#D7301F"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "A. Heterozygosity context",
       subtitle = "MAT-B shows no heterozygosity (likely collapsed in primary assembly); compare to genome-wide mean and ptg000010l hot-spot",
       x = NULL, y = "Variants per 100 kb") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey30", size = 11),
        axis.text.x = element_text(size = 11, color = "black"))

# =============================================================================
# PANEL B: Gene arrows with strand
# =============================================================================
# Use geom_polygon for arrow shapes
arrow_polygons <- do.call(rbind, lapply(seq_len(nrow(gff)), function(i) {
  g <- gff[i, ]
  tip_frac <- 0.3
  body_h <- 0.25
  head_h <- 0.4
  gene_len <- g$end - g$start
  tip_len <- min(gene_len * tip_frac, 800)

  if (g$strand == "+") {
    tip_x <- g$end - tip_len
    coords <- data.frame(
      x = c(g$start, tip_x, tip_x, g$end,   tip_x, tip_x, g$start),
      y = c(1-body_h, 1-body_h, 1-head_h, 1, 1+head_h, 1+body_h, 1+body_h)
    )
  } else {
    tip_x <- g$start + tip_len
    coords <- data.frame(
      x = c(g$end, tip_x, tip_x, g$start, tip_x, tip_x, g$end),
      y = c(1-body_h, 1-body_h, 1-head_h, 1, 1+head_h, 1+body_h, 1+body_h)
    )
  }
  coords$gene_id <- g$gene_id
  coords$category <- g$category
  coords
}))

# Layer order: other genes first, highlighted last
arrow_polygons$category <- factor(arrow_polygons$category,
                                   levels = c("Other gene",
                                              "Pheromone precursor candidate",
                                              "STE3 pheromone receptor"))

# Label data: STE3 + precursor
label_data <- gff[gff$category != "Other gene", ]
label_data$label_short <- gsub("FUN_00", "", label_data$gene_id)
label_data$y_label <- 1.75

pB <- ggplot() +
  geom_polygon(data = arrow_polygons,
               aes(x = x/1000, y = y, group = gene_id, fill = category),
               color = "black", linewidth = 0.2) +
  # Connector lines
  geom_segment(data = label_data,
               aes(x = (start+end)/2/1000, xend = (start+end)/2/1000,
                   y = 1.42, yend = 1.65),
               color = "grey40", linewidth = 0.3) +
  # Gene labels
  geom_text(data = label_data,
            aes(x = (start+end)/2/1000, y = y_label, label = label_short,
                color = category),
            size = 3.2, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = cat_cols, name = NULL,
                    breaks = c("STE3 pheromone receptor",
                               "Pheromone precursor candidate",
                               "Other gene")) +
  scale_color_manual(values = cat_cols, guide = "none") +
  scale_x_continuous(labels = scales::comma,
                     limits = c(REGION_START/1000, REGION_END/1000),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 2.1), expand = c(0, 0),
                     breaks = 1, labels = "Genes") +
  labs(title = "B. Gene content and orientation",
       subtitle = sprintf("Five STE3 receptors span %.1f kb (%.0f-%.0f kb); arrow direction indicates strand",
                          (max(gff$end[gff$category=="STE3 pheromone receptor"]) -
                           min(gff$start[gff$category=="STE3 pheromone receptor"]))/1000,
                          min(gff$start[gff$category=="STE3 pheromone receptor"])/1000,
                          max(gff$end[gff$category=="STE3 pheromone receptor"])/1000),
       x = paste0("Position on ", CONTIG, " (kb)"), y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey30"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

# =============================================================================
# COMBINE
# =============================================================================
combined <- pA / pB + plot_layout(heights = c(1, 1.3))

ggsave(OUT_FIG, combined, width = 13, height = 7, dpi = 300, bg = "white")
cat("\nSaved figure:", OUT_FIG, "\n")
