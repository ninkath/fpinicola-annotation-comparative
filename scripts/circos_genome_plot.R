#!/usr/bin/env Rscript
# =============================================================================
# circos_landscape.R
#
# Circos-style genome landscape plot for Fomitopsis pinicola.
# Shows the N largest contigs with tracks for GC content, coverage,
# gene density, grouped repeats, and heterozygosity. Telomere-positive
# contig ends are highlighted as red caps on the ideogram.
#
# Usage:
#   Rscript scripts/circos_landscape.R
#
# Adjust INPUT PATHS below to match your pipeline output locations.
# =============================================================================

suppressPackageStartupMessages({
  library(circlize)
})

# =============================================================================
# INPUT PATHS
# =============================================================================
# These paths assume:
#   - Outputs from the assembly/QC pipeline are linked or copied into
#     external_inputs/assembly_pipeline/ in this repository
#   - Outputs from this annotation/comparative pipeline are in results/

REFERENCE <- "fpindikaryon"

# From the assembly/QC pipeline (fpinicola-assembly-qc)
FAI_FILE      <- file.path("external_inputs/assembly_pipeline/results/medaka", REFERENCE, "consensus.fasta.fai")
GC_FILE       <- file.path("external_inputs/assembly_pipeline/results/circos", REFERENCE, "gc_content_10kb.tsv")
COVERAGE_FILE <- file.path("external_inputs/assembly_pipeline/results/circos", REFERENCE, "coverage.regions.bed.gz")
HET_FILE      <- file.path("external_inputs/assembly_pipeline/results/qc/het_density_hifiasm", REFERENCE, "het_windows.tsv")
TELOMERE_FILE <- file.path("external_inputs/assembly_pipeline/results/circos", REFERENCE, "tidk_ttagg_telomeric_repeat_windows.tsv")

# From this pipeline (annotation/comparative)
GENE_FILE   <- file.path("results/circos", REFERENCE, "gene_density_10kb.tsv")
REPEAT_FILE <- file.path("results/circos", REFERENCE, "repeat_bp_density_10kb.tsv")

OUTPUT_FILE <- "genome_circos_landscape.png"
OUTPUT_WIDTH    <- 8.27   # A4 portrait width (inches)
OUTPUT_HEIGHT   <- 11.69  # A4 portrait height (inches)
N_CONTIGS       <- 15

REPEAT_BIN_BP   <- 100000
TELOMERE_THR    <- 20    # 99.9th percentile of internal background

# =============================================================================
# COLOURS
# =============================================================================
IDEO_COLS <- c("#D5E3F0", "#EBF2F8")

GC_COL         <- "#1A9850"   # dark green
GC_MEAN_COL    <- "#C51B7D"   # magenta
COVERAGE_COL   <- "#4D4D4D"   # charcoal grey
GENE_COL       <- "#2C7FB8"   # blue
HET_COL_DARK   <- "#3F007D"   # deep purple

REPEAT_COLS <- c(
  LTR          = "#00A9B5",   # cyan/teal (changed from red)
  DNA_TE       = "#FF7F00",   # orange
  LINE         = "#252525",   # near-black
  Unclassified = "#FEB24C"    # amber
)

TELO_COL <- "#D7301F"   # bright red, reserved for telomere caps

# =============================================================================
# LOAD DATA
# =============================================================================
cat("Reading input files...\n")

fai <- read.delim(FAI_FILE, header = FALSE,
                  col.names = c("chr", "length", "offset", "bases", "bytes"))[, 1:2]
fai <- fai[order(-fai$length), ]
fai <- head(fai, N_CONTIGS)

genome_df <- data.frame(chr = fai$chr, start = 0, end = fai$length)

# GC content
gc <- read.delim(GC_FILE, header = FALSE,
                 col.names = c("chr", "start", "end", "gc"))
gc <- gc[gc$chr %in% fai$chr, ]
if (max(gc$gc, na.rm = TRUE) > 1) gc$gc <- gc$gc / 100
genome_mean_gc <- mean(gc$gc, na.rm = TRUE)
gc_ylim <- c(max(0.25, min(gc$gc, na.rm = TRUE) - 0.02),
             min(0.75, max(gc$gc, na.rm = TRUE) + 0.02))

# Coverage
cov <- read.delim(gzfile(COVERAGE_FILE), header = FALSE,
                  col.names = c("chr", "start", "end", "depth"))
cov <- cov[cov$chr %in% fai$chr, ]
cov_median <- median(cov$depth, na.rm = TRUE)
cov_ylim <- c(0, quantile(cov$depth, 0.995, na.rm = TRUE))

# Gene density
genes <- read.delim(GENE_FILE, header = FALSE,
                    col.names = c("chr", "start", "end", "count"))
genes <- genes[genes$chr %in% fai$chr, ]
gene_max <- max(1, quantile(genes$count, 0.99, na.rm = TRUE))

# Repeats (wide format)
repeats <- read.delim(REPEAT_FILE, header = TRUE)
repeats <- repeats[repeats$chrom %in% fai$chr, ]
for (cls in c("LTR", "DNA_TE", "LINE", "Unclassified", "Other")) {
  repeats[[cls]] <- repeats[[cls]] / 10000
}

aggregate_repeat_track <- function(df, classes, bin_bp) {
  df$bin <- floor(df$start / bin_bp)
  out <- aggregate(df[, classes],
                   by = list(chrom = df$chrom, bin = df$bin),
                   FUN = mean, na.rm = TRUE)
  out$start <- out$bin * bin_bp
  out$end   <- out$start + bin_bp
  contig_len <- setNames(fai$length, fai$chr)
  out$end <- mapply(function(chrom, e) min(e, contig_len[[chrom]]),
                    out$chrom, out$end)
  out[, c("chrom", "start", "end", classes)]
}

repeat_classes <- c("LTR", "DNA_TE", "LINE", "Unclassified")
repeat_agg <- aggregate_repeat_track(repeats, repeat_classes, REPEAT_BIN_BP)
rep_max <- min(1.0, max(0.05,
                max(sapply(repeat_classes,
                    function(c) quantile(repeat_agg[[c]], 0.99, na.rm = TRUE)))))

# Heterozygosity
het <- read.delim(HET_FILE, header = TRUE)
colnames(het)[1:4] <- c("chr", "start", "end", "count")
het <- het[het$chr %in% fai$chr, ]
het_max <- max(1, quantile(het$count, 0.99, na.rm = TRUE))

HET_COL_FUN <- colorRamp2(
  c(0, het_max * 0.2, het_max * 0.5, het_max),
  c("#DADAEB", "#9E9AC8", "#6A51A3", HET_COL_DARK)
)

# Telomere terminal markers (only first/last window above threshold)
telo <- read.delim(TELOMERE_FILE, header = TRUE)
colnames(telo)[1:5] <- c("chr", "window", "forward", "reverse", "motif")
telo$total <- telo$forward + telo$reverse
telo <- telo[telo$chr %in% fai$chr, ]

telo_markers <- do.call(rbind, lapply(fai$chr, function(ctg) {
  d <- telo[telo$chr == ctg, ]
  if (nrow(d) == 0) return(NULL)
  contig_end <- fai$length[fai$chr == ctg]
  first_win <- d[1, ]
  last_win  <- d[nrow(d), ]
  out <- data.frame(chr = character(), end_type = character())
  if (first_win$total >= TELOMERE_THR) {
    out <- rbind(out, data.frame(chr = ctg, end_type = "5'"))
  }
  if (last_win$total >= TELOMERE_THR) {
    out <- rbind(out, data.frame(chr = ctg, end_type = "3'"))
  }
  out
}))

cat(sprintf("Contigs: %d\n", nrow(fai)))
cat(sprintf("Mean GC: %.3f\n", genome_mean_gc))
cat(sprintf("Median coverage: %.1fx\n", cov_median))
cat(sprintf("Telomere markers passing threshold (>=%d): %d\n",
            TELOMERE_THR, nrow(telo_markers)))

# =============================================================================
# PLOT
# =============================================================================
cat("Rendering plot...\n")

ideo_cols <- rep(IDEO_COLS, length.out = nrow(fai))
names(ideo_cols) <- fai$chr

png(OUTPUT_FILE, width = OUTPUT_WIDTH, height = OUTPUT_HEIGHT,
    units = "in", res = 300)

# Top margin larger for title + multi-line description
par(mai = c(1.4, 0, 2.5, 0))

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = rep(1.5, nrow(fai)),
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.006, 0.006),
  points.overflow.warning = FALSE
)

circos.initialize(factors = factor(genome_df$chr, levels = fai$chr),
                  xlim = as.matrix(genome_df[, c("start", "end")]))

# -----------------------------------------------------------------------------
# Track 1: Ideogram with contig labels, Mb ticks, and telomere caps
# -----------------------------------------------------------------------------
circos.track(
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    chr_len <- xlim[2]

    # Base contig rectangle
    circos.rect(xlim[1], 0, xlim[2], 1,
                col = ideo_cols[chr], border = "black", lwd = 0.6)

    # Telomere caps - paint red over the ends that have terminal TTAGG signal
    tmark <- telo_markers[telo_markers$chr == chr, ]
    if (nrow(tmark) > 0) {
      cap_len <- max(chr_len * 0.025, min(chr_len * 0.04, 80000))
      for (i in seq_len(nrow(tmark))) {
        if (tmark$end_type[i] == "5'") {
          circos.rect(0, -0.1, cap_len, 1.1,
                      col = TELO_COL, border = TELO_COL, lwd = 0.8)
        } else {
          circos.rect(chr_len - cap_len, -0.1, chr_len, 1.1,
                      col = TELO_COL, border = TELO_COL, lwd = 0.8)
        }
      }
    }

    # Contig label
    short_name <- gsub("ptg0*([0-9]+)[lc]$", "ptg\\1", chr)
    circos.text(mean(xlim), 0.5, short_name,
                facing = "bending.inside", niceFacing = TRUE,
                cex = 0.7, font = 2, col = "black")

    # Mb ticks
    max_mb <- floor(chr_len / 1e6)
    tick_pos <- c(0, seq(1e6, max_mb * 1e6, by = 1e6))
    tick_labels <- c("0", as.character(seq_len(max_mb)))
    for (i in seq_along(tick_pos)) {
      p <- tick_pos[i]
      circos.segments(p, 1.00, p, 1.28, col = "#333333", lwd = 0.8)
      circos.text(p, 1.42, tick_labels[i],
                  facing = "clockwise", niceFacing = TRUE,
                  adj = c(0.5, 0), cex = 0.55)
    }
  }
)

# -----------------------------------------------------------------------------
# Track 2: GC content
# -----------------------------------------------------------------------------
circos.track(
  ylim = gc_ylim,
  track.height = 0.08,
  bg.border = "black", bg.lwd = 0.4,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    d <- gc[gc$chr == chr, ]
    if (nrow(d) == 0) return()
    mid <- (d$start + d$end) / 2
    circos.lines(mid, d$gc, col = GC_COL, lwd = 1.0)
    circos.segments(CELL_META$xlim[1], genome_mean_gc,
                    CELL_META$xlim[2], genome_mean_gc,
                    col = GC_MEAN_COL, lty = 2, lwd = 1.3)
  }
)

# -----------------------------------------------------------------------------
# Track 3: Coverage
# -----------------------------------------------------------------------------
circos.track(
  ylim = cov_ylim,
  track.height = 0.07,
  bg.border = "black", bg.lwd = 0.4,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    d <- cov[cov$chr == chr, ]
    if (nrow(d) == 0) return()
    mid <- (d$start + d$end) / 2
    circos.lines(mid, pmin(d$depth, cov_ylim[2]),
                 col = COVERAGE_COL, lwd = 0.9)
    circos.segments(CELL_META$xlim[1], cov_median,
                    CELL_META$xlim[2], cov_median,
                    col = "#999999", lty = 3, lwd = 0.8)
  }
)

# -----------------------------------------------------------------------------
# Track 4: Gene density
# -----------------------------------------------------------------------------
circos.track(
  ylim = c(0, gene_max),
  track.height = 0.08,
  bg.border = "black", bg.lwd = 0.4,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    d <- genes[genes$chr == chr, ]
    if (nrow(d) == 0) return()
    circos.rect(d$start, 0, d$end, pmin(d$count, gene_max),
                col = GENE_COL, border = NA)
  }
)

# -----------------------------------------------------------------------------
# Tracks 5-8: Repeats by class
# -----------------------------------------------------------------------------
for (cls in repeat_classes) {
  local({
    current <- cls
    circos.track(
      ylim = c(0, rep_max),
      track.height = 0.065,
      bg.border = "black", bg.lwd = 0.4,
      panel.fun = function(x, y) {
        chr <- CELL_META$sector.index
        d <- repeat_agg[repeat_agg$chrom == chr, ]
        if (nrow(d) == 0) return()
        vals <- d[[current]]
        keep <- !is.na(vals) & vals > 0
        if (!any(keep)) return()
        circos.rect(d$start[keep], 0, d$end[keep],
                    pmin(vals[keep], rep_max),
                    col = REPEAT_COLS[current],
                    border = REPEAT_COLS[current], lwd = 0.4)
      }
    )
  })
}

# -----------------------------------------------------------------------------
# Track 9 (innermost): Heterozygosity
# -----------------------------------------------------------------------------
circos.track(
  ylim = c(0, het_max),
  track.height = 0.05,
  bg.border = "black", bg.lwd = 0.4,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    d <- het[het$chr == chr, ]
    if (nrow(d) == 0) return()
    vals <- pmin(d$count, het_max)
    circos.rect(d$start, 0, d$end, vals,
                col = HET_COL_FUN(vals), border = NA)
  }
)

# =============================================================================
# TITLE, DESCRIPTION, LEGEND
# =============================================================================
title(
  main = expression(bolditalic("Fomitopsis pinicola") ~
                    bold("core genome landscape")),
  cex.main = 1.2, line = 6.0, xpd = TRUE
)


desc <- paste0(
  "Fifteen largest contigs of the primary assembly are shown. Outer ticks ",
  "indicate 1 Mb intervals. Red caps at contig ends mark terminal windows ",
  "with ≥", TELOMERE_THR, " TTAGG repeats per 10 kb, corresponding to recovered ",
  "telomeric sequence. Tracks from outer to inner: GC content (genome mean ",
  sprintf("%.3f", genome_mean_gc),
  " in magenta), read coverage (median ", round(cov_median), "×, grey dotted), ",
  "gene density, grouped repeats (LTR, DNA TE, LINE, Unclassified; mean bp ",
  "coverage per ", REPEAT_BIN_BP / 1000, " kb bin), and heterozygosity ",
  "(variants per 100 kb, scaled to ", round(het_max), ")."
)
desc_lines <- strwrap(desc, width = 105)

for (i in seq_along(desc_lines)) {
  mtext(desc_lines[i], side = 3, line = 5.0 - (i * 0.95),
        cex = 0.78, col = "grey15", xpd = TRUE)
}


# Single unified legend - lines for line-tracks, filled squares for density tracks
legend_items <- c(
  "GC content", "Genome mean GC", "Coverage", "Gene density",
  "Heterozygosity", "Telomere (TTAGG cap)",
  "LTR", "DNA TE", "LINE", "Unclassified"
)
legend_colors <- c(
  GC_COL, GC_MEAN_COL, COVERAGE_COL, GENE_COL,
  HET_COL_DARK, TELO_COL,
  REPEAT_COLS["LTR"], REPEAT_COLS["DNA_TE"],
  REPEAT_COLS["LINE"], REPEAT_COLS["Unclassified"]
)
# Point character: NA for lines, 22 (filled square) for boxes
legend_pch <- c(NA, NA, NA, 22, 22, 22, 22, 22, 22, 22)
# Line type: 1 (solid) for GC/coverage, 2 (dashed) for GC mean, NA for boxes
legend_lty <- c(1, 2, 1, NA, NA, NA, NA, NA, NA, NA)
# Line width: set for lines, NA for boxes
legend_lwd <- c(1.5, 1.5, 1.5, NA, NA, NA, NA, NA, NA, NA)

legend(
  x = 0, y = -1.05, xjust = 0.5, yjust = 1, xpd = TRUE,
  legend = legend_items,
  col = legend_colors,
  pt.bg = legend_colors,
  pch = legend_pch,
  lty = legend_lty,
  lwd = legend_lwd,
  pt.cex = 1.8,
  seg.len = 2.5,
  cex = 0.85, bty = "n", ncol = 3
)

circos.clear()
invisible(dev.off())

cat(sprintf("Saved: %s\n", OUTPUT_FILE))
