library(readr)
library(pheatmap)
library(dplyr)

infile <- snakemake@input[["table"]]
outfile <- snakemake@output[["plot"]]

df <- read_tsv(infile, show_col_types = FALSE)

# Ensure numeric columns
df <- df %>%
  mutate(
    F_pinicola = as.numeric(F_pinicola),
    F_schrenkii = as.numeric(F_schrenkii),
    F_rosea = as.numeric(F_rosea)
  )

# Calculate variability and abundance
df$max_diff <- apply(
  df[, c("F_pinicola", "F_schrenkii", "F_rosea")],
  1,
  function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
)

df$total <- rowSums(
  df[, c("F_pinicola", "F_schrenkii", "F_rosea")],
  na.rm = TRUE
)

# Remove empty families, then select the 25 most variable
top_df <- df %>%
  filter(total > 0) %>%
  arrange(desc(max_diff), desc(total)) %>%
  slice_head(n = 25)

mat <- as.data.frame(top_df[, c("F_pinicola", "F_schrenkii", "F_rosea")])
rownames(mat) <- top_df$Family
colnames(mat) <- c("F. pinicola", "F. schrenkii", "F. rosea")
mat <- as.matrix(mat)

# Integer labels
number_mat <- matrix(
  sprintf("%d", as.integer(mat)),
  nrow = nrow(mat),
  ncol = ncol(mat),
  dimnames = dimnames(mat)
)

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

png(outfile, width = 1450, height = 1700, res = 180)

pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = number_mat,
  number_color = "black",
  border_color = "grey75",
  fontsize_row = 10,
  fontsize_col = 13,
  angle_col = 0,
  labels_col = c(expression(italic("F. pinicola")),
                 expression(italic("F. schrenkii")),
                 expression(italic("F. rosea"))),
  color = colorRampPalette(c("white", "#fee08b", "#f46d43", "#a50026"))(100),
  main = expression(paste("Most variable secreted CAZyme families across ", italic("Fomitopsis"), " species"))
)

dev.off()
