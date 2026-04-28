library(readr)
library(dplyr)
library(ComplexUpset)
library(ggplot2)

infile <- snakemake@input[["presence"]]
outfile <- snakemake@output[["plot"]]

df <- read_tsv(infile, show_col_types = FALSE)

df <- df %>%
  mutate(
    F_pinicola = as.integer(F_pinicola),
    F_schrenkii = as.integer(F_schrenkii),
    F_rosea = as.integer(F_rosea)
  ) %>%
  rename(
    "F. pinicola" = F_pinicola,
    "F. schrenkii" = F_schrenkii,
    "F. rosea" = F_rosea
  )

p <- upset(
  df,
  intersect = c("F. pinicola", "F. schrenkii", "F. rosea"),
  name = "Orthogroups"
) +
  labs(
    title = expression(paste("Orthogroup sharing among ", italic("Fomitopsis"), " species"))
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    text = element_text(size = 12)
  )

ggsave(outfile, p, width = 8.4, height = 5.8, dpi = 300)
