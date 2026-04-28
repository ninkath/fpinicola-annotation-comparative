library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)

infile <- snakemake@input[["table"]]
outfile <- snakemake@output[["plot"]]

df <- read_tsv(infile, show_col_types = FALSE)

# Keep only actual COG rows, remove summary rows, and drop tiny categories W and N
df_plot <- df %>%
  filter(!is.na(COG_category)) %>%
  filter(!COG_category %in% c("TOTAL_PROTEINS", "COG_ANNOTATED", "ANNOTATION_RATE", "W", "N")) %>%
  filter(COG_category != "") %>%
  mutate(
    F_pinicola = as.numeric(F_pinicola),
    F_schrenkii = as.numeric(F_schrenkii),
    F_rosea = as.numeric(F_rosea)
  ) %>%
  pivot_longer(
    cols = c("F_pinicola", "F_schrenkii", "F_rosea"),
    names_to = "Species",
    values_to = "Count"
  ) %>%
  mutate(
    Species = recode(
      Species,
      "F_pinicola" = "F. pinicola",
      "F_schrenkii" = "F. schrenkii",
      "F_rosea" = "F. rosea"
    )
  ) %>%
  group_by(COG_category, Description) %>%
  mutate(Total = sum(Count, na.rm = TRUE)) %>%
  ungroup()

# Order categories by total count across species
category_order <- df_plot %>%
  distinct(COG_category, Total) %>%
  arrange(desc(Total)) %>%
  pull(COG_category)

df_plot$COG_category <- factor(df_plot$COG_category, levels = category_order)
df_plot$Species <- factor(df_plot$Species, levels = c("F. pinicola", "F. schrenkii", "F. rosea"))

p <- ggplot(df_plot, aes(x = COG_category, y = Count, fill = Species)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  labs(
    title = expression(paste("COG category profiles across ", italic("Fomitopsis"), " species")),
    x = "COG category",
    y = "Protein count",
    fill = NULL
  ) +
  scale_fill_manual(
    values = c(
      "F. pinicola" = "#1f77ff",
      "F. schrenkii" = "#e31a1c",
      "F. rosea" = "#f4a300"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(face = "italic"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave(outfile, plot = p, width = 8.8, height = 5.6, dpi = 300)
