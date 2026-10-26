library(ggplot2)
library(dplyr)
library(patchwork)

# Import data
data <- read.csv2(
  "Data/Suf/suf_all_mutants.csv"
)

# Process data
data_con <- data %>%
  group_by(Mutant, NO_ext, H2O2_ext, O2, Fe_ext) %>%
  summarise(
    Suf_avg = mean(
      ifelse(Suf == "*", 9, as.numeric(Suf)),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Change_Status = case_when(
      Suf_avg == 9 ~ "~",
      Suf_avg %in% c(0, 1, 2) ~ "",
      TRUE ~ "~"
    )
  ) %>%
  mutate(
    Suf_avg = ifelse(Suf_avg == 9, 1, Suf_avg)
  )

# Function to create individual heatmaps
create_heatmap <- function(subset_data, no_label, h2o2_label) {
  
  ggplot(
    subset_data,
    aes(
      x = factor(O2, levels = c(0, 1)),
      y = factor(Fe_ext, levels = c(1, 0)),
      fill = Suf_avg
    )
  ) +
    geom_tile(color = "black") +
    geom_text(
      aes(label = Change_Status),
      color = "black",
      size = 4
    ) +
    scale_fill_gradientn(
      colors = c("white", "turquoise"),
      limits = c(0, 2),
      na.value = "gray",
      breaks = c(0, 1, 2)
    ) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    coord_fixed(ratio = 1) +
    labs(
      title = paste("(", no_label, ";", h2o2_label, ")"),
      x = "O2",
      y = "Fe_ext",
      fill = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        size = 8,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      plot.margin = margin(3, 3, 3, 3),
      legend.position = "none"
    )
}

# Function to create one compact block per mutant
create_mutant_plot <- function(mutant_name) {
  
  mutant_data <- data_con %>%
    filter(Mutant == mutant_name)
  
  p1 <- create_heatmap(
    filter(mutant_data, NO_ext == 0, H2O2_ext == 0),
    "NO = 0", "H2O2 = 0"
  )
  
  p2 <- create_heatmap(
    filter(mutant_data, NO_ext == 1, H2O2_ext == 0),
    "NO = 1", "H2O2 = 0"
  )
  
  p3 <- create_heatmap(
    filter(mutant_data, NO_ext == 0, H2O2_ext == 1),
    "NO = 0", "H2O2 = 1"
  )
  
  p4 <- create_heatmap(
    filter(mutant_data, NO_ext == 1, H2O2_ext == 1),
    "NO = 1", "H2O2 = 1"
  )
  
  block <- wrap_plots(
    p1, p3,
    p2, p4,
    ncol = 2,
    nrow = 2
  )
  
  block <- block +
    plot_annotation(
      title = mutant_name,
      theme = theme(
        plot.title = element_text(
          size = 13,
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 8)
        )
      )
    )
  
  wrap_elements(full = block)
}

# Get mutant names
mutants <- unique(data_con$Mutant)

# Place wild-type first
mutants <- c(
  "wild-type",
  setdiff(mutants, "wild-type")
)

# Generate all mutant plots
all_plots <- lapply(mutants, create_mutant_plot)

# Function to assemble five mutant blocks with spacing
make_row <- function(plots) {
  
  wrap_plots(
    plots[[1]],
    plot_spacer(),
    plots[[2]],
    plot_spacer(),
    plots[[3]],
    plot_spacer(),
    plots[[4]],
    plot_spacer(),
    plots[[5]],
    ncol = 9,
    widths = c(
      1, 0.25,
      1, 0.25,
      1, 0.25,
      1, 0.25,
      1
    )
  )
}

# Create three rows
row1 <- make_row(all_plots[1:5])
row2 <- make_row(all_plots[6:10])
row3 <- make_row(all_plots[11:15])

# Combine rows with vertical spacing
main_plot <- wrap_plots(
  row1,
  plot_spacer(),
  row2,
  plot_spacer(),
  row3,
  ncol = 1,
  heights = c(1, 0.20, 1, 0.20, 1)
)

# Create one shared legend
legend_plot <- ggplot(
  data.frame(
    x = 1,
    y = 1,
    value = 1
  ),
  aes(x, y, fill = value)
) +
  geom_tile(alpha = 0) +
  scale_fill_gradientn(
    colors = c("white", "turquoise"),
    limits = c(0, 2),
    breaks = c(0, 1, 2),
    name = "Suf"
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth = grid::unit(5, "cm"),
      barheight = grid::unit(0.4, "cm")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Final figure
final_plot <- wrap_plots(
  main_plot,
  legend_plot,
  ncol = 1,
  heights = c(1, 0.08)
)

# Display
print(final_plot)

# Export vector PDF with Unicode support
ggsave(
  filename = "Supplementary/Figures/Suf_all_mutants.pdf",
  plot = final_plot,
  device = cairo_pdf,
  width = 32,
  height = 20,
  limitsize = FALSE
)
