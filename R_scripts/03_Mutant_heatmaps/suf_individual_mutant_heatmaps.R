library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
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
      x = factor(O2),
      y = factor(Fe_ext),
      fill = Suf_avg
    )
  ) +
    geom_tile(color = "black") +
    geom_text(
      aes(label = Change_Status),
      color = "black",
      size = 6
    ) +
    scale_fill_gradientn(
      colors = c("white", "turquoise"),
      limits = c(0, 2),
      na.value = "gray",
      breaks = c(0, 1, 2)
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(
      limits = rev(levels(factor(subset_data$Fe_ext)))
    ) +
    labs(
      title = paste("(", no_label, ";", h2o2_label, ")"),
      x = "O2",
      y = "Fe_ext",
      fill = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      strip.text = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
}


# Function to generate four heatmaps for each mutant
create_mutant_plot <- function(mutant_name) {
  
  mutant_data <- data_con %>%
    filter(Mutant == mutant_name)
  
  plot1 <- create_heatmap(
    filter(mutant_data, NO_ext == 0, H2O2_ext == 0),
    "NO = 0", "H2O2 = 0"
  )
  
  plot2 <- create_heatmap(
    filter(mutant_data, NO_ext == 1, H2O2_ext == 0),
    "NO = 1", "H2O2 = 0"
  )
  
  plot3 <- create_heatmap(
    filter(mutant_data, NO_ext == 0, H2O2_ext == 1),
    "NO = 0", "H2O2 = 1"
  )
  
  plot4 <- create_heatmap(
    filter(mutant_data, NO_ext == 1, H2O2_ext == 1),
    "NO = 1", "H2O2 = 1"
  )
  
  (plot1 | plot3) /
    (plot2 | plot4) +
    plot_layout(
      guides = "collect",
      heights = c(2, 2)
    ) +
    plot_annotation(
      title = mutant_name,
      theme = theme(
        plot.title = element_text(
          size = 14,
          face = "bold",
          hjust = 0.5
        )
      )
    )
}


# Get all mutant names
mutants <- unique(data_con$Mutant)


# Output directory
output_dir <- "Supplementary/Figures/Suf"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Generate and save all mutant figures
for (mutant in mutants) {
  
  p <- create_mutant_plot(mutant)
  
  filename <- paste0("Suf_", mutant, ".pdf")
  
  ggsave(
    filename = file.path(output_dir, filename),
    plot = p,
    device = cairo_pdf,
    width = 8,
    height = 6
  )
}
