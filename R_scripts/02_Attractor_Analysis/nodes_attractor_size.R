# --- Packages ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# =========================================================
# Input file
# =========================================================

file <- "Data/wild-type.csv"
Data <- read.csv2(file, check.names = FALSE)

# =========================================================
# Data preparation
# =========================================================

Data <- Data %>%
  mutate(
    NO_ext   = as.integer(as.character(NO_ext)),
    H2O2_ext = as.integer(as.character(H2O2_ext)),
    O2       = as.integer(as.character(O2)),
    Fe_ext   = as.integer(as.character(Fe_ext)),
    input_group = paste(NO_ext, H2O2_ext, O2, Fe_ext, sep = "_")
  ) %>%
  mutate(
    across(
      -c(NO_ext, H2O2_ext, O2, Fe_ext, input_group),
      as.character
    )
  )

Data_long <- Data %>%
  pivot_longer(
    cols = -c(input_group, NO_ext, H2O2_ext, O2, Fe_ext),
    names_to = "Node",
    values_to = "Value"
  )

# =========================================================
# Node variation detection
# =========================================================

result_transformed <- Data_long %>%
  group_by(input_group, Node) %>%
  summarise(
    Status = case_when(
      length(unique(Value)) == 1 & !("*" %in% Value) ~
        paste0("Stable (", unique(Value), ")"),
      "*" %in% Value ~
        "(*)",
      TRUE ~
        paste(sort(unique(Value)), collapse = "_")
    ),
    .groups = "drop"
  ) %>%
  mutate(
    is_oscillating = Status %in% c("0_1", "1_2", "0_1_2", "(*)"),
    symbol = ifelse(is_oscillating, "~", "")
  )

# =========================================================
# Node order and labels
# =========================================================

node_order <- c(
  "Fe_free", "Fur", "RyhB",
  "FNR", "FnrS",
  "H2O2", "OxOxyR", "OxyS", "Hpx",
  "NO", "NsrR", "SNOOxyR",
  "ErpA", "NfuA", "IscR", "Isc", "Suf"
)

node_labels <- c(
  "bold(Fe[free])",
  "bold(Fur)",
  "bold(RyhB)",
  "bold(FNR)",
  "bold(FnrS)",
  "bold(H[2]*O[2])",
  "bold(phantom()[Ox]*OxyR)",
  "bold(OxyS)",
  "bold(Hpx)",
  "bold(NO)",
  "bold(NsrR)",
  "bold(phantom()[SNO]*OxyR)",
  "bold(ErpA)",
  "bold(NfuA)",
  "bold(IscR)",
  "bold(Isc)",
  "bold(Suf)"
)

node_pos <- data.frame(
  Node = node_order,
  x_id = seq_along(node_order)
)

# =========================================================
# Module colors
# =========================================================

module_background <- data.frame(
  Module = c(
    "Fe module",
    "FNR module",
    "ROS module",
    "NO module",
    "Fe-S module"
  ),
  xmin = c(0.5, 3.5, 5.5, 9.5, 12.5),
  xmax = c(3.5, 5.5, 9.5, 12.5, 17.5),
  fill = alpha(
    c(
      "#66cc66",
      "#f1c232",
      "#e85c5c",
      "#32b8c8",
      "#9467bd"
    ),
    0.22
  )
)

# =========================================================
# Point fill palette
# =========================================================

fill_vals <- c(
  "Stable (0)" = "#ffffff",
  "Stable (1)" = "#b0b0b0",
  "Stable (2)" = "#000000",
  "0_1"        = "#d9d9d9",
  "0_1_2"      = "#7f7f7f",
  "(*)"        = "#7f7f7f",
  "1_2"        = "#404040"
)

fill_breaks <- c(
  "Stable (0)",
  "Stable (1)",
  "Stable (2)",
  "0_1",
  "0_1_2",
  "(*)",
  "1_2"
)

# =========================================================
# Merge node information
# =========================================================

node_result <- result_transformed %>%
  left_join(node_pos, by = "Node") %>%
  filter(!is.na(x_id))

# =========================================================
# Input order
# =========================================================

row_spacing <- 1.70
panel_gap <- 0.80

input_labels <- Data %>%
  distinct(NO_ext, H2O2_ext, O2, Fe_ext) %>%
  mutate(
    Stress = case_when(
      NO_ext == 0 & H2O2_ext == 0 ~ "NO_0_H2O2_0",
      NO_ext == 0 & H2O2_ext == 1 ~ "NO_0_H2O2_1",
      NO_ext == 1 & H2O2_ext == 0 ~ "NO_1_H2O2_0",
      NO_ext == 1 & H2O2_ext == 1 ~ "NO_1_H2O2_1"
    ),
    Stress = factor(
      Stress,
      levels = c(
        "NO_1_H2O2_1",
        "NO_1_H2O2_0",
        "NO_0_H2O2_1",
        "NO_0_H2O2_0"
      )
    ),
    Fe_label = factor(Fe_ext, levels = c(1, 0)),
    O2_label = factor(O2, levels = c(1, 0))
  ) %>%
  arrange(Stress, Fe_label, O2_label) %>%
  mutate(
    input_group = paste(NO_ext, H2O2_ext, O2, Fe_ext, sep = "_"),
    row_id = n() - row_number() + 1,
    block_id = case_when(
      NO_ext == 0 & H2O2_ext == 0 ~ 1,
      NO_ext == 0 & H2O2_ext == 1 ~ 2,
      NO_ext == 1 & H2O2_ext == 0 ~ 3,
      NO_ext == 1 & H2O2_ext == 1 ~ 4
    ),
    y_plot = row_id * row_spacing + (block_id - 1) * panel_gap
  )

node_result <- node_result %>%
  left_join(
    input_labels %>% select(input_group, row_id, y_plot, block_id),
    by = "input_group"
  )

# =========================================================
# Right-side labels linked to inputs
# =========================================================

attractor_size_original <- Data %>%
  distinct(NO_ext, H2O2_ext, O2, Fe_ext) %>%
  mutate(
    Stress_original = case_when(
      H2O2_ext == 1 & NO_ext == 0 ~ "H2O2_1_NO_0",
      H2O2_ext == 1 & NO_ext == 1 ~ "H2O2_1_NO_1",
      H2O2_ext == 0 & NO_ext == 1 ~ "H2O2_0_NO_1",
      H2O2_ext == 0 & NO_ext == 0 ~ "H2O2_0_NO_0"
    ),
    Stress_original = factor(
      Stress_original,
      levels = c(
        "H2O2_1_NO_0",
        "H2O2_1_NO_1",
        "H2O2_0_NO_1",
        "H2O2_0_NO_0"
      )
    ),
    Fe_label_original = factor(Fe_ext, levels = c(1, 0)),
    O2_label_original = factor(O2, levels = c(1, 0))
  ) %>%
  arrange(Stress_original, Fe_label_original, O2_label_original) %>%
  mutate(
    input_group = paste(NO_ext, H2O2_ext, O2, Fe_ext, sep = "_"),
    row_id_original = n() - row_number() + 1
  ) %>%
  arrange(row_id_original) %>%
  mutate(
    right_label = c(
      "1", "1", "14", "192",
      "1536", "64", "36864", "3072",
      "18432", "768", "442368", "36864",
      "576", "192", "13824", "9216"
    )
  ) %>%
  select(input_group, right_label)

attractor_size <- input_labels %>%
  left_join(attractor_size_original, by = "input_group")

# =========================================================
# Input lines
# =========================================================

make_input_lines <- function(df, var, label, x, color) {
  df2 <- df %>%
    arrange(y_plot) %>%
    filter(.data[[var]] == 1) %>%
    mutate(
      run = cumsum(c(1, diff(row_id) != 1))
    )
  
  df2 %>%
    group_by(run) %>%
    summarise(
      x = x,
      x_text = x - 0.08,
      ymin = min(y_plot) - row_spacing / 2,
      ymax = max(y_plot) + row_spacing / 2,
      y_text = mean(c(min(y_plot), max(y_plot))),
      label = label,
      line_col = color,
      .groups = "drop"
    )
}

input_lines <- bind_rows(
  make_input_lines(
    input_labels,
    "H2O2_ext",
    "bold(H[2]*O[2][ext] == 1)",
    -2.7,
    "#e85c5c"
  ),
  make_input_lines(
    input_labels,
    "NO_ext",
    "bold(NO[ext] == 1)",
    -1.90,
    "#32b8c8"
  ),
  make_input_lines(
    input_labels,
    "Fe_ext",
    "bold(Fe[ext] == 1)",
    -1.0,
    "#66cc66"
  ),
  make_input_lines(
    input_labels,
    "O2",
    "bold(O[2] == 1)",
    -0.14,
    "#f1c232"
  )
)

# =========================================================
# Plot
# =========================================================

p <- ggplot(
  node_result,
  aes(x = x_id, y = y_plot)
) +
  
  geom_rect(
    data = module_background,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -3.0,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = module_background$fill,
    color = NA
  ) +
  
  geom_segment(
    data = input_lines,
    aes(
      x = x,
      xend = x,
      y = ymin,
      yend = ymax
    ),
    inherit.aes = FALSE,
    color = input_lines$line_col,
    linewidth = 1.2
  ) +
  
  geom_text(
    data = input_lines,
    aes(
      x = x_text,
      y = y_text,
      label = label
    ),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 0.5,
    vjust = -0.2,
    size = 3.0,
    color = "black",
    fontface = "bold",
    parse = TRUE
  ) +
  
  geom_text(
    data = node_pos %>%
      mutate(label = node_labels),
    aes(
      x = x_id,
      y = 0.0,
      label = label
    ),
    inherit.aes = FALSE,
    angle = 70,
    hjust = 1,
    vjust = 0.5,
    size = 3.0,
    color = "black",
    fontface = "bold",
    parse = TRUE
  ) +
  
  geom_point(
    aes(fill = Status),
    size = 8,
    shape = 21,
    color = "black",
    stroke = 1
  ) +
  
  scale_fill_manual(
    values = fill_vals,
    breaks = fill_breaks,
    drop = FALSE
  ) +
  
  geom_text(
    aes(label = symbol),
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  
  geom_text(
    data = attractor_size,
    aes(
      x = length(node_order) + 1.0,
      y = y_plot,
      label = right_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    color = "black",
    fontface = "bold"
  ) +
  
  # =========================================================
# Custom legend at bottom
# =========================================================

annotate("text", x = 4.0, y = -4.6, label = "Stable",
         hjust = 1, fontface = "bold", size = 3.2) +
  
  annotate("text", x = 4.0, y = -6.4, label = "Oscillating",
           hjust = 1, fontface = "bold", size = 3.2) +
  
  # Stable row
  annotate("point", x = 5.0, y = -4.6,
           shape = 21, size = 6,
           fill = "white", color = "black", stroke = 1.2) +
  annotate("text", x = 5.65, y = -4.6,
           label = "(0)", hjust = 0,
           size = 2.8, fontface = "bold") +
  
  annotate("point", x = 8.3, y = -4.6,
           shape = 21, size = 6,
           fill = "#b0b0b0", color = "black", stroke = 1.2) +
  annotate("text", x = 8.95, y = -4.6,
           label = "(1)", hjust = 0,
           size = 2.8, fontface = "bold") +
  
  annotate("point", x = 11.6, y = -4.6,
           shape = 21, size = 6,
           fill = "#000000", color = "black", stroke = 1.2) +
  annotate("text", x = 12.25, y = -4.6,
           label = "(2)", hjust = 0,
           size = 2.8, fontface = "bold") +
  
  # Oscillating row
  annotate("point", x = 5.0, y = -6.4,
           shape = 21, size = 6,
           fill = "#d9d9d9", color = "black", stroke = 1.2) +
  annotate("text", x = 5.0, y = -6.4,
           label = "~", size = 5,
           fontface = "bold", color = "black") +
  annotate("text", x = 5.65, y = -6.4,
           label = "(0_1)", hjust = 0,
           size = 2.8, fontface = "bold") +
  
  annotate("point", x = 8.3, y = -6.4,
           shape = 21, size = 6,
           fill = "#7f7f7f", color = "black", stroke = 1.2) +
  annotate("text", x = 8.3, y = -6.4,
           label = "~", size = 5,
           fontface = "bold", color = "black") +
  annotate("text", x = 8.95, y = -6.4,
           label = "(0_max)", hjust = 0,
           size = 2.8, fontface = "bold") +
  
  annotate("point", x = 11.6, y = -6.4,
           shape = 21, size = 6,
           fill = "#404040", color = "black", stroke = 1.2) +
  annotate("text", x = 11.6, y = -6.4,
           label = "~", size = 5,
           fontface = "bold", color = "black") +
  annotate("text", x = 12.25, y = -6.4,
           label = "(1_2)", hjust = 0,
           size = 2.8, fontface = "bold") +
  
  scale_x_continuous(
    breaks = NULL,
    limits = c(-3.0, length(node_order) + 3.0),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 15,
      l = 15
    )
  ) +
  
  coord_cartesian(
    ylim = c(-7.4, max(node_result$y_plot, na.rm = TRUE) + 1),
    clip = "off"
  )

print(p)
ggsave(
  filename = "Figures/Figure_S1.pdf",
  plot = p,
  device = cairo_pdf,
  width = 20,
  height = 19,
  units = "cm",
  bg = "white"
)
ggsave(
  "Figures/Figure_S1.png",
  plot = p,
  width = 20,
  height = 19,
  units = "cm",
  dpi = 600,
  bg = "white"
)