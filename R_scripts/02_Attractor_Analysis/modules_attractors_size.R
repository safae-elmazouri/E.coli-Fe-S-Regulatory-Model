# --- Packages ---
library(ggplot2)
library(dplyr)
library(tidyr)

file <- file.path("Data", "wild-type.csv")
Data <- read.csv2(file, check.names = FALSE)

Data <- Data %>%
  mutate(
    NO_ext   = as.integer(as.character(NO_ext)),
    H2O2_ext = as.integer(as.character(H2O2_ext)),
    O2       = as.integer(as.character(O2)),
    Fe_ext   = as.integer(as.character(Fe_ext)),
    input_group = paste(NO_ext, H2O2_ext, O2, Fe_ext, sep = "_")
  ) %>%
  mutate(across(-c(NO_ext, H2O2_ext, O2, Fe_ext, input_group), as.character))

Data_long <- Data %>%
  pivot_longer(
    cols = -c(input_group, NO_ext, H2O2_ext, O2, Fe_ext),
    names_to = "Node",
    values_to = "Value"
  )

result_transformed <- Data_long %>%
  group_by(input_group, Node) %>%
  summarise(
    Status = case_when(
      length(unique(Value)) == 1 & !("*" %in% Value) ~ paste0("Stable (", unique(Value), ")"),
      "*" %in% Value ~ "Stable (*)",
      TRUE ~ paste(sort(unique(Value)), collapse = "_")
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Status = ifelse(Status == "Stable (*)", "(*)", Status),
    is_oscillating = Status %in% c("0_1", "1_2", "(*)", "0_1_2")
  )

# =========================================================
# Spacing parameters
# =========================================================

column_spacing <- 2.80
row_spacing <- 1.40
panel_gap <- 0.80

# =========================================================
# Modules
# =========================================================

module_order <- c(
  "Fe module",
  "FNR module",
  "ROS module",
  "NO module",
  "Fe-S module"
)

node_to_module <- data.frame(
  Node = c(
    "Fe_free", "Fur", "RyhB",
    "FNR", "FnrS",
    "H2O2", "OxOxyR", "OxyS", "Hpx",
    "NO", "NsrR", "SNOOxyR",
    "Suf", "NfuA", "IscR", "ErpA", "Isc"
  ),
  Module = c(
    rep("Fe module", 3),
    rep("FNR module", 2),
    rep("ROS module", 4),
    rep("NO module", 3),
    rep("Fe-S module", 5)
  )
)

module_pos <- data.frame(
  Module = module_order,
  x_id = c(1.60, 3.40, 5.20, 7.00, 8.80) * column_spacing
)

x_attractor <- max(module_pos$x_id) + 6

module_labels <- c(
  "Iron\nhomeostasis",
  "Anaerobic\nadaptation",
  "Oxidative\nstress\nresponse",
  "Nitrosative\nstress\nresponse",
  "Fe-S\nbiogenesis"
)

module_colors <- data.frame(
  Module = module_order,
  circle_fill = c("#bdf4b8", "#ffe39a", "#ff8f8f", "#9ee5ff", "#e6b4ff")
)

module_text_colors <- data.frame(
  Module = module_order,
  text_col = c("#2ca02c", "#c99700", "#d62728", "#1f77b4", "#9467bd")
)

module_result <- result_transformed %>%
  left_join(node_to_module, by = "Node") %>%
  filter(!is.na(Module)) %>%
  group_by(input_group, Module) %>%
  summarise(
    any_oscillation = any(is_oscillating),
    symbol = ifelse(any_oscillation, "~", ""),
    .groups = "drop"
  ) %>%
  left_join(module_pos, by = "Module") %>%
  left_join(module_colors, by = "Module")

# =========================================================
# Input order
# =========================================================

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

module_result <- module_result %>%
  left_join(
    input_labels %>% select(input_group, row_id, y_plot, block_id),
    by = "input_group"
  )

# =========================================================
# Attractor size
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

make_input_lines <- function(df, var, label, x, color, x_text) {
  df2 <- df %>%
    arrange(y_plot) %>%
    filter(.data[[var]] == 1) %>%
    mutate(run = cumsum(c(1, diff(row_id) != 1)))
  
  df2 %>%
    group_by(run) %>%
    summarise(
      x = x,
      x_text = x_text,
      ymin = min(y_plot) - row_spacing / 2,
      ymax = max(y_plot) + row_spacing / 2,
      y_text = mean(c(min(y_plot), max(y_plot))),
      label = label,
      line_col = color,
      .groups = "drop"
    )
}

input_lines <- bind_rows(
  make_input_lines(input_labels, "H2O2_ext", "H[2]*O[2][ext] == 1", -3.60, "#e85c5c", -4.10),
  make_input_lines(input_labels, "NO_ext", "NO[ext] == 1", -2.40, "#32b8c8", -2.90),
  make_input_lines(input_labels, "Fe_ext", "Fe[ext] == 1", -1.20, "#66cc66", -1.70),
  make_input_lines(input_labels, "O2", "O[2] == 1", 0.00, "#f1c232", -0.50)
)

# =========================================================
# Top labels positions
# =========================================================

y_top <- max(module_result$y_plot) + 2.50
y_module_names <- y_top
y_block_titles <- y_top - 0.10

inputs_title <- data.frame(
  x = mean(c(min(input_lines$x), max(input_lines$x))),
  y = y_block_titles,
  label = ""
)

attractor_title <- data.frame(
  x = x_attractor,
  y = y_block_titles,
  label = "Attractor\nsize"
)

# =========================================================
# Plot
# =========================================================

p <- ggplot(module_result, aes(x = x_id, y = y_plot)) +
  
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
    aes(x = x_text, y = y_text, label = label),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,
    size = 2.8,
    color = "black",
    fontface = "bold",
    parse = TRUE
  ) +
  
  geom_text(
    data = module_pos %>%
      left_join(module_text_colors, by = "Module") %>%
      mutate(label = module_labels),
    aes(
      x = x_id,
      y = y_module_names,
      label = label,
      color = text_col
    ),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    size = 3.4,
    fontface = "bold",
    lineheight = 0.9
  ) +
  
  geom_text(
    data = inputs_title,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    size = 3.6,
    fontface = "bold",
    color = "black"
  ) +
  
  geom_text(
    data = attractor_title,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    size = 3.6,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_color_identity() +
  
  geom_point(
    aes(fill = circle_fill),
    size = 9,
    shape = 21,
    color = "black",
    stroke = 1
  ) +
  
  scale_fill_identity() +
  
  geom_text(
    aes(label = symbol),
    size = 7,
    color = "black"
  ) +
  
  geom_text(
    data = attractor_size,
    aes(
      x = x_attractor,
      y = y_plot,
      label = right_label
    ),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    size = 3.5,
    color = "black",
    fontface = "bold"
  ) +
  
  scale_x_continuous(
    breaks = NULL,
    limits = c(-5.0, x_attractor + 3.0),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = NULL,
    limits = c(-0.3, y_top + 1.0),
    expand = c(0, 0)
  ) +
  
  labs(x = "", y = "", title = "") +
  
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25)
  ) +
  
  coord_cartesian(clip = "off")

print(p)

ggsave(
  filename = file.path("Figures", "Figure_2.pdf"),
  plot = p,
  device = cairo_pdf,
  width = 17,
  height = 18,
  units = "cm",
  limitsize = FALSE
)