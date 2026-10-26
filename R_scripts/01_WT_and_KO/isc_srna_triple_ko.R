# Figure S4 — Isc activity in wild-type and OxyS/FnrS/RyhB knockout
# An asterisk (*) is interpreted as allowing both Isc levels (0 and 1).

library(ggplot2)
library(dplyr)
library(patchwork)

base_path_isc <- file.path("Data", "Isc")
output_dir <- "Figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

file_wt  <- file.path(base_path_isc, "isc (wild_type).csv")
file_mut <- file.path(base_path_isc, "isc (oxys_fnrs_ryhb_ko).csv")

col_H2O2  <- "#E41A1C"
col_NO    <- "#1F5FBF"
col_combo <- "#9C27B0"
col_Fe    <- "#2E9D44"
col_O2    <- "#F2B500"

rings_isc <- data.frame(
  level = c("0", "1"),
  r_inner = c(0.00, 1.55),
  r_outer = c(1.55, 3.10)
)

make_ring_sector <- function(id, theta1, theta2, r_inner, r_outer, level) {
  data.frame(
    id = id,
    level = level,
    x = c(
      r_inner * cos(theta1),
      r_inner * cos(theta2),
      r_outer * cos(theta2),
      r_outer * cos(theta1)
    ),
    y = c(
      r_inner * sin(theta1),
      r_inner * sin(theta2),
      r_outer * sin(theta2),
      r_outer * sin(theta1)
    )
  )
}

make_arc <- function(theta1, theta2, r, n = 100) {
  theta <- seq(theta1, theta2, length.out = n)
  data.frame(
    x = r * cos(theta),
    y = r * sin(theta)
  )
}

make_label_angle <- function(theta_mid) {
  angle <- theta_mid * 180 / pi + 90
  angle <- angle %% 360
  angle <- ifelse(angle > 90 & angle < 270, angle + 180, angle)
  angle %% 360
}

make_hjust <- function(theta_mid) {
  0.5
}

make_radial_plot <- function(file, node, title, rings, fill_color) {
  
  quarter_gap <- 0.06
  quarter_shift <- 0.12
  
  df <- read.csv(file, sep = ";")
  
  df <- df %>%
    mutate(value_chr = trimws(as.character(.data[[node]])))
  
  df_plot <- df %>%
    group_by(NO_ext, H2O2_ext, O2, Fe_ext) %>%
    summarise(
      has_star = any(value_chr == "*"),
      has_0_raw = any(value_chr == "0"),
      has_1_raw = any(value_chr == "1"),
      .groups = "drop"
    ) %>%
    mutate(
      has_0 = has_0_raw | has_star,
      has_1 = has_1_raw | has_star,
      
      Stress = case_when(
        H2O2_ext == 1 & NO_ext == 0 ~ "H[2]*O[2][ext]==1",
        H2O2_ext == 0 & NO_ext == 1 ~ "NO[ext]==1",
        H2O2_ext == 1 & NO_ext == 1 ~ "H[2]*O[2][ext]~~'='~~NO[ext]~~'='~~1",
        TRUE ~ ""
      ),
      
      Stress_color = case_when(
        H2O2_ext == 1 & NO_ext == 0 ~ col_H2O2,
        H2O2_ext == 0 & NO_ext == 1 ~ col_NO,
        H2O2_ext == 1 & NO_ext == 1 ~ col_combo,
        TRUE ~ NA_character_
      ),
      
      Stress_order = case_when(
        H2O2_ext == 1 & NO_ext == 0 ~ 1,
        H2O2_ext == 1 & NO_ext == 1 ~ 2,
        H2O2_ext == 0 & NO_ext == 1 ~ 3,
        H2O2_ext == 0 & NO_ext == 0 ~ 4
      ),
      
      Fe_order = ifelse(Fe_ext == 1, 1, 2),
      O2_order = ifelse(O2 == 1, 1, 2)
    ) %>%
    arrange(Stress_order, Fe_order, O2_order) %>%
    mutate(
      id = row_number(),
      quarter = ceiling(id / 4),
      
      theta1 = pi / 2 - 2 * pi * (id - 1) / 16,
      theta2 = pi / 2 - 2 * pi * id / 16,
      
      theta1 = theta1 - quarter_gap / 2,
      theta2 = theta2 + quarter_gap / 2,
      theta_mid = (theta1 + theta2) / 2,
      
      shift_angle = case_when(
        quarter == 1 ~ pi / 4,
        quarter == 2 ~ -pi / 4,
        quarter == 3 ~ -3 * pi / 4,
        quarter == 4 ~ 3 * pi / 4
      ),
      
      dx = quarter_shift * cos(shift_angle),
      dy = quarter_shift * sin(shift_angle)
    )
  
  # Build the two-level radial grid for each environmental condition.
  grid_df <- bind_rows(lapply(1:nrow(df_plot), function(i) {
    bind_rows(lapply(1:nrow(rings), function(j) {
      make_ring_sector(
        id = df_plot$id[i],
        theta1 = df_plot$theta1[i],
        theta2 = df_plot$theta2[i],
        r_inner = rings$r_inner[j],
        r_outer = rings$r_outer[j],
        level = rings$level[j]
      ) %>%
        mutate(
          x = x + df_plot$dx[i],
          y = y + df_plot$dy[i]
        )
    }))
  }))
  

    plot_df <- bind_rows(lapply(1:nrow(df_plot), function(i) {
    active_levels <- c()
    
    if (df_plot$has_0[i]) active_levels <- c(active_levels, "0")
    if (df_plot$has_1[i]) active_levels <- c(active_levels, "1")
    
    bind_rows(lapply(active_levels, function(lv) {
      ring <- rings[rings$level == lv, ]
      
      make_ring_sector(
        id = df_plot$id[i],
        theta1 = df_plot$theta1[i],
        theta2 = df_plot$theta2[i],
        r_inner = ring$r_inner,
        r_outer = ring$r_outer,
        level = lv
      ) %>%
        mutate(
          x = x + df_plot$dx[i],
          y = y + df_plot$dy[i]
        )
    }))
  }))
  

  o2_arc_df <- df_plot %>%
    filter(O2 == 1) %>%
    group_by(Stress_order, Fe_order) %>%
    summarise(
      theta_start = max(theta1),
      theta_end = min(theta2),
      theta_mid = mean(c(theta_start, theta_end)),
      dx = mean(dx),
      dy = mean(dy),
      .groups = "drop"
    ) %>%
    mutate(
      gap = 0.02,
      theta_start = theta_start - gap,
      theta_end = theta_end + gap,
      arc_label = "O[2]==1",
      arc_color = col_O2,
      r_arc = 3.30,
      r_text = 3.60,
      x_text = r_text * cos(theta_mid) + dx,
      y_text = r_text * sin(theta_mid) + dy,
      angle_text = make_label_angle(theta_mid),
      hjust_text = make_hjust(theta_mid)
    )
  
  o2_arc_lines <- bind_rows(lapply(1:nrow(o2_arc_df), function(i) {
    make_arc(o2_arc_df$theta_start[i], o2_arc_df$theta_end[i], o2_arc_df$r_arc[i], 80) %>%
      mutate(
        x = x + o2_arc_df$dx[i],
        y = y + o2_arc_df$dy[i],
        group = paste(o2_arc_df$Stress_order[i], o2_arc_df$Fe_order[i]),
        arc_color = o2_arc_df$arc_color[i]
      )
  }))
  
  fe_arc_df <- df_plot %>%
    filter(Fe_ext == 1) %>%
    group_by(Stress_order) %>%
    summarise(
      theta_start = max(theta1),
      theta_end = min(theta2),
      theta_mid = mean(c(theta_start, theta_end)),
      dx = mean(dx),
      dy = mean(dy),
      .groups = "drop"
    ) %>%
    mutate(
      gap = 0.06,
      theta_start = theta_start - gap,
      theta_end = theta_end + gap,
      arc_label = "Fe[ext]==1",
      arc_color = col_Fe,
      r_arc = 4.00,
      r_text = 4.30,
      x_text = r_text * cos(theta_mid) + dx,
      y_text = r_text * sin(theta_mid) + dy,
      angle_text = make_label_angle(theta_mid),
      hjust_text = make_hjust(theta_mid)
    )
  
  fe_arc_lines <- bind_rows(lapply(1:nrow(fe_arc_df), function(i) {
    make_arc(fe_arc_df$theta_start[i], fe_arc_df$theta_end[i], fe_arc_df$r_arc[i], 80) %>%
      mutate(
        x = x + fe_arc_df$dx[i],
        y = y + fe_arc_df$dy[i],
        group = fe_arc_df$Stress_order[i],
        arc_color = fe_arc_df$arc_color[i]
      )
  }))
  
  stress_arc_df <- df_plot %>%
    filter(NO_ext == 1 | H2O2_ext == 1) %>%
    group_by(Stress_order, Stress, Stress_color) %>%
    summarise(
      theta_start = max(theta1),
      theta_end = min(theta2),
      theta_mid = mean(c(theta_start, theta_end)),
      dx = mean(dx),
      dy = mean(dy),
      .groups = "drop"
    ) %>%
    mutate(
      gap = 0.06,
      theta_start = theta_start - gap,
      theta_end = theta_end + gap,
      arc_label = as.character(Stress),
      arc_color = Stress_color,
      r_arc = 4.80,
      r_text = 5.10,
      x_text = r_text * cos(theta_mid) + dx,
      y_text = r_text * sin(theta_mid) + dy,
      angle_text = make_label_angle(theta_mid),
      hjust_text = make_hjust(theta_mid)
    )
  
  stress_arc_lines <- bind_rows(lapply(1:nrow(stress_arc_df), function(i) {
    make_arc(stress_arc_df$theta_start[i], stress_arc_df$theta_end[i], stress_arc_df$r_arc[i], 160) %>%
      mutate(
        x = x + stress_arc_df$dx[i],
        y = y + stress_arc_df$dy[i],
        group = stress_arc_df$Stress_order[i],
        arc_color = stress_arc_df$arc_color[i]
      )
  }))
  
  p <- ggplot() +
    geom_polygon(
      data = grid_df,
      aes(x = x, y = y, group = interaction(id, level)),
      fill = "white",
      color = "grey75",
      linewidth = 0.30
    ) +
    geom_polygon(
      data = plot_df,
      aes(x = x, y = y, group = interaction(id, level)),
      fill = fill_color,
      color = "grey55",
      linewidth = 0.30
    ) +
    geom_path(
      data = o2_arc_lines,
      aes(x = x, y = y, group = group, color = arc_color),
      linewidth = 0.80
    ) +
    geom_text(
      data = o2_arc_df,
      aes(
        x = x_text,
        y = y_text,
        label = arc_label,
        angle = angle_text,
        color = arc_color,
        hjust = hjust_text
      ),
      size = 2.5,
      parse = TRUE
    ) +
    geom_path(
      data = fe_arc_lines,
      aes(x = x, y = y, group = group, color = arc_color),
      linewidth = 0.80
    ) +
    geom_text(
      data = fe_arc_df,
      aes(
        x = x_text,
        y = y_text,
        label = arc_label,
        angle = angle_text,
        color = arc_color,
        hjust = hjust_text
      ),
      size = 2.5,
      parse = TRUE
    ) +
    geom_path(
      data = stress_arc_lines,
      aes(x = x, y = y, group = group, color = arc_color),
      linewidth = 0.80
    ) +
    geom_text(
      data = stress_arc_df,
      aes(
        x = x_text,
        y = y_text,
        label = arc_label,
        angle = angle_text,
        color = arc_color,
        hjust = hjust_text
      ),
      size = 2.5,
      parse = TRUE
    ) +
    scale_color_identity() +
    coord_equal(xlim = c(-5.0, 5.0), ylim = c(-5.0, 5.0)) +
    ggtitle(title) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 10,
        face = "bold",
        margin = margin(b = 2)
      ),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  label_shift_x <- quarter_shift * cos(pi / 4)
  label_shift_y <- quarter_shift * sin(pi / 4)
  
  p <- p +
    annotate("text", x = 0.77 + label_shift_x, y = -0.07 + label_shift_y, label = "0", color = "grey60", size = 1.8) +
    annotate("text", x = 2.32 + label_shift_x, y = -0.10 + label_shift_y, label = "1", color = "grey60", size = 1.8)
  
  return(p)
}

p_wt <- make_radial_plot(
  file = file_wt,
  node = "Isc",
  title = "Wild-type",
  rings = rings_isc,
  fill_color = "#f6c1a9bd"
)

p_mut <- make_radial_plot(
  file = file_mut,
  node = "Isc",
  title = "OxyS/FnrS/RyhB KO",
  rings = rings_isc,
  fill_color = "#f6c1a9bd"
)

combined_plot <- p_wt + p_mut +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.tag = element_text(size = 13, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

combined_plot
ggsave(
  filename = file.path(output_dir, "Figure_S4.pdf"),
  plot = combined_plot,
  device = cairo_pdf,
  width = 17,
  height = 11.5,
  units = "cm",
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Figure_S4.png"),
  plot = combined_plot,
  width = 17,
  height = 11.5,
  units = "cm",
  dpi = 600,
  bg = "white"
)
