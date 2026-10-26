library(ggplot2)
library(dplyr)
library(patchwork)
library(geomtextpath)

base_path <- "Data/Suf"

files <- c(
  "suf (wild-type).csv",
  "suf (ΔfnrS ΔoxyS).csv",
  "suf (ΔryhB ΔoxyS).csv",
  "suf (ΔryhB ΔfnrS).csv"
)

panel_titles <- c(
  "Wild-type",
  "FnrS/OxyS KO",
  "RyhB/OxyS KO",
  "RyhB/FnrS KO"
)

col_H2O2  <- "#E41A1C"
col_NO    <- "#1F5FBF"
col_combo <- "#9C27B0"
col_Fe    <- "#2E9D44"
col_O2    <- "#F2B500"

rings <- data.frame(
  level = c("0", "1", "2"),
  r_inner = c(0.00, 1.00, 2.05),
  r_outer = c(1.00, 2.05, 3.10)
)

make_ring_sector <- function(id, theta1, theta2, r_inner, r_outer, level) {
  data.frame(
    id = id,
    level = level,
    x = c(r_inner*cos(theta1), r_inner*cos(theta2), r_outer*cos(theta2), r_outer*cos(theta1)),
    y = c(r_inner*sin(theta1), r_inner*sin(theta2), r_outer*sin(theta2), r_outer*sin(theta1))
  )
}

make_arc <- function(theta1, theta2, r, n = 120) {
  theta <- seq(theta1, theta2, length.out = n)
  data.frame(x = r * cos(theta), y = r * sin(theta))
}

make_suf_plot <- function(file, title) {
  
  quarter_gap <- 0.06
  quarter_shift <- 0.12
  
  df <- read.csv(file.path(base_path, file), sep = ";") %>%
    mutate(Suf_chr = trimws(as.character(Suf)))
  
  df_plot <- df %>%
    group_by(NO_ext, H2O2_ext, O2, Fe_ext) %>%
    summarise(
      has_star = any(Suf_chr == "*"),
      has_0_raw = any(Suf_chr == "0"),
      has_1_raw = any(Suf_chr == "1"),
      has_2_raw = any(Suf_chr == "2"),
      .groups = "drop"
    ) %>%
    mutate(
      has_0 = has_0_raw | has_star,
      has_1 = has_1_raw | has_star | (has_0_raw & has_2_raw),
      has_2 = has_2_raw | has_star,
      
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
      shift_angle = case_when(
        quarter == 1 ~ pi / 4,
        quarter == 2 ~ -pi / 4,
        quarter == 3 ~ -3 * pi / 4,
        quarter == 4 ~ 3 * pi / 4
      ),
      dx = quarter_shift * cos(shift_angle),
      dy = quarter_shift * sin(shift_angle)
    )
  
  grid_df <- bind_rows(lapply(1:nrow(df_plot), function(i) {
    bind_rows(lapply(1:nrow(rings), function(j) {
      make_ring_sector(
        df_plot$id[i], df_plot$theta1[i], df_plot$theta2[i],
        rings$r_inner[j], rings$r_outer[j], rings$level[j]
      ) %>%
        mutate(x = x + df_plot$dx[i], y = y + df_plot$dy[i])
    }))
  }))
  
  plot_df <- bind_rows(lapply(1:nrow(df_plot), function(i) {
    active_levels <- c()
    if (df_plot$has_0[i]) active_levels <- c(active_levels, "0")
    if (df_plot$has_1[i]) active_levels <- c(active_levels, "1")
    if (df_plot$has_2[i]) active_levels <- c(active_levels, "2")
    
    bind_rows(lapply(active_levels, function(lv) {
      ring <- rings[rings$level == lv, ]
      make_ring_sector(
        df_plot$id[i], df_plot$theta1[i], df_plot$theta2[i],
        ring$r_inner, ring$r_outer, lv
      ) %>%
        mutate(x = x + df_plot$dx[i], y = y + df_plot$dy[i])
    }))
  }))
  
  o2_arc_df <- df_plot %>%
    filter(O2 == 1) %>%
    group_by(Stress_order, Fe_order) %>%
    summarise(
      theta_start = max(theta1),
      theta_end = min(theta2),
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
      r_text = 3.60
    )
  
  o2_arc_lines <- bind_rows(lapply(1:nrow(o2_arc_df), function(i) {
    make_arc(o2_arc_df$theta_start[i], o2_arc_df$theta_end[i], o2_arc_df$r_arc[i], 120) %>%
      mutate(
        x = x + o2_arc_df$dx[i],
        y = y + o2_arc_df$dy[i],
        group = paste(o2_arc_df$Stress_order[i], o2_arc_df$Fe_order[i]),
        arc_color = o2_arc_df$arc_color[i]
      )
  }))
  
  o2_text_lines <- bind_rows(lapply(1:nrow(o2_arc_df), function(i) {
    make_arc(o2_arc_df$theta_start[i], o2_arc_df$theta_end[i], o2_arc_df$r_text[i], 120) %>%
      mutate(
        x = x + o2_arc_df$dx[i],
        y = y + o2_arc_df$dy[i],
        group = paste(o2_arc_df$Stress_order[i], o2_arc_df$Fe_order[i]),
        arc_color = o2_arc_df$arc_color[i],
        arc_label = o2_arc_df$arc_label[i]
      )
  }))
  
  fe_arc_df <- df_plot %>%
    filter(Fe_ext == 1) %>%
    group_by(Stress_order) %>%
    summarise(
      theta_start = max(theta1),
      theta_end = min(theta2),
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
      r_text = 4.30
    )
  
  fe_arc_lines <- bind_rows(lapply(1:nrow(fe_arc_df), function(i) {
    make_arc(fe_arc_df$theta_start[i], fe_arc_df$theta_end[i], fe_arc_df$r_arc[i], 140) %>%
      mutate(
        x = x + fe_arc_df$dx[i],
        y = y + fe_arc_df$dy[i],
        group = fe_arc_df$Stress_order[i],
        arc_color = fe_arc_df$arc_color[i]
      )
  }))
  
  fe_text_lines <- bind_rows(lapply(1:nrow(fe_arc_df), function(i) {
    make_arc(fe_arc_df$theta_start[i], fe_arc_df$theta_end[i], fe_arc_df$r_text[i], 140) %>%
      mutate(
        x = x + fe_arc_df$dx[i],
        y = y + fe_arc_df$dy[i],
        group = fe_arc_df$Stress_order[i],
        arc_color = fe_arc_df$arc_color[i],
        arc_label = fe_arc_df$arc_label[i]
      )
  }))
  
  stress_arc_df <- df_plot %>%
    filter(NO_ext == 1 | H2O2_ext == 1) %>%
    group_by(Stress_order, Stress, Stress_color) %>%
    summarise(
      theta_start = max(theta1),
      theta_end = min(theta2),
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
      r_text = 5.10
    )
  
  stress_arc_lines <- bind_rows(lapply(1:nrow(stress_arc_df), function(i) {
    make_arc(stress_arc_df$theta_start[i], stress_arc_df$theta_end[i], stress_arc_df$r_arc[i], 180) %>%
      mutate(
        x = x + stress_arc_df$dx[i],
        y = y + stress_arc_df$dy[i],
        group = stress_arc_df$Stress_order[i],
        arc_color = stress_arc_df$arc_color[i]
      )
  }))
  
  stress_text_lines <- bind_rows(lapply(1:nrow(stress_arc_df), function(i) {
    make_arc(stress_arc_df$theta_start[i], stress_arc_df$theta_end[i], stress_arc_df$r_text[i], 180) %>%
      mutate(
        x = x + stress_arc_df$dx[i],
        y = y + stress_arc_df$dy[i],
        group = stress_arc_df$Stress_order[i],
        arc_color = stress_arc_df$arc_color[i],
        arc_label = stress_arc_df$arc_label[i]
      )
  }))
  
  p <- ggplot() +
    geom_polygon(
      data = grid_df,
      aes(x = x, y = y, group = interaction(id, level)),
      fill = "white",
      color = "grey75",
      linewidth = 0.25
    ) +
    geom_polygon(
      data = plot_df,
      aes(x = x, y = y, group = interaction(id, level)),
      fill = "#afddc1bd",
      color = "grey60",
      linewidth = 0.25
    ) +
    geom_path(
      data = o2_arc_lines,
      aes(x = x, y = y, group = group, color = arc_color),
      linewidth = 0.55
    ) +
    geom_textpath(
      data = o2_text_lines,
      aes(x = x, y = y, group = group, label = arc_label, color = arc_color),
      size = 1.8,
      parse = TRUE,
      upright = TRUE,
      text_only = TRUE
    ) +
    geom_path(
      data = fe_arc_lines,
      aes(x = x, y = y, group = group, color = arc_color),
      linewidth = 0.55
    ) +
    geom_textpath(
      data = fe_text_lines,
      aes(x = x, y = y, group = group, label = arc_label, color = arc_color),
      size = 1.8,
      parse = TRUE,
      upright = TRUE,
      text_only = TRUE
    ) +
    geom_path(
      data = stress_arc_lines,
      aes(x = x, y = y, group = group, color = arc_color),
      linewidth = 0.65
    ) +
    geom_textpath(
      data = stress_text_lines,
      aes(x = x, y = y, group = group, label = arc_label, color = arc_color),
      size = 2.0,
      parse = TRUE,
      upright = TRUE,
      text_only = TRUE
    ) +
    scale_color_identity() +
    coord_equal(xlim = c(-5.0, 5.0), ylim = c(-5.0, 5.0)) +
    ggtitle(title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold", margin = margin(b = 4)),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 8, 8, 8)
    )
  
  label_shift_x <- quarter_shift * cos(pi / 4)
  label_shift_y <- quarter_shift * sin(pi / 4)
  
  p <- p +
    annotate("text", x = 0.52 + label_shift_x, y = -0.07 + label_shift_y, label = "0", color = "grey60", size = 1.8) +
    annotate("text", x = 1.55 + label_shift_x, y = -0.10 + label_shift_y, label = "1", color = "grey60", size = 1.8) +
    annotate("text", x = 2.58 + label_shift_x, y = -0.10 + label_shift_y, label = "2", color = "grey60", size = 1.8)
  
  return(p)
}

plots <- lapply(seq_along(files), function(k) {
  make_suf_plot(
    file = files[k],
    title = panel_titles[[k]]
  )
})

combined_plot <- wrap_plots(
  plots,
  ncol = 2,
  heights = c(1.3, 1.3)
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.tag = element_text(size = 13, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

combined_plot
ggsave(
  filename = "Figures/Figure_S5.pdf",
  plot = combined_plot,
  device = cairo_pdf,
  width = 17,
  height = 16,
  units = "cm",
  bg = "white"
)
ggsave(
  filename = "Figures/Figure_S5.png",
  plot = combined_plot,
  width = 17,
  height = 16,
  units = "cm",
  dpi = 600,
  bg = "white"
)