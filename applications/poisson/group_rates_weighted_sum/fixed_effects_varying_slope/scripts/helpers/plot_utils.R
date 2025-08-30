plot_base <- function() {
  
  ggplot() +   
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "#444444", color = NA),  # Dark panel
      plot.background = element_rect(fill = "#2A2A2A", color = NA),  # Dark background
      panel.grid.major = element_line(color = "#4C4C4C"),  # Darker grid
      panel.grid.minor = element_line(color = "#333333"),  
      axis.ticks = element_line(color = "white"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      strip.text = element_text(color = "white"),
      plot.title = element_text(color = "white", face = "bold"),
      legend.background = element_rect(fill = "#444444"),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"))
}

get_LL_plot <- function(df) {
  
  df |>
    ggplot(aes(x = psi, y = .data[[names(df)[2]]])) +
    geom_point(color = "cyan", size = 3, alpha = 0.7) +
    theme_minimal(base_size = 15) +  # Minimal theme with a larger base font size
    theme(
      plot.background = element_rect(fill = "#2E2E2E", color = NA),  # Dark background for the whole plot
      panel.background = element_rect(fill = "#3A3A3F", color = "#1A1A1A", linewidth = 2),  # Lighter panel with a border
      axis.text = element_text(color = "white"),  # White axis labels
      axis.title = element_text(color = "white"),  # White axis titles
      plot.title = element_text(color = "white", size = 18, face = "bold"),  # White title
      plot.caption = element_text(color = "gray", size = 10),  # Gray caption
      panel.grid = element_line(color = "gray30", linetype = "dashed")  # Subtle grid lines
    ) +
    labs(
      title = paste(names(df)[[2]], "Log-Likelihood"),
      x = "\u03C8",
      y = expression("log L("*psi*")")
    )
}

get_branches_plot <- function(mat) {
  
  df <- as.data.frame(mat)
  df$CurveID <- paste0("Curve_", seq_len(nrow(df)))
  
  df_long <- df |> 
    pivot_longer(-CurveID, names_to = "X", values_to = "Y") |> 
    mutate(X = as.numeric(X))
  
  plot_base() + 
    theme(legend.position = "none") +
    geom_line(data = df_long, aes(x = X, y = Y, group = CurveID, color = CurveID), linewidth = 1) +
    labs(title = "Integrated Log-Likelihood Branches", 
         x = "\u03C8",
         y = expression("log L("*psi*")"))
}

get_LL_comparison_plot <- function(stat_fns, 
                                   LL_df_long,
                                   MLE_data,
                                   alpha_levels,
                                   psi_0) {
  
  c(stat_fn_IL, stat_fn_PL) %<-% stat_fns
  
  crit_df <- tibble(
    alpha = alpha_levels,
    crit = qchisq(1 - alpha_levels, df = 1) / 2,
    label = paste0(100 * (1 - alpha_levels), "% CI"),
    color = RColorBrewer::brewer.pal(length(alpha_levels), "Dark2")
  )
  
  palette_colors <- c("Integrated" = "#E41A1C",
                      "Profile" = "#377EB8",
                      "Truth" = "#4DAF4A") |> 
    c(get_ci_colors(crit_df))
  
  psi_endpoints <- LL_df_long |> 
    pull(psi) |> 
    (\(x) c(head(x, 1), tail(x, 1)))()
  
  y_min <- -max(crit_df$crit) - 0.5
  
  label_data <- MLE_data |>
    add_row(Source = "Truth",
            MLE = psi_0,
            MLE_label = "psi[0]")
  
  plot_base() +
    stat_fn_IL +
    stat_fn_PL +
    
    # Baseline
    geom_hline(yintercept = 0, linetype = 5) +
    
    # Confidence interval lines (not mapped to color)
    geom_hline(data = crit_df,
               aes(yintercept = -crit, color = label),
               linewidth = 1.2,
               show.legend = TRUE) +
    
    # Vertical lines at MLEs and Truth
    geom_vline(
      aes(xintercept = MLE, color = Source),
      data = label_data,
      show.legend = FALSE
    ) +
    
    # Repelled MLE/truth labels
    geom_label_repel(
      aes(x = MLE, y = y_min / 2, label = MLE_label, color = Source),
      data = label_data,
      direction = "y",
      parse = TRUE,
      show.legend = FALSE,
      seed = 7835
    ) +
    
    # Log-likelihood points
    geom_point(
      aes(x = psi, y = value - max(value, na.rm = TRUE)),
      data = LL_df_long |> filter(pseudolikelihood == "Integrated"),
      size = 1
    ) +
    
    geom_point(
      aes(x = psi, y = value - max(value, na.rm = TRUE)),
      data = LL_df_long |> filter(pseudolikelihood == "Profile"),
      size = 1
    ) +
    
    # Axes, labels
    ylab(expression("log L("*psi*")")) +
    xlab(expression(psi)) +
    ggtitle("Pseudo Log-Likelihood Comparison Plot") +
    scale_x_continuous(expand = c(0, 0), limits = psi_endpoints) +
    scale_y_continuous(expand = c(0, 0), limits = c(y_min, 0.1)) +
    scale_color_manual(name = NULL, values = palette_colors, breaks = crit_df$label) +
    theme(legend.position = "inside",
          legend.position.inside = c(1, 1),
          legend.justification = c(1, 1))
}

