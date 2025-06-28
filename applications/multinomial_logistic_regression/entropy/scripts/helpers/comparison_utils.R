
get_LL_df_long <- function(LL_df) {
  
  LL_df |>
    pivot_longer(cols = -psi,
                 names_to = "pseudolikelihood",
                 values_to = "value") |>
    mutate(pseudolikelihood = pseudolikelihood |>
             as_factor())
}

get_spline_models <- function(LL_df_long) {
  
  LL_df_long |>
    drop_na(value) |>
    group_by(pseudolikelihood) |>
    group_map(~ smooth.spline(.x$psi, .x$value, spar = 0.3)) |>
    set_names(levels(LL_df_long$pseudolikelihood))
}

get_MLE_data <- function(spline_models, LL_df_long) {
  
  spline_models |>
    imap(
      \(mod, pseudolikelihood) {
        optimize(\(psi) predict(mod, psi)$y,
                 lower = LL_df_long |>
                   filter(pseudolikelihood == pseudolikelihood) |> 
                   select(psi) |>
                   min(),
                 upper = LL_df_long |>
                   filter(pseudolikelihood == pseudolikelihood) |> 
                   select(psi) |>
                   max(),
                 maximum = TRUE) |>
          data.frame() |>
          mutate(MLE = as.numeric(maximum),
                 Maximum = as.numeric(objective)) |>
          select(MLE, Maximum)
      }) |>
    do.call(rbind, args = _) |>
    mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[PL]")) |>
    rownames_to_column("Source")
}

get_pseudolikelihoods <- function(spline_models, MLE_data) {
  
  spline_models |> 
    map2(MLE_data$Maximum,
         \(mod, maximum) \(psi) predict(mod, psi)$y - maximum)
}

get_confidence_intervals <- function(pseudolikelihoods, 
                                     LL_df_long,
                                     MLE_data, 
                                     alpha_levels, 
                                     J) {
  
  map2_dfr(
    pseudolikelihoods,
    names(pseudolikelihoods),
    \(pseudolikelihood, name) {
      
      MLE <- MLE_data |> 
        filter(Source == name) |> 
        pull(MLE)
      
      psi_endpoints <- LL_df_long |> 
        filter(pseudolikelihood == name) |> 
        drop_na() |> 
        pull(psi) |> 
        (\(x) c(head(x, 1), tail(x, 1)))()
      
      map_dfr(alpha_levels, \(alpha) {
        crit <- qchisq(1 - alpha, df = 1) / 2
        conf_label <- paste0(100 * (1 - alpha), "%")
        
        lower_bound <- tryCatch(
          uniroot(\(psi) pseudolikelihood(psi) + crit,
                  interval = c(psi_endpoints[1], MLE))$root,
          error = function(e) return(0)
        ) |> round(3)
        
        upper_bound <- tryCatch( 
          uniroot(\(psi) pseudolikelihood(psi) + crit,
                  interval = c(MLE, psi_endpoints[2]))$root,
          error = function(e) return(log(J))
        ) |> round(3)
        
        tibble(
          pseudolikelihood = name,
          confidence = conf_label,
          alpha = alpha,
          lower = lower_bound,
          upper = upper_bound
        )
      })
    }
  )
}

make_kable_caption <- function(path) {

  path_parts <- strsplit(path, "/")[[1]]
  
  experiment_version <- path_parts[grep("^exp_v", path_parts)]
  sim_num <- path_parts[grep("^sim_", path_parts)]
  iter_num <- path_parts[grep("^iter_", path_parts)]

  subtitle <- paste("Experiment", experiment_version,
                    "— Simulation", sim_num,
                    "— Iteration", iter_num)
  title = "Pseudolikelihood Estimate Comparison Tables"
  
  caption = paste0("<center><span style='font-size:100%'>",
                   title,
                   "</span><br><span style='font-size:50%'>", 
                   subtitle,
                   "</span></center>")
}

render_LL_comparison_tables <- function(MLE_data, conf_ints, psi_0) {
  
  # Format MLE table
  mle_table <- MLE_data |>
    select(Source, MLE) |> 
    mutate(MLE = sprintf("%.3f", MLE)) |>
    bind_rows(tibble(Source = "Truth", MLE = sprintf("%.3f", psi_0))) |>
    kbl(
      caption = "MLE Table",
      align = "c",
      escape = FALSE
    ) |>
    kable_material_dark(
      lightable_options = c("hover"),
      position = "center",
      font_size = 18
    ) |>
    row_spec(nrow(MLE_data) + 1, color = "green", bold = TRUE)
  
  # Format Confidence Interval table
  ci_table <- conf_ints |>
    mutate(
      contains_truth = (lower <= psi_0 & upper >= psi_0),
      Status = ifelse(contains_truth, "✅", "❌"),
      Interval = sprintf("(%.3f, %.3f)", lower, upper),
      Length = sprintf("%.3f", upper - lower),
      `Confidence Level` = confidence
    ) |>
    rename(Source = pseudolikelihood) |>
    select(`Confidence Level`, Source, Interval, Length, Status)
  
  # Add "Truth" rows for each confidence level
  confidence_levels <- unique(ci_table$`Confidence Level`)
  truth_rows <- tibble(
    `Confidence Level` = confidence_levels,
    Source = "Truth",
    Interval = sprintf("%.3f", psi_0),
    Length = "",
    Status = ""
  )
  
  ci_table <- bind_rows(ci_table, truth_rows) |>
    arrange(`Confidence Level`, Source)
  
  ci_kable <- ci_table |>
    kbl(
      caption = "Confidence Interval Table",
      escape = FALSE,
      align = "c"
    ) |>
    kable_material_dark(
      lightable_options  = c("hover"),
      position = "center",
      font_size = 18,
      html_font = "Arial"
    ) |>
    collapse_rows(columns = 1, valign = "top") |>
    row_spec(which(ci_table$Source == "Truth"), bold = TRUE, color = "green")
  
  list(MLEs = mle_table, CIs = ci_kable)
}

get_stat_fns <- function(pseudolikelihoods, LL_df) {
  
  psi_endpoints <- LL_df |> 
    pull(psi) |> 
    (\(x) c(head(x, 1), tail(x, 1)))()
  
  pseudolikelihoods |> 
    imap(
      \(pseudolikelihood, name) {
        
        stat_function(fun = pseudolikelihood,
                      geom = "line",
                      aes(color = name),
                      linewidth = 1.5,
                      show.legend = FALSE,
                      xlim = psi_endpoints)
      }
    )
}

get_ci_colors <- function(crit_df) {
  labels <- crit_df$label
  n <- length(labels)
  
  # Get 'n' colors from the Dark2 palette
  colors <- RColorBrewer::brewer.pal(n, "Dark2")
  
  # Assign names
  names(colors) <- labels
  
  return(colors)
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

summarize_confidence_intervals <- function(ci_list, true_value) {
  
  combined <- bind_rows(ci_list, .id = "source_id")
  
  # Add a coverage indicator
  combined <- combined |>
    mutate(covers_true = lower <= true_value & upper >= true_value)
  
  # Group and summarize
  summary_df <- combined |>
    group_by(pseudolikelihood, confidence) |>
    summarise(
      length = mean(upper - lower),
      coverage_rate = mean(covers_true),
      .groups = "drop"
    )
  
  return(summary_df)
}

summarize_mle_performance <- function(mle_list, true_value) {
  
  combined <- bind_rows(mle_list, .id = "source_id")
  
  summary_df <- combined |>
    group_by(Source) |>
    summarise(
      bias = mean(MLE - true_value),
      sd = sqrt(mean((MLE - mean(MLE))^2)),
      rmse = sqrt(mean((MLE - true_value)^2)),
      .groups = "drop"
    )
  
  return(summary_df)
}



