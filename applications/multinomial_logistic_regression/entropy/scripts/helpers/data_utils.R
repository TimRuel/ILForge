# applications/multinomial_logistic_regression/entropy/scripts/helpers/data_utils.R

get_X1 <- function(X1_levels) {
  
  X1_level_names <- names(X1_levels)
  
  X1_ref_level <- X1_levels |> 
    map_lgl(\(x) x$ref_level) |> 
    which() |> 
    names()
  
  m <- map_dbl(X1_levels, \(x) x$m)
  
  X1_level_names |>
    rep(times = m) |>
    factor() |>
    relevel(X1_ref_level)
}

get_X2 <- function(X1_levels) {
  
  X1_levels |> 
    imap(\(level, h) {
      m <- level$m
      list2env(level$X2, environment())
      args <- as.list(c(m, dist$params))
      samples <- do.call(dist$name, args)
      support <- unlist(level$X2$support)
      X2_obs <- support[1] + diff(support) * samples
      set_names(X2_obs, rep(h, m))
    }) |> 
    unname() |>
    unlist()
}

get_X1_ref_level <- function(X1_levels) {

  X1_levels |>
    map_lgl(\(x) x$ref_level) |>
    which() |>
    names()
}

get_X1_level_of_interest <- function(X1_levels) {
  
  X1_levels |>
    map_lgl(\(x) x$level_of_interest) |>
    which() |>
    names()
}

get_X_design <- function(formula, ...) {
  
  formula <- as.formula(formula)

  df <- data.frame(...)
  mf <- model.frame(formula, data = df)
  X_design <- model.matrix(attr(mf, "terms"), data = mf)
  
  attr(X_design, "original_model_frame") <- mf
  attr(X_design, "terms") <- terms(formula, data = df)
  attr(X_design, "formula") <- formula
  attr(X_design, "contrasts") <- attr(X_design, "contrasts")
  
  return(X_design)
}

get_X_h_design <- function(X_design, X1_levels) {
  
  X1_ref_level <- get_X1_ref_level(X1_levels)
  
  h <- get_X1_level_of_interest(X1_levels)
  
  X1_level_names <- names(X1_levels)
  
  X1_main_effect_names <- paste0("X1", X1_level_names)
  
  X1X2_interaction_names <- paste0("X1", setdiff(X1_level_names, X1_ref_level)) |>
    paste0(":X2")
  
  colnames(X_design) <- c(X1_main_effect_names, "X2", X1X2_interaction_names)
  
  rows_to_keep <- X_design[, paste0("X1", h)] == 1
  
  X_h_design <- X_design[rows_to_keep,]
  
  return(X_h_design)
}

get_Y_design <- function(model_df) {
  
  model_df |>
    pull(Y) |>
    (\(Y) model.matrix(~ Y)[,-1])()
}

recover_model_frame <- function(X_design) {
  
  mf <- attr(X_design, "original_model_frame")
  
  if (is.null(mf)) {
    stop("No model frame stored in the design matrix.")
  }
  
  as.data.frame(mf)
}

get_Y_probs <- function(X_design, Beta_0) {
  
  df <- recover_model_frame(X_design)
  
  Y_probs <- X_design %*% cbind(0, Beta_0) |>
    apply(1, softmax) |>
    t() |>
    data.frame() |>
    rename_with( ~ paste0("Y", 1:(ncol(Beta_0) + 1)))
  
  cbind(df, Y_probs)
}

get_Y <- function(Y_probs) {
  
  df <- Y_probs |>
    select(starts_with("Y"))
  
  J <- ncol(df)
  
  df |>
    apply(1, \(prob) sample(1:J, size = 1, prob = prob)) |>
    unlist() |>
    unname() |>
    factor(levels = 1:J)
}

get_data <- function(X1_levels, formula, Beta_0) {
  
  while (TRUE) {
    
    X1 <- get_X1(X1_levels)
    
    X2 <- get_X2(X1_levels)
    
    X_design <- get_X_design(formula, X1, X2)
    
    Y_probs <- get_Y_probs(X_design, Beta_0)
    
    Y <- get_Y(Y_probs)
    
    if (all(table(Y) > 0)) break
  }

  model_df <- data.frame(X1 = X1,
                         X2 = X2,
                         Y)
  
  data <- list(model_df = model_df,
               X_design = X_design,
               Y_probs = Y_probs)
  
  return(data)
}
