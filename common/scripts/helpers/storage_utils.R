save_list_objects <- function(object_list, dir_path) {
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  for (object in names(object_list)) saveRDS(object_list[[object]], file = file.path(dir_path, paste0(object, ".rds")))
}

save_list_plots <- function(plots_list, dir_path) {
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  for (plot in names(plots_list)) {
    
    ggsave(filename = here(dir_path, paste0(plot, ".png")), 
           plot = plots_list[[plot]],
           width = 10,
           height = 4.5,
           dpi = 120,           
           units = "in")
  }
}

write_strict_yaml <- function(object, file) {
  yaml_str <- yaml::as.yaml(
    object,
    handlers = list(
      logical = function(x) {
        if (x) {
          structure("true", class = "verbatim")
        } else {
          structure("false", class = "verbatim")
        }
      }
    )
  )
  writeLines(yaml_str, file)
}

