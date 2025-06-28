`%||%` <- function(a, b) if (!is.null(a)) a else b

get_seed_for_iter <- function(base_seed, iter_id) {
  
  return(as.integer(abs(Reduce(`+`, utf8ToInt(iter_id))) + base_seed))
}
