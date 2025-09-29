## utils_io.R
library(data.table)
library(purrr)
library(glue)
#' Create directories if they do not exist
#' 
#' @param dir_list A character vector (or list) of directory paths
#' @return None (side effect: creates directories on disk)
create_directories <- function(dir_list) {
  purrr::walk(dir_list, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))
}