#' Ensure a Python environment and module are available
#'
#' Shared setup used by the Python-backed wrapper functions
#' (`findMotif()`, `getSeqFromBed()`, `CalculatePeaksRelativeDist()`):
#' configures reticulate, optionally points at a specific Python
#' installation, and installs `module` if it isn't already available.
#'
#' @param module name of the Python module required by the caller.
#' @param pythonPath optional path to a Python interpreter to use.
#'
#' @return `TRUE` (invisibly) if Python and `module` are ready to use,
#' `FALSE` (invisibly) if no Python installation could be found.
#' @keywords internal
#' @noRd
ensure_python_module <- function(module,pythonPath = NULL){
  reticulate::py_config()

  if(!reticulate::py_available()){
    message("Please install python first!")
    return(invisible(FALSE))
  }

  if(!is.null(pythonPath)){
    reticulate::use_python(pythonPath)
  }

  if(!reticulate::py_module_available(module)){
    message("Installing ",module," ...")
    reticulate::py_install(module)
  }

  invisible(TRUE)
}
