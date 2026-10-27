#' Ensure a Python environment and module are available
#'
#' Shared setup used by the Python-backed wrapper functions
#' (`findMotif()`, `getSeqFromBed()`, `CalculatePeaksRelativeDist()`):
#' optionally points reticulate at a specific Python installation and installs
#' `module` if it isn't already available.
#'
#' @param module name of the Python module required by the caller.
#' @param pythonPath optional path to a Python interpreter to use.
#'
#' @return `TRUE` (invisibly) if Python and `module` are ready to use,
#' `FALSE` (invisibly) if no Python installation could be found.
#' @keywords internal
#' @noRd
ensure_python_module <- function(module,pythonPath = NULL){
  # reticulate can only be pointed at an interpreter *before* Python is
  # initialised, so this has to run before py_available() touches the runtime
  if(!is.null(pythonPath)){
    reticulate::use_python(pythonPath,required = TRUE)
  }

  # py_available(initialize = TRUE) returns FALSE instead of erroring when no
  # interpreter can be found, which keeps the message below reachable
  if(!reticulate::py_available(initialize = TRUE)){
    message("Please install python first!")
    return(invisible(FALSE))
  }

  if(!reticulate::py_module_available(module)){
    # standard-library modules always ship with Python and have no PyPI package
    if(module %in% c("re","sys","os","math","json","collections","itertools")){
      stop("Python standard library module `",module,"` is not available in the ",
           "interpreter selected by reticulate.",call. = FALSE)
    }

    message("Installing ",module," ...")
    reticulate::py_install(module)
  }

  invisible(TRUE)
}
