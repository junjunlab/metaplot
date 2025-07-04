#' calculate relative position for motif
#'
#' @param peaksfa path of fasta file for motif search.
#' @param motif motif sequence for searching, regular expression is accepted,
#' default `[G|A][G|A]AC[A|C|T]|[T|G|A]GT[C|T][C|T]` which means "RRACH".
#' @param mid the position for calculation the relative position, default 3.
#' @param pythonPath python path.
#'
#' @return the relative position.
#' @export
CalculatePeaksRelativeDist <- function(peaksfa = NULL,
                                       motif = "[G|A][G|A]AC[A|C|T]|[T|G|A]GT[C|T][C|T]",
                                       mid = 3,
                                       pythonPath = NULL){
  reticulate::py_config()

  # check python
  if(reticulate::py_available() == FALSE){
    message("Please install python first!")
  }else{
    if(!is.null(pythonPath)){
      reticulate::use_python(pythonPath)
    }

    # check re
    if (!reticulate::py_module_available("re")) {
      cat("Installing re ...\n")
      reticulate::py_install("re")
    }

    # run code
    pyscript.path = system.file("extdata", "CalculatePeaksRelativeDist.py", package = "metaplot")
    reticulate::source_python(pyscript.path)

    suppressMessages(
      reticulate::py$CalculatePeaksRelativeDist(peaksfa = peaksfa,
                                                motif = motif,
                                                mid = as.integer(mid))
    )
  }
}
