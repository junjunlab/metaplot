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
  if(!ensure_python_module("re",pythonPath)){
    return(invisible(NULL))
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
