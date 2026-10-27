#' calculate relative position for motif
#'
#' @param peaksfa path of fasta file for motif search.
#' @param motif motif sequence for searching, regular expression is accepted,
#' default `"[GA][GA]AC[ACT]|[TGA]GT[CT][CT]"` which means "RRACH".
#' @param mid the position for calculation the relative position, default 3.
#' @param pythonPath python path.
#'
#' @return the relative position.
#' @export
CalculatePeaksRelativeDist <- function(peaksfa = NULL,
                                       motif = "[GA][GA]AC[ACT]|[TGA]GT[CT][CT]",
                                       mid = 3,
                                       pythonPath = NULL){
  .check_file_arg(peaksfa,"peaksfa")

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
