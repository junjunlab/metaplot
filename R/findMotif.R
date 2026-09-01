#' find motif with given fasta files
#'
#' @param peaksfa path of fasta file for motif search
#' @param motif motif sequence for searching, regular expression is accepted,
#' default `"[G|A][G|A]AC[A|C|T]|[T|G|A]GT[C|T][C|T]"` which means "RRACH".
#' @param pythonPath python path.
#'
#' @return target motif
#' @export
findMotif <- function(peaksfa = NULL,
                      motif = "[G|A][G|A]AC[A|C|T]|[T|G|A]GT[C|T][C|T]",
                      pythonPath = NULL){
  if(!ensure_python_module("re",pythonPath)){
    return(invisible(NULL))
  }

  # run code
  pyscript.path = system.file("extdata", "findMotif.py", package = "metaplot")
  reticulate::source_python(pyscript.path)

  suppressMessages(
    reticulate::py$findMotif(peaksfa = peaksfa,motif = motif)
  )
}
