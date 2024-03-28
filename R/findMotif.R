#' find motif with given fasta files
#'
#' @param peaksfa path of fasta file for motif search
#' @param motif motif sequence for searching, regular expression is accepted,
#' default "[G|A][G|A]AC[A|C|T]|[T|G|A]GT[C|T][C|T]" which means "RRACH".
#' @param pythonPath python path.
#'
#' @return target motif
#' @export
findMotif <- function(peaksfa = NULL,
                      motif = "[G|A][G|A]AC[A|C|T]|[T|G|A]GT[C|T][C|T]",
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
    pyscript.path = system.file("extdata", "findMotif.py", package = "metaplot")
    reticulate::source_python(pyscript.path)

    suppressMessages(
      reticulate::py$findMotif(peaksfa = peaksfa,motif = motif)
    )
  }
}
