#' extract sequence from peaks file
#'
#' @param peak bed format file or narrow peaks file.
#' @param genomefie the path for genome file.
#' @param outfasta the path for output fasta file.
#' @param type input data file, "bed" or "narrowpeak", default "bed".
#' @param a the upstream extend bp for peak, default 0.
#' @param b the downstream extend bp for peak, default 0.
#' @param minlen the minimum peak length to extract sequence, default 5.
#' @param pythonPath python path.
#'
#' @return fasta file
#' @export
getSeqFromBed <- function(peak = NULL,
                          genomefie = NULL,
                          outfasta = NULL,
                          type = "bed",
                          a = 0,b = 0,minlen = 5,
                          pythonPath = NULL){
  if(!ensure_python_module("pyfaidx",pythonPath)){
    return(invisible(NULL))
  }

  # run code
  pyscript.path = system.file("extdata", "getSeqFromBed.py", package = "metaplot")
  reticulate::source_python(pyscript.path)

  suppressMessages(
    reticulate::py$GetPeaksFa(peak = peak,
                              genomefie = genomefie,
                              outfasta = outfasta,
                              type = type,
                              a = as.integer(a),
                              b = as.integer(b),
                              minlen = as.integer(minlen))
  )
}
