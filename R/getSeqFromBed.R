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
  .check_file_arg(peak,"peak")
  .check_file_arg(genomefie,"genomefie")

  if(is.null(outfasta) || !is.character(outfasta) || length(outfasta) != 1 ||
     is.na(outfasta) || !nzchar(outfasta)){
    stop("`outfasta` must be a single output file path.",call. = FALSE)
  }

  type <- match.arg(type,c("bed","narrowpeak"))

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

#' Check that an argument refers to an existing file
#'
#' @param path Value supplied by the caller.
#' @param arg Name of the argument, used in the error message.
#'
#' @return `TRUE` (invisibly) when the file exists, otherwise an error.
#' @keywords internal
#' @noRd
.check_file_arg <- function(path,arg){
  if(is.null(path)){
    stop("`",arg,"` must be supplied.",call. = FALSE)
  }

  if(!is.character(path) || length(path) != 1 || is.na(path)){
    stop("`",arg,"` must be a single file path.",call. = FALSE)
  }

  if(!file.exists(path)){
    stop("`",arg,"` does not exist: ",path,call. = FALSE)
  }

  invisible(TRUE)
}
