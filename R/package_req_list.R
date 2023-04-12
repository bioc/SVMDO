#' @title SVMDO
#' @name package_req_list
#' @return List of packages involved in SVMDO

if (!requireNamespace("nortest", quietly = TRUE))
  stop("Install 'nortest' to use this package.")

if (!requireNamespace("e1071", quietly = TRUE))
  stop("Install 'e1071' to use this package.")

if (!requireNamespace("BSDA", quietly = TRUE))
  stop("Install 'BSDA' to use this package.")

if (!requireNamespace("data.table", quietly = TRUE))
  stop("Install 'data.table' to use this package.")

if (!requireNamespace("sjmisc", quietly = TRUE))
  stop("Install 'sjmisc'  to use this package.")

if (!requireNamespace("klaR", quietly = TRUE))
  stop("Install 'klaR' to use this package.")

if (!requireNamespace("caTools", quietly = TRUE))
  stop("Install 'caTools' to use this package.")

if (!requireNamespace("caret", quietly = TRUE))
  stop("Install 'caret' to use this package.")

if (!requireNamespace("survival", quietly = TRUE))
  stop("Install 'survival' to use this package.")

if (!requireNamespace("DOSE", quietly = TRUE))
  stop("Install 'DOSE' to use this package.")

if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  stop("Install 'AnnotationDbi' to use this package.")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  stop("Install 'org.Hs.eg.db' to use this package.")

if (!requireNamespace("dplyr", quietly = TRUE))
  stop("Install 'dplyr' to use this package.")

if (!requireNamespace("shiny", quietly = TRUE))
  stop("Install 'shiny' to use this package.")

if (!requireNamespace("shinytitle", quietly = TRUE))
  stop("Install 'shinytitle' to use this package.")

if (!requireNamespace("shinyFiles", quietly = TRUE))
  stop("Install 'shinyFiles' to use this package.")

if (!requireNamespace("golem", quietly = TRUE))
  stop("Install 'golem' to use this package.")