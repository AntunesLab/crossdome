#' HLA Database
#'
#' Description
#'
#' @format A data frame with 901893 rows and 9 variables:
#' \describe{
#'   \item{peptide_sequence}{Description}
#'   \item{hla_allele}{Description}
#'   \item{peptide_length}{Description}
#'   \item{resource}{Description}
#'   \item{n_resource}{Description}
#' }
"hla_database"

#' HPA Expression
#'
#' Description
#'
#' @format A data frame with 20090 rows and 42 variables:
#' \describe{
#'   \item{ensembl_id}{Description}
#'   \item{gene_donor}{Description}
#'   \item{healthy_tissues}{Description}
#'   \item{Group}{Description}
#'   \item{spec_degree}{Description}
#'   \item{tissues}{Description}
#' }
#' @source \url{https://www.proteinatlas.org/}
"hpa_database"

#' MAGE3A Off-targets (Gee et. al. 2018)
#'
#' Description
#' Facile method for screening clinical T cell receptors for off-target peptide-HLA reactivity
#'
#' @format A data frame with 61 rows and 5 variables:
#' \describe{
#'   \item{peptide_sequence}{Description}
#'   \item{hla_allele}{Description}
#'   \item{peptide_length}{Description}
#'   \item{resource}{Description}
#'   \item{n_resource}{Description}
#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/472480v1.full}
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6002776/}
"mage_off_targets"

#' Peptide Annotation
#'
#' Description
#'
#' @format A data frame with 12158833 rows and 3 variables:
#' \describe{
#'   \item{peptide_sequence}{Description}
#'   \item{ensembl_id}{Description}
#'   \item{gene_donor}{Description}
#' }
"peptide_annotation"
