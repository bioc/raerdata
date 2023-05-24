#' raerdata
#'
#' @description A collection of datasets and databases to demonstrate identification of RNA-editing sites using the raer package.
#'
#' @docType package
#' @keywords internal
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom utils read.csv
#' @rdname raerdata
#' @name raerdata
NULL

#' Databases of known RNA editing sites
#'
#'
#' @details
#' `rediportal_full_hg38()` will download the human rediportal database for hg38 which has been
#'  coverted into a GRanges object. The Granges is supplemented with additional columns of information
#'  provided by the rediportal database, including gene location, repeat type, dbSNP annotation, and
#'  potential for amino-acid recoding.
#'
#' `rediportal_coords_hg38()` will download the human rediportal database for hg38 which has been
#'  coverted into a GRanges object, which only contains the coordinates of the editing site.
#'
#' `rediportal_full_mm10()` will download the mouse rediportal database for mm10 which has been
#'  coverted into a GRanges object. The Granges is supplemented with additional columns of information
#'  provided by the rediportal database, including gene location, repeat type, dbSNP annotation, and
#'  potential for amino-acid recoding.
#'
#' `rediportal_coords_mm10()`will download the mouse rediportal database for mm10 which has been
#'  coverted into a GRanges object, which only contains the coordinates of the editing site.
#'
#'  `gabay_sites_hg38()` will download high-confidence human CDS editing sites (hg38).
#'
#'  `gabay_sites_mm10()` will download high-confidence mouse CDS editing sites (lifted-over from hg38 to mm10).
#'
#' @rdname atlases
#' @examples
#' rediportal_full_mm10()
NULL

#' @rdname rediportal
"rediportal_full_mm10"
#' @rdname rediportal
"gabay_sites_hg38"
#' @rdname rediportal
"rediportal_coords_mm10"


#' Whole genome and RNA sequencing data from NA12877 cell line
#'
#'
#' @details
#' `NA12877_wgs()` will download a BAM and BAM index file from whole genome sequencing
#'  of the NA12877 cell line.
#'
#' `NA12877_rnaseq()` will download a BAM and BAM index file from RNA-seq
#'  of the NA12877 cell line.
#' @returns A data.frame containing:
#' \itemize{
#' \item \code{Name}, the name of the BAM file without the \code{*.bam} extension.
#' \item \code{Description}, a short string containing the experimental condition and replicate number.
#' \item \code{Path}, a \linkS4class{BamFile} containing the path to each BAM file.
#'
#' @rdname NA12877
#'
NULL

#' RNA sequencing data from study GSE99249
#'
#' @description
#' Study [GSE99249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99249) examined RNA editing in
#' WT, ADAR1KO, and ADAR1-p150 HEK293T cells treated with and without Interferon Beta.
#'
#' @details
#' `GSE99249()` will download a BAM and BAM index files from 6 RNA-seq libraries
#'
#' @returns A data.frame containing:
#' \itemize{
#' \item \code{Name}, the name of the BAM file without the \code{*.bam} extension.
#' \item \code{Description}, a short string containing the experimental condition and replicate number.
#' \item \code{Path}, a \linkS4class{BamFile} containing the path to each BAM file.
#'
#' @rdname GSE99249
#'
NULL

#' RNA sequencing data from study GSE99249
#'
#' @description
#' Study [GSE99249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99249) examined RNA editing in
#' WT, ADAR1KO, and ADAR1-p150 HEK293T cells treated with and without Interferon Beta.
#'
#' @details
#' `GSE99249()` will download a BAM and BAM index files from 6 RNA-seq libraries
#'
#' @returns A data.frame containing:
#' \itemize{
#' \item \code{Name}, the name of the BAM file without the \code{*.bam} extension.
#' \item \code{Description}, a short string containing the experimental condition and replicate number.
#' \item \code{Path}, a \linkS4class{BamFile} containing the path to each BAM file.
#'
#' @rdname GSE99249
#'
NULL

#' RNA sequencing data from study GSE99249
#'
#' @description
#' Study [GSE99249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99249) examined RNA editing in
#' WT, ADAR1KO, and ADAR1-p150 HEK293T cells treated with and without Interferon Beta.
#'
#' @details
#' `GSE99249()` will download a BAM and BAM index files from 6 RNA-seq libraries
#'
#' @returns A data.frame containing:
#' \itemize{
#' \item \code{Name}, the name of the BAM file without the \code{*.bam} extension.
#' \item \code{Description}, a short string containing the experimental condition and replicate number.
#' \item \code{Path}, a \linkS4class{BamFile} containing the path to each BAM file.
#'
#' @rdname GSE99249
#'
NULL
