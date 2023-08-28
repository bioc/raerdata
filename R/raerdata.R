#' raerdata
#'
#' @description A collection of datasets and databases to demonstrate RNA-editing 
#' analysis approaches using the raer package.
#' 
#' @details
#' [`atlases`][rediportal_full_hg38] a collection of RNA editing databases  
#' 
#' [`NA12878`] Whole genome and RNA sequencing from the NA12878 cell line
#' 
#' [`GSE99249`]  RNA sequencing data from a study that examined RNA editing in
#' WT, ADAR1KO, and ADAR1-p150 HEK293T cells treated with and without interferon beta.
#' 
#' [`pbmc_10x`]  single cell RNA sequencing data from human PBMCs from 10x Genomics
#' 
#' 
#' @docType package
#' @importFrom ExperimentHub ExperimentHub
#' @import Rsamtools 
#' @rdname raerdata
#' @name raerdata
NULL

#' Databases of known RNA editing sites
#'
#'
#' @details
#' `rediportal_full_hg38()` will download the human REDIportal database for hg38 which has been
#'  converted into a GRanges object. The Granges is supplemented with additional columns of information
#'  provided by the REDIportal database, including gene location, repeat type, dbSNP annotation, and
#'  potential for amino-acid recoding.
#'
#' `rediportal_coords_hg38()` will download the human REDIportal database for hg38 which has been
#'  converted into a GRanges object, which only contains the coordinates of the editing site.
#'
#' `rediportal_full_mm10()` will download the mouse REDIportal database for mm10 which has been
#'  converted into a GRanges object. The Granges is supplemented with additional columns of information
#'  provided by the REDIportal database, including gene location, repeat type, dbSNP annotation, and
#'  potential for amino-acid recoding.
#'
#' `rediportal_coords_mm10()`will download the mouse REDIportal database for mm10 which has been
#'  converted into a GRanges object, which only contains the coordinates of the editing site.
#'
#'  `gabay_sites_hg38()` will download high-confidence human CDS editing sites (hg38).
#'
#'  `gabay_sites_mm10()` will download high-confidence mouse CDS editing sites (lifted-over from hg38 to mm10).
#'
#' @rdname atlases
#' @family atlases
#' @export
rediportal_full_mm10 <- function(){
    eh <- ExperimentHub()
    eh[["EH8233"]]
}

#' @rdname atlases
#' @family atlases
#' @export
rediportal_coords_mm10 <- function(){
    eh <- ExperimentHub()
    eh[["EH8234"]]
}

#' @rdname atlases
#' @family atlases
#' @export
rediportal_full_hg38 <- function(){
    eh <- ExperimentHub()
    eh[["EH8235"]]
}

#' @rdname atlases
#' @family atlases
#' @export
rediportal_coords_hg38 <- function(){
    eh <- ExperimentHub()
    eh[["EH8236"]]
}

#' @rdname atlases
#' @export
gabay_sites_mm10 <- function(){
    eh <- ExperimentHub()
    eh[["EH8237"]]
}

#' @rdname atlases
#' @export
gabay_sites_hg38 <- function(){
    eh <- ExperimentHub()
    eh[["EH8238"]]
}

#' Whole genome and RNA sequencing data from NA12878 cell line
#'
#'
#' @details
#' Will download BAM and BAM index files from whole genome and RNA sequencing of
#' the NA12878 cell line, The data is from the first megabase of chromosome 4. Additionally 
#' a fasta file and a database of known SNPS will be downloaded. 
#' 
#' @returns A list containing:
#' \itemize{
#' \item \code{bams} A [BamFileList] object, indicating the bam file paths and BAI indexes. 
#' \item \code{fasta} the path to a fasta file containing the first megabase of chr4
#' \item \code{snps} a GRanges object containing snps from the first megabase of chr4
#' }
#' @rdname NA12878
#' @import Rsamtools 
#' @importFrom BiocGenerics path
#' @import GenomicRanges
#' @export
NA12878 <- function() {
    eh <- ExperimentHub::ExperimentHub()
    
    # bam files
    eh_ids <- c("EH8256", "EH8257", "EH8258", "EH8259")
    paths <- unlist(lapply(eh_ids, function(x) eh[[x]]))
    nms <- c("NA12878_RNASEQ", "NA12878_WGS")
    desc <- eh[eh_ids]$title
    is_bam <- grepl("(BAM)", desc)
    desc <- desc[is_bam]
    bams <- paths[is_bam]
    idxs <- paths[!is_bam]
    bams <- BamFileList(bams, idxs)
    names(bams) <- nms

    # snps
    snps <- eh[["EH8260"]]
    
    # small genome seq
    fasta <- path(eh[["EH8261"]])
    list(
        bams = bams,
        fasta = fasta,
        snps = snps
    )
}

#' RNA sequencing data from study GSE99249
#'
#' @description
#' Study [GSE99249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99249) examined RNA editing in
#' WT, ADAR1KO, and ADAR1-p150 HEK293T cells treated with and without interferon beta.
#'
#' @details
#' `GSE99249()` will download BAM and BAM index files from 6 RNA-seq libraries. 
#' 3 libraries are ADAR1 knockout cells treated with interferon beta. and 3 libraries
#' are wild type cells treated with interferon beta. The BAM files contain alignments
#' from chromosome 18. 
#'
#' @returns A list containing:
#' \itemize{
#' \item \code{bams} A [BamFileList] object, indicating the bam file paths and BAI indexes. 
#' \item \code{fasta} the path to a fasta file from chr18 of hg38
#' \item \code{snps} a [GRanges] object containing known snps from the Rediportal database (hg38)
#' }
#' @rdname GSE99249
#' @import rtracklayer
#' @export
GSE99249 <- function() {
    eh <- ExperimentHub::ExperimentHub()
    
    # get bam files
    eh_ids <- paste0("EH", 8239:8250)
    paths <- unlist(lapply(eh_ids, function(x) eh[[x]]))
    nms <- c("NA12877_WGS", "NA12877_RNASEQ")
    desc <- eh[eh_ids]$title
    is_bam <- grepl(".bam$", desc)
    desc <- desc[is_bam]
    bams <- paths[is_bam]
    idxs <- paths[!is_bam]
    bams <- BamFileList(bams, idxs)
    
    names(bams) <- sub("_chr18.bam", "", desc)
    
    # rediportal sites
    sites <- eh[["EH8236"]]
    
    # small genome seq
    fasta <- path(eh[["EH8251"]])
    list(
        bams = bams,
        fasta = fasta,
        sites = sites
    )
}

#' single cell RNA sequencing data from human PBMCs
#'
#' @description
#' A 10x Genomics 3' single cell RNA-seq library from 10k PBMCs. The BAM file contains
#' alignments from chr16. A [SingleCellExperiment] is also provided with pre-processed
#' gene expression data, a UMAP projection and cell type annotations. The dataset was obtained
#' from [10x Genomics](https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-x-with-intronic-reads-3-1-high)
#'
#' @details
#' `pbmc_10x()` will download a BAM, BAM index file, REDIportal RNA editing sites, and 
#' a SingleCellExperiment object from the [ExperimentHub]. 
#'
#' @returns A list containing:
#' \itemize{
#' \item \code{bam} a [BamFile] object indicating the BAM and BAI file paths. 
#' Contains alignments from only chr16 (hg38).
#' \item \code{sites} a [GRanges] object containing known RNA editing sites from 
#' the REDIportal database (hg38).
#' \item \code{sce} a [SingleCellExperiment] object containing gene expression data,
#'  a UMAP projection and cell type annotations. 
#' }
#' @rdname pbmc_10x
#' @import SingleCellExperiment
#' @export
pbmc_10x <- function(){
    eh <- ExperimentHub::ExperimentHub()
    bam <- unname(eh[["EH8252"]])
    bai <- eh[["EH8253"]]
    
    sites <- eh[["EH8236"]]
    sce <- eh[["EH8254"]]
    list(bam = BamFile(bam, bai),
         sites = sites,
         sce = sce)
}

