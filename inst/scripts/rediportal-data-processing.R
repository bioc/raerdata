library(data.table)
library(here)
library(tools)
library(Rsamtools)
library(GenomicRanges)
library(readxl)
library(AnnotationHub)
library(Biostrings)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)

# Known editing sites from the Rediportal database:
# Mansi L, Tangaro MA, Lo Giudice C, et al. REDIportal: millions of novel A-to-I RNA editing events from thousands of RNAseq experiments. Nucleic Acids Res [Internet]. 2021;49:D1012–D1019. Available from: http://dx.doi.org/10.1093/nar/gkaa916.

data_dir <- here("inst/scripts/rediportal")
dir.create(data_dir, showWarnings = FALSE)

# consider using fasta files directly instead of BSgenome objects
# the getSeq method is unworkably slow on our linux server for some reason but fast when using a fasta file as input
# extracting from bsgenome is very fast (a few seconds for 15 million sites) locally on macOS, so might be a drive networking issue?

# hg38 database of known RNA editing sites from Rediportal (v2)
url <- "http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz"

rediportal_flat_file <- file.path(data_dir, "TABLE1_hg38.txt.gz")
if(!file.exists(rediportal_flat_file)){
    download.file(url, rediportal_flat_file)
    md5 <- md5sum(rediportal_flat_file)
    stopifnot(unname(md5) == "cd809b839ccb38d1b6b60184a2d190e2")
}

dat <- fread(rediportal_flat_file, data.table = FALSE)
dat <- dat[order(dat$Region, dat$Position), ]

gr <- makeGRangesFromDataFrame(dat,
                               seqnames.field = "Region",
                               start.field = "Position",
                               end.field = "Position",
                               strand.field = "Strand",
                               keep.extra.columns = TRUE)
genome(gr) <- "hg38"

# double check coordinates are correct (should be all A bases w.r.t strand)
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)

# 43,223 sites are not A in reference, not sure why.
rediportal_full_hg38 <- gr[seqs == "A"]
rediportal_coords_hg38 <- rediportal_full_hg38
mcols(rediportal_coords_hg38) <- NULL

save(rediportal_full_hg38, file = file.path(data_dir, "rediportal_full_hg38.rda"))
save(rediportal_coords_hg38, file = file.path(data_dir, "rediportal_coords_hg38.rda"))


#mm10 database of known RNA editing sites from Rediportal (v2)
url <-  "http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload//TABLE1_mm10.txt.gz"

rediportal_flat_file <- file.path(data_dir, "TABLE1_mm10.txt.gz")
if(!file.exists(rediportal_flat_file)){
    download.file(url, rediportal_flat_file)
    md5 <- md5sum(rediportal_flat_file)
    stopifnot(unname(md5) == "4f0e8a3b47cee70f7bab855236e4c4e4")
}

dat <- fread(rediportal_flat_file, data.table = FALSE)
dat <- dat[order(dat$Region, dat$Position), ]
gr <- makeGRangesFromDataFrame(dat,
                               seqnames.field = "Region",
                               start.field = "Position",
                               end.field = "Position",
                               strand.field = "Strand",
                               keep.extra.columns = TRUE)
genome(gr) <- "mm10"

# double check coordinates are correct (should be all A bases w.r.t strand)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)

# all sites are A in reference, unnecessary, but may be needed in future versions
rediportal_full_mm10 <- gr[seqs == "A"]
rediportal_coords_mm10 <- rediportal_full_mm10
mcols(rediportal_coords_mm10) <- NULL

save(rediportal_full_mm10, file = file.path(data_dir, "rediportal_full_mm10.rda"))
save(rediportal_coords_mm10, file = file.path(data_dir, "rediportal_coords_mm10.rda"))

# Sites from:
# Gabay O, Shoshan Y, Kopel E, et al. Landscape of adenosine-to-inosine RNA recoding across human tissues. Nat Commun [Internet]. 2022;13:1184. Available from: http://dx.doi.org/10.1038/s41467-022-28841-4.

if(!file.exists(file.path(data_dir, "gabay_et_al_2022_sup_files.xlsx"))){
    download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28841-4/MediaObjects/41467_2022_28841_MOESM3_ESM.xlsx",
                  file.path(data_dir, "gabay_et_al_2022_sup_files.xlsx"))
}

gabay_sites <- read_excel(file.path(data_dir, "gabay_et_al_2022_sup_files.xlsx"),
                          sheet = 3)
gabay_sites_hg38 <- makeGRangesFromDataFrame(gabay_sites,
                             seqnames = "Chromosome",
                             start.field = "Postion",
                             end.field = "Postion",
                             strand.field = "Strand",
                             keep.extra.columns = TRUE)
genome(gabay_sites_hg38) <- "hg38"
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gabay_sites_hg38)
# all sites are A in hg38
stopifnot(all(seqs == "A"))
save(gabay_sites_hg38, file = file.path(data_dir, "gabay_sites_hg38.rda"))

# liftover to mm10
ah <- AnnotationHub()
chainfile <- ah[["AH14109"]]
gabay_sites_mm10 <- liftOver(gabay_sites_hg38, chainfile)
# 1448 of 1517 hg38 sites have 1-to-1 mapping,
gabay_sites_mm10 <- unlist(gabay_sites_mm10[elementNROWS(gabay_sites_mm10) == 1])

# ensure that site is an A in mm10
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gabay_sites_mm10)
# 1160 of 1148 are A
gabay_sites_mm10 <- gabay_sites_mm10[seqs == "A"]
genome(gabay_sites_mm10) <- "mm10"

save(gabay_sites_mm10, file = file.path(data_dir, "gabay_sites_mm10.rda"))

## clean up
unlink(c(file.path(data_dir, "TABLE1_hg38.txt.gz"),
         file.path(data_dir, "TABLE1_mm10.txt.gz"),
         file.path(data_dir, "gabay_et_al_2022_sup_files.xlsx")))
