library(here)
rediportal <- data.frame(
    Title = c(
        "Rediportal RNA editing sites, all data, (mm10)",
        "Rediportal RNA editing sites, coordinates only, (mm10)",
        "Rediportal RNA editing sites, all data, (hg38)",
        "Rediportal RNA editing sites, coordinates only, (hg38)",
        "Mouse adenosine-to-inosine RNA recoding sites",
        "Human adenosine-to-inosine RNA recoding sites"
    ),
    Description = c(
        "GRanges object containing 107,094 mouse RNA editing sites, mm10, additional columns describe editing site characteristics, as described by REDIPORTAL http://srv00.recas.ba.infn.it/atlas/help.html",
        "GRanges object containing 107,094 mouse RNA editing sites",
        "GRanges object containing 15,638,648 human RNA editing sites, hg38, additional columns describe editing site characteristics, as described by REDIPORTAL http://srv00.recas.ba.infn.it/atlas/help.html",
        "GRanges object containing 15,638,648 human RNA editing sites, hg38",
        "GRanges object containing 1,160 mouse CDS recoding RNA editing sites, mm10",
        "GRanges object containing 1,517 human CDS recoding RNA editing sites, hg38"
    ),
    RDataPath = file.path("raerdata",
                          "1.0.0",
                          c(
                              "rediportal_full_mm10.rda",
                              "rediportal_coords_mm10.rda",
                              "rediportal_full_hg38.rda",
                              "rediportal_coords_hg38.rda",
                              "gabay_sites_mm10.rds",
                              "gabay_sites_hg38.rds"
                           )),
    BiocVersion="3.18",
    Genome=c(
        "mm10",
        "mm10",
        "hg38",
        "hg38",
        "mm10",
        "hg38"
    ),
    SourceType="RDA",
    SourceUrl=c(
        "http://srv00.recas.ba.infn.it/atlas/",
        "http://srv00.recas.ba.infn.it/atlas/",
        "http://srv00.recas.ba.infn.it/atlas/",
        "http://srv00.recas.ba.infn.it/atlas/",
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28841-4/MediaObjects/41467_2022_28841_MOESM3_ESM.xlsx",
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28841-4/MediaObjects/41467_2022_28841_MOESM3_ESM.xlsx"
  ),
    SourceVersion="1",
    Species=c(
        "Mus musculus",
        "Mus musculus",
        "Homo sapiens",
        "Homo sapiens",
        "Mus musculus",
        "Homo sapiens"
    ),
    TaxonomyId=c(
        "10090",
        "10090",
        "9606",
        "9606",
        "10090",
        "9606"
    ),
    Coordinate_1_based = TRUE,
    DataProvider=c(
        "REDIportal",
        "REDIportal",
        "REDIportal",
        "REDIportal",
        "https://doi.org/10.1038/s41467-022-28841-4",
        "https://doi.org/10.1038/s41467-022-28841-4"
    ),
    Maintainer="Kent Riemondy <kent.riemondy@cuanschutz.com>",
    RDataClass="GRanges",
    DispatchClass="Rda",
    stringsAsFactors = FALSE
)

adar_ko <- data.frame(

)

sc_data <- data.frame(

)

rna_dna <- data.frame(

)

misaligned_sites <- data.frame(

)

main.data <- rbind(rediportal, adar_ko, pbmc_sc, rna_dna, misaligned_sites)
write.csv(file=here("inst/extdata/metadata.csv"), main.data, row.names=FALSE)
