main.data <- data.frame(
    Title = c(
        "ref_MCA",
        "ref_tabula_muris_drop",
        "ref_tabula_muris_facs",
        "ref_mouse.rnaseq",
        "ref_moca_main",
        "ref_immgen",
        "ref_hema_microarray",
        "ref_cortex_dev",
        "ref_pan_indrop",
        "ref_pan_smartseq2",
        "ref_mouse_atlas"
    ),
    Description = c(
        "Mouse Cell Atlas",
        "Tabula Muris (10X)",
        "Tabula Muris (SmartSeq2)",
        "Mouse RNA-seq from 28 cell types",
        "Mouse Organogenesis Cell Atlas (main cell types)",
        "Mouse sorted immune cells",
        "Human hematopoietic cell microarray",
        "Human cortex development scRNA-seq",
        "Human pancreatic cell scRNA-seq (inDrop)",
        "Human pancreatic cell scRNA-seq (SmartSeq2)",
        "Mouse Atlas scRNA-seq from 321 cell types"
    ),
    RDataPath = file.path("clustifyrdatahub",
                          c("ref_MCA.rda",
                            "ref_tabula_muris_drop.rda",
                            "ref_tabula_muris_facs.rda",
                            "ref_mouse.rnaseq.rda",
                            "ref_moca_main.rda",
                            "ref_immgen.rda",
                            "ref_hema_microarray.rda",
                            "ref_cortex_dev.rda",
                            "ref_pan_indrop.rda",
                            "ref_pan_smartseq2.rda",
                            "ref_mouse_atlas.rda"
                          )),
    BiocVersion="3.18",
    Genome=c(
        "mm10",
        "hg38",
    ),
    SourceType=c(
        "FASTQ",
        "RDA"
    ),
    SourceUrl=c(
        "https://ndownloader.figshare.com/files/10756795",
  ),
    SourceVersion=c(
        "7",
        "3",
        "3",
        "mouse.rnaseq.rda",
        "1",
        "immgen.rda",
        "1",
        "1",
        "1",
        "1",
        "1"
    ),
    Species=c(
        "Mus musculus",
        "Homo sapiens",
    ),
    TaxonomyId=c(
        "10090",
        "9606",
    ),
    Coordinate_1_based=FALSE,
    DataProvider=c(
        "GEO",
        "GEO",
    ),
    Maintainer="Kent Riemondy <kent.riemondy@cuanschutz.com>",
    RDataClass="data.frame",
    DispatchClass="Rda",
    stringsAsFactors = FALSE
)

write.csv(file="metadata.csv", main.data, row.names=FALSE)
