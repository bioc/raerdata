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
                          "rediportal",
                          "1.0.0",
                          c(
                              "rediportal_full_mm10.rda",
                              "rediportal_coords_mm10.rda",
                              "rediportal_full_hg38.rda",
                              "rediportal_coords_hg38.rda",
                              "gabay_sites_mm10.rda",
                              "gabay_sites_hg38.rda"
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
    SourceType="TSV",
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
    Maintainer="kent.riemondy@cuanschutz.com",
    RDataClass="GRanges",
    DispatchClass="Rda"
)

sc_data <- data.frame(
    Title = c(
        "10k PBMC scRNA-seq library BAM",
        "10k PBMC scRNA-seq library BAM index",
        "10k PBMC scRNA-seq library SingleCellExperiment",
        "Rediportal RNA editing sites, chr16, coordinates only, (hg38)"
    ),
    Description = c(
        "BAM alignments on chr16 from 10x Genomics 3' scRNA-seq library from 10k PBMCs",
        "BAM index of alignments on chr16 from 10x Genomics 3' scRNA-seq library from 10k PBMCs",
        "Preprocessed SingleCellExperiment containing UMI counts from 10x Genomics 3' scRNA-seq library from 10k PBMCs",
        "GRanges object containing human RNA editing sites from chr16, hg38"
    ),
    RDataPath = file.path("raerdata",
                          "10x",
                          "1.0.0",
                          c(
                              "10k_PBMC_3p_nextgem_Chromium_X_intron_possorted_chr16_rp.bam",
                              "10k_PBMC_3p_nextgem_Chromium_X_intron_possorted_chr16_rp.bam.bai",
                              "human_pbmc_sce.rda",
                              "rediportal_chr16.bed.gz"
                          )
    ),
    BiocVersion="3.18",
    Genome="hg38",
    SourceType=c(
        "BAM",
        "BAM",
        "MTX",
        "TSV"
    ),
    SourceUrl= c(    
        "https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-x-with-intronic-reads-3-1-high",
        "https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-x-with-intronic-reads-3-1-high",
        "https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-x-with-intronic-reads-3-1-high",
        "http://srv00.recas.ba.infn.it/atlas/"
    ),
    SourceVersion="1",
    Species="Homo sapiens",
    TaxonomyId="9606",
    Coordinate_1_based = c(
        FALSE,
        FALSE,
        NA,
        TRUE
    ),
    DataProvider=c(
        "10x Genomics",
        "10x Genomics",
        "10x Genomics",
        "REDIportal"
    ),
    Maintainer="kent.riemondy@cuanschutz.com",
    RDataClass=c(
        "character",
        "character",
        "SingleCellExperiment",
        "GRanges"
    ),
    DispatchClass=c(
        "FilePath",
        "FilePath",
        "Rda",
        "BEDFile"
    )
)

rna_dna <- data.frame(
    Title = c(
        "WGS of NA12878 cell line, first megabase of chr4, (BAM)",
        "WGS of NA12878 cell line, first megabase of chr4, (BAI)",
        "RNA-seq of NA12878 cell line, first megabase of chr4, (BAM)",
        "RNA-seq of NA12878 cell line, first megabase of chr4, (BAI)",
        "SNPs from chr4 in dbSNP155",
        "Genome sequence of first megabase of chr4 (hg38)"
    ),
    Description = c(
        "Whole genome sequencing of NA12878 cell line, alignments from first megabase of chr4, (BAM)",
        "Whole genome sequencing of NA12878 cell line, alignments from first megabase of chr4, (BAI)",
        "RNA sequencing of NA12878 cell line, alignments from first megabase of chr4, (BAM)",
        "RNA sequencing of NA12878 cell line, alignments from first megabase of chr4, (BAI)",
        "SNPs from dbSNP155 within first megabase of chr4 (hg38)",
        "Genome sequence of first megabase of chr4 (hg38)"
    ),
    RDataPath = file.path("raerdata",
                          "NA12878",
                          "1.0.0",
                          c(
                              "NA12878.wgs.sub.bam",
                              "NA12878.wgs.sub.bam.bai",
                              "NA12878.rnaseq.sub.bam",
                              "NA12878.rnaseq.sub.bam.bai",
                              "chr4snps.bed.gz",
                              "hg38_chr4.fa.bgz"
                          )
    ),
    BiocVersion="3.18",
    Genome="hg38",
    SourceType=c(
        "FASTQ",
        "FASTQ",
        "FASTQ",
        "FASTQ",
        "RDS",
        "FASTA"
    ),
    SourceUrl=c(
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR262/ERR262997",
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR262/ERR262997",
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/008/SRR1258218",
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/008/SRR1258218",
        "https://bioconductor.org/packages/3.17/data/annotation/html/SNPlocs.Hsapiens.dbSNP155.GRCh38.html",
        "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz"
    ),
    SourceVersion="1",
    Species="Homo sapiens",
    TaxonomyId="9606",
    Coordinate_1_based = c(
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        NA
    ),
    DataProvider=c(
        "ENA",
        "ENA",
        "ENA",
        "ENA",
        "Bioconductor",
        "Gencode"
    ),
    Maintainer="kent.riemondy@cuanschutz.com",
    RDataClass=c(
        "character",
        "character",
        "character",
        "character",
        "GRanges",
        "FaFile"
    ),
    DispatchClass=c(
        "FilePath",
        "FilePath",
        "FilePath",
        "FilePath",
        "BEDFile",
        "FaFile"
    )
)

adar_ko <- data.frame(
    Title = c(
        "SRR5564260_chr18.bam",
        "SRR5564260_chr18.bam.bai",
        "SRR5564261_chr18.bam",
        "SRR5564261_chr18.bam.bai",
        "SRR5564269_chr18.bam",
        "SRR5564269_chr18.bam.bai",
        "SRR5564270_chr18.bam",
        "SRR5564270_chr18.bam.bai",
        "SRR5564271_chr18.bam",
        "SRR5564271_chr18.bam.bai",
        "SRR5564277_chr18.bam",
        "SRR5564277_chr18.bam.bai",
        "Genome sequence of chr18 (hg38)"
    ),
    Description = c(
            "ADAR1KO Interferon beta biological replicate 2, BAM",
            "ADAR1KO Interferon beta biological replicate 2, BAI",
            "ADAR1KO Interferon beta biological replicate 3, BAM",
            "ADAR1KO Interferon beta biological replicate 3, BAI",
            "ADAR1KO Interferon beta biological replicate 1, BAM",
            "ADAR1KO Interferon beta biological replicate 1, BAI",
            "Wildtype Interferon beta biological replicate 2, BAM",
            "Wildtype Interferon beta biological replicate 2, BAI",
            "Wildtype Interferon beta biological replicate 3, BAM",
            "Wildtype Interferon beta biological replicate 3, BAI",
            "Wildtype Interferon beta biological replicate 1, BAM",
            "Wildtype Interferon beta biological replicate 1, BAI",
            "Genome sequence of chr18 (hg38)"
    ),
    RDataPath = file.path("raerdata",
                          "GSE99249",
                          "1.0.0",
                          c("SRR5564260_chr18.bam",
                            "SRR5564260_chr18.bam.bai",
                            "SRR5564261_chr18.bam",
                            "SRR5564261_chr18.bam.bai",
                            "SRR5564269_chr18.bam",
                            "SRR5564269_chr18.bam.bai",
                            "SRR5564270_chr18.bam",
                            "SRR5564270_chr18.bam.bai",
                            "SRR5564271_chr18.bam",
                            "SRR5564271_chr18.bam.bai",
                            "SRR5564277_chr18.bam",
                            "SRR5564277_chr18.bam.bai",
                            "chr18.fasta.bgz")
    ),
    BiocVersion="3.18",
    Genome="hg38",
    SourceType=c(rep("FASTQ", 12), "FASTA"),
    SourceUrl=c(
        rep("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99249", 12),
        "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz"
    ),
    SourceVersion="1",
    Species="Homo sapiens",
    TaxonomyId="9606",
    Coordinate_1_based = c(rep(FALSE, 12), NA),
    DataProvider=c(rep("GEO", 12), "Gencode"),
    Maintainer="kent.riemondy@cuanschutz.com",
    RDataClass=c(rep("character", 12),"FaFile"),
    DispatchClass=c(rep("FilePath", 12), "FaFile")
)

main.data <- rbind(rediportal, adar_ko, sc_data, rna_dna)
write.csv(file=here("inst/extdata/metadata.csv"), main.data, row.names=FALSE)
