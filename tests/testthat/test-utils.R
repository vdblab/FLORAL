test_that("encoding covariates works", {
  dat = data.frame(
    colA=1:5, 
    colB=c("a", "b", "a", "c", "b"), 
    colC=factor(c("a", "b", "a", "c", "b")))
  res = data.frame(
    colA=1:5,
    colB_a=c(1, 0, 1, 0, 0),
    colB_b=c(0, 1, 0, 0, 1),
    colB_c=c(0, 0, 0, 1, 0),
    colC_a=c(1, 0, 1, 0, 0),
    colC_b=c(0, 1, 0, 0, 1),
    colC_c=c(0, 0, 0, 1, 0)
  )
  expect_equal(res, FLORAL:::clean_covariate_columns(dat, cols = colnames(dat), drop_first = FALSE))
}
)


test_that("FLORAL() works from phy", {
  require(phyloseq)
  testphy <- new(
    "phyloseq",
    otu_table = new(
      "otu_table", 
      .Data = structure(
        c(14L, 
          150L, 260L, 434L, 79L, 287L, 361L, 8L, 42L, 501L, 89L, 8L, 29L, 
          2629L, 154L, 200L, 865L, 694L, 217L, 144L, 0L, 0L, 0L, 0L, 0L, 
          10L, 0L, 0L, 521L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 46L, 0L, 0L, 
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 436L, 0L, 0L, 0L, 0L, 
          0L, 291L, 4L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
          0L, 0L, 0L, 0L, 0L, 0L, 30L, 3L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
          0L, 0L, 37L, 65L, 59L, 0L, 106L, 388L, 61L, 26L, 0L, 0L, 0L, 
          0L, 0L, 24L, 0L, 0L, 0L), dim = c(20L, 6L),
        dimnames = list(c("ASV_1216", 
                          "ASV_12580", "ASV_12691", "ASV_135", "ASV_147", "ASV_15", "ASV_16", 
                          "ASV_184", "ASV_19", "ASV_20", "ASV_21", "ASV_2202", "ASV_239", 
                          "ASV_25", "ASV_253", "ASV_260", "ASV_29", "ASV_3005", "ASV_302", 
                          "ASV_3048"), c("1000A", "1000B", "1000C", "1000D", "1000E", "1001"
                          ))), taxa_are_rows = TRUE), 
    tax_table = new("taxonomyTable", 
                    .Data = structure(
                      c("Bacteria", "Bacteria", "Bacteria", "Bacteria", 
                        "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", 
                        "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", 
                        "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", 
                        "Bacteria", "Firmicutes", "Firmicutes", "Firmicutes", "Firmicutes", 
                        "Firmicutes", "Firmicutes", "Firmicutes", "Firmicutes", "Firmicutes", 
                        "Firmicutes", "Firmicutes", "Firmicutes", "Firmicutes", "Firmicutes", 
                        "Firmicutes", "Firmicutes", "Firmicutes", "Firmicutes", "Firmicutes", 
                        "Firmicutes", "Clostridia", "Clostridia", "Clostridia", "Clostridia", 
                        "Clostridia", "Clostridia", "Bacilli", "Clostridia", "Erysipelotrichia", 
                        "Clostridia", "Clostridia", "Clostridia", "Bacilli", "Clostridia", 
                        "Erysipelotrichia", "Clostridia", "Clostridia", "Erysipelotrichia", 
                        "Clostridia", "Clostridia", "Clostridiales", "Clostridiales", 
                        "Clostridiales", "Clostridiales", "Clostridiales", "Clostridiales", 
                        "Lactobacillales", "Clostridiales", "Erysipelotrichales", 
                        "Clostridiales", "Clostridiales", "Clostridiales", "Lactobacillales", 
                        "Clostridiales", "Erysipelotrichales", "Clostridiales", "Clostridiales", 
                        "Erysipelotrichales", "Clostridiales", "Clostridiales", "Lachnospiraceae", 
                        "Lachnospiraceae", "Lachnospiraceae", "Lachnospiraceae", 
                        "Lachnospiraceae", "Lachnospiraceae", "Streptococcaceae", 
                        "Lachnospiraceae", "Erysipelotrichaceae", "Lachnospiraceae", 
                        "Lachnospiraceae", "Ruminococcaceae", "Streptococcaceae", 
                        "Ruminococcaceae", "Erysipelotrichaceae", "Lachnospiraceae", 
                        "Lachnospiraceae", "Erysipelotrichaceae", "Lachnospiraceae", 
                        "Lachnospiraceae", "<not present>", "Blautia", "<not present>", 
                        "Blautia", "Anaerostipes", "Sellimonas", "Streptococcus", 
                        "Eisenbergiella", "[Clostridium] innocuum group", "<not present>", 
                        "Lachnoclostridium", "<not present>", "Lactococcus", "Subdoligranulum", 
                        "Candidatus Stoquefichus", "Blautia", "[Ruminococcus] gnavus group", 
                        "Faecalitalea", "Blautia", "<not present>"), dim = c(20L, 6L), 
                      dimnames = list(
                        c("ASV_1216", "ASV_12580", "ASV_12691", 
                          "ASV_135", "ASV_147", "ASV_15", "ASV_16", "ASV_184", "ASV_19", 
                          "ASV_20", "ASV_21", "ASV_2202", "ASV_239", "ASV_25", "ASV_253", 
                          "ASV_260", "ASV_29", "ASV_3005", "ASV_302", "ASV_3048"), 
                        c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"
                        )))), 
    sam_data = new(
      "sample_data",
      .Data = list(
        c("1000", 
          "1000", "1000", "1000", "1000", "pt_with_samples_1001_1002_1003_1004_1005_1006_1007_1008_1048_1121_152"
        ), 
        c(0L, 5L, 15L, 18L, 22L, 15L),
        c("formed", "liquid", "liquid", "semi-formed", "formed", "formed"),
        c("SRR11414397", "SRR11414992", "SRR11414991", "SRR11414990", "SRR11414989", "SRR11414988"), 
        c("PRJNA545312", "PRJNA545312", "PRJNA545312", "PRJNA545312", 
          "PRJNA545312", "PRJNA545312"), c(-9, -4, 6, 9, 13, -6),
        c("", "", "", "", "", "")),
      names = c("PatientID", "Timepoint", "Consistency", "Accession", "BioProject", "DayRelativeToNearestHCT", 
                "AccessionShotgun"),
      row.names = c("1000A", "1000B", "1000C", "1000D", "1000E", "1001"), .S3Class = "data.frame"), phy_tree = NULL, 
    refseq = NULL)
  
  
  dat <- FLORAL:::phy_to_floral_data(
    testphy, covariates=c("Consistency"), y = "DayRelativeToNearestHCT")
  
  expect_no_error(
    fit <- FLORAL(dat$xcount,dat$y,family="gaussian",progress=FALSE,step2=TRUE, ncv = 1)
  )

})
