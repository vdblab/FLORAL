if (interactive()){
  cols = c("SampleID", "PatientID", "Timepoint", "Consistency", "Accession", 
           "BioProject", "DayRelativeToNearestHCT")
  #this file has duplicate rows, and has multiple rows per pool
  samples <- read.csv("https://figshare.com/ndownloader/files/33076496")[, cols] 
  samples <- samples[!duplicated(samples),]
  samples <- samples[1:100,] # Using the first 100 samples only.
  counts <- read.csv("https://figshare.com/ndownloader/files/26393788") 
  counts <- counts[counts$SampleID %in% samples$SampleID, ]
  taxonomy <- read.csv("https://figshare.com/ndownloader/files/26770997")
  
  usethis::use_data(samples,counts,taxonomy, internal = TRUE)
}