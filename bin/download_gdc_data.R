library(TCGAbiolinks)

source("global.R")

main = function() {
  
  for (tumor_type in TUMOR_TYPES) {
    download_TCGA_project(tumor_type)
  }
  
}

download_TCGA_project = function(tumor_type) {
  
  project = paste0('TCGA-', tumor_type)
  
  query = GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal"),
    access = "open"
  )
  GDCdownload(query, directory = GDC_DATA)
  
  query <- GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Gene Level Copy Number",
    workflow.type = "ASCAT3",
    access="open"
  )
  GDCdownload(query, directory = GDC_DATA)

}

main()