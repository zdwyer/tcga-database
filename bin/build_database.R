library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(duckdb)
library(tidyr)

source("global.R")

main = function() {
  
  prepare_directory_structure()
  
  con = dbConnect(duckdb::duckdb(), dbdir = file.path(DB_PATH, DB_NAME_TEMP))
  prepare_database(con)
  gene_id_key = preprocess_gene_metadata("ACC", con)
  
  for (tumor_type in TUMOR_TYPES) {
    preprocess_tumor_type(tumor_type)
  }

  barcode_key = build_sample_metadata(con)
  build_temp_gene_data(con, barcode_key, gene_id_key)
  sort_gene_data(con)
  cleanup_database(con)
  dbDisconnect(con, shutdown = TRUE)
  
  con = dbConnect(duckdb::duckdb(), dbdir = file.path(DB_PATH, DB_NAME))
  DBI::dbExecute(con, "IMPORT DATABASE 'temp_db'")
  dbDisconnect(con, shutdown = TRUE)
  
  cleanup_directory()
}

create_temp_directory = function(dir) {
  directory = file.path("temp", dir)
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

prepare_directory_structure = function() {
  
  temp_dirs = c('sample_metadata', 'expected_count', 'tpm', 'copy_number')
  for (dir in temp_dirs) {
    create_temp_directory(dir)
  }
  
  if (!dir.exists(DB_PATH)) {
    dir.create(DB_PATH, recursive = TRUE)
  }
}

prepare_database = function(con) {
  
  dbExecute(con, "
    CREATE TABLE gene_type_key (
      gene_type_id USMALLINT,
      gene_type VARCHAR
    );
  ")
  
  dbExecute(con, "
    CREATE TABLE gene_metadata (
      gene_id USMALLINT,
      ensembl_id VARCHAR,
      gene_name VARCHAR,
      gene_type_id USMALLINT
    );
  ")

  dbExecute(con, "
    CREATE TABLE sample_metadata (
      tumor_type VARCHAR,
      barcode_id USMALLINT,
      barcode VARCHAR,
      sample VARCHAR,
      patient VARCHAR,
      tcga_sample_id VARCHAR,
      sample_type VARCHAR,
      overall_survival INTEGER,
      deceased BOOLEAN,
      preservation_method VARCHAR,
      age_at_diagnosis INTEGER,
      site_of_resection_or_biopsy VARCHAR,
      race VARCHAR,
      gender VARCHAR,
      ethnicity VARCHAR,
    )
  ")
  
  dbExecute(con, "
    CREATE TABLE gene_data_temp (
      gene_id USMALLINT,
      barcode_id USMALLINT,
      expected_count INTEGER,
      tpm FLOAT,
      copy_number UTINYINT
    );
  ")
  
  dbExecute(con, "
    CREATE TABLE gene_data (
      gene_id USMALLINT,
      barcode_id USMALLINT,
      expected_count INTEGER,
      tpm FLOAT,
      copy_number UTINYINT
    );
  ")
}

preprocess_gene_metadata = function(tumor_type, con) {
  
  project = paste0('TCGA-', tumor_type)
  
  query = GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    sample.type = c("Primary Tumor", "Solid Tissue Normal"),
    access = "open"
  )
  
  data = GDCprepare(query, directory = GDC_DATA, summarizedExperiment = TRUE)
  
  gene_metadata = as.data.frame(rowData(data)) %>% 
    tibble::rownames_to_column("rownames") %>% 
    select(ensembl_id=gene_id, gene_type, gene_name) %>%
    arrange(ensembl_id) %>%
    mutate(gene_id = row_number()) %>%
    select(gene_id, ensembl_id, gene_name, gene_type)
  
  gene_type = data.frame(gene_type = unique(gene_metadata$gene_type)) %>%
    mutate(gene_type_id = row_number())
  
  dbWriteTable(con, "gene_type_key", gene_type, append = TRUE, overwrite = FALSE)
  
  gene_metadata_final = gene_metadata %>%
    left_join(gene_type, by='gene_type') %>%
    select(-gene_type)
  
  dbWriteTable(con, "gene_metadata", gene_metadata_final, append = TRUE, overwrite = FALSE)
  
  gene_id_key = gene_metadata_final %>%
    select(gene_id, ensembl_id)
  
  return(gene_id_key)
}

preprocess_tumor_type = function(tumor_type) {
  
  project = paste0('TCGA-', tumor_type)
  
  query = GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    sample.type = c("Primary Tumor", "Solid Tissue Normal"),
    access = "open"
  )
  
  data = GDCprepare(query, directory = GDC_DATA, summarizedExperiment = TRUE)
  
  sample_metadata = as.data.frame(colData(data)) %>% 
    tibble::rownames_to_column("rownames") %>% 
    mutate(overall_survival = ifelse(vital_status == "Alive", days_to_last_follow_up, days_to_death)) %>% 
    mutate(deceased = ifelse(vital_status == "Alive", FALSE, TRUE)) %>% 
    select(tcga_sample_id=sample_id, barcode, sample, patient, sample_type, overall_survival, deceased, preservation_method, age_at_diagnosis, site_of_resection_or_biopsy, race, gender, ethnicity) %>%
    mutate(tumor_type = tumor_type)
  
  saveRDS(sample_metadata, file.path("temp", "sample_metadata", paste0(tumor_type, ".rds")))
  
  expected_count = as.data.frame(assay(data, "unstranded")) %>%
    tibble::rownames_to_column('gene_id') %>%
    pivot_longer(!gene_id, names_to='barcode', values_to='expected_count')
  
  saveRDS(expected_count, file.path("temp", "expected_count", paste0(tumor_type, ".rds")))
  remove(expected_count)
  
  tpm = as.data.frame(assay(data, "tpm_unstrand")) %>%
    tibble::rownames_to_column('gene_id') %>%
    pivot_longer(!gene_id, names_to='barcode', values_to='tpm')
  
  saveRDS(tpm, file.path("temp", "tpm", paste0(tumor_type, ".rds")))
  remove(tpm)
  remove(data)
  gc()
  
  query = GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Gene Level Copy Number",
    workflow.type = "ASCAT3",
    access="open"
  )
  
  data = GDCprepare(query, directory = GDC_DATA, summarizedExperiment = TRUE)
  
  cnv = data.frame(data@assays@data$copy_number)
  colnames(cnv) = colData(data)@listData$sample_id
  rownames(cnv) = data@rowRanges@elementMetadata@listData$gene_id
  
  remove(data)
  gc()
  
  cnv_final = cnv %>%
    tibble::rownames_to_column('gene_id') %>%
    pivot_longer(!gene_id, names_to='tcga_sample_id', values_to='copy_number') %>%
    select(gene_id, tcga_sample_id, copy_number)
  
  saveRDS(cnv_final, file.path("temp", "copy_number", paste0(tumor_type, ".rds")))
  rm(cnv_final)
  gc()
}

build_sample_metadata = function(con) {
  
  first_tumor_type = TUMOR_TYPES[[1]]
  other_tumor_types = TUMOR_TYPES[-1]
  
  df = readRDS(file.path("temp", "sample_metadata", paste0(first_tumor_type, ".rds")))
  
  for (tumor_type in other_tumor_types) {
    df = df %>% add_row(readRDS(file.path("temp", "sample_metadata", paste0(tumor_type, ".rds"))))
  }
  
  sample_metadata = df %>% 
    arrange(tumor_type, barcode) %>%
    mutate(barcode_id = row_number()) %>%
    relocate(barcode_id)
    
  dbWriteTable(con, "sample_metadata", sample_metadata, append = TRUE, overwrite = FALSE)
  
  barcode_key = sample_metadata %>%
    select(barcode_id, barcode)
  
  return(barcode_key)
  
}

build_temp_gene_data = function(con, barcodes, gene_id_key) {
  
  for (tumor_type in TUMOR_TYPES) {
    print(tumor_type)
    expected_count = readRDS(file.path("temp", "expected_count", paste0(tumor_type, ".rds")))
    tpm = readRDS(file.path("temp", "tpm", paste0(tumor_type, ".rds")))
    cn = readRDS(file.path("temp", "copy_number", paste0(tumor_type, ".rds")))
    sample_metadata = readRDS(file.path("temp", "sample_metadata", paste0(tumor_type, ".rds")))
    
    gene_data = expected_count %>% 
      left_join(tpm, by = c('gene_id', 'barcode')) %>%
      left_join(sample_metadata %>% select(barcode, tcga_sample_id), by = 'barcode') %>%
      left_join(cn, by=c('gene_id', 'tcga_sample_id')) %>%
      select(ensembl_id=gene_id, barcode, expected_count, tpm, copy_number) %>%
      left_join(barcodes, by = 'barcode') %>%
      left_join(gene_id_key, by = 'ensembl_id') %>%
      select(gene_id, barcode_id, expected_count, tpm, copy_number) %>%
      arrange(gene_id, barcode_id)

    dbWriteTable(con, "gene_data_temp", gene_data, append = TRUE, overwrite = FALSE)
    
  }
  
}

sort_gene_data = function(con) {
  
  DBI::dbExecute(con, "
    INSERT INTO gene_data
    SELECT gene_id, barcode_id, expected_count, tpm, copy_number
    FROM gene_data_temp
    ORDER BY gene_id, barcode_id
  ")
  DBI::dbExecute(con, "DROP TABLE gene_data_temp")

}

cleanup_database = function(con) {
  DBI::dbExecute(con, "VACUUM;")
  DBI::dbExecute(con, "EXPORT DATABASE 'temp_db' (FORMAT PARQUET)")
}

cleanup_directory = function() {
  
  temp_db = file.path(DB_PATH, DB_NAME_TEMP)
  if (file.exists(temp_db)) {
    file.remove(temp_db)
  }
  
  unlink("temp", recursive = TRUE)
  unlink("temp_db", recursive = TRUE)
  
  if (file.exists("results.rds")) {
    file.remove("results.rds")
  }
  
  if (file.exists("df.rds")) {
    file.remove("df.rds")
  }
  
}

main()