# TCGA Database

Download publicly available data from TCGA and build a DuckDB database that stores select sample metadata as well as gene expression and copy number information.

## Quick Start

1. Configure paths and tumor types of interest in `global.R`

2. Download raw TCGA data (~82 Gb, may take hours)

    ```
    Rscript bin/download_gdc_data.R
    ```
    
3. Build DuckDB database (~2.6 Gb final)

    ```
    Rscript bin/build_database.R
    ```

## Requirements

- R (>= 4.x)
- Packages:
  - TCGAbiolinks
  - SummarizedExperiment
  - dplyr
  - tidyr
  - duckdb

## Configuration

All configuration is defined in `global.R`

| Variable     | Description                                                 | Default            |
|--------------|-------------------------------------------------------------|--------------------|
| GDC_DATA     | Directory for downloaded raw TCGA data                      | `data/raw/GDC`     |
| DB_PATH      | Final database directory                                    | `data/db/`         |
| DB_NAME      | Final database name                                         | `tcga.duckdb`      |
| DB_NAME_TEMP | Intermediate database name (automatically removed)          | `tcga_temp.duckdb` |
| TUMOR_TYPES  | Vector of TCGA  tumor types to include (e.g. "ACC", "BRCA") | All                |

## Data Download

The following script will download raw TCGA data for tumor types listed in global.R config (default all).

```
Rscript bin/download_gdc_data.R
```

Due to TCGA upload capacity, this script will take several hours to complete. Occasionally TCGA will drop the connection. Restarting the script will not attempt to re-download files that were successfully downloaded prior to the error. At the time of writing the total size of the raw TCGA data is about 82 Gb. Once the database has been built, the raw data can be removed. However, if the raw data is kept, updating the database will not require re-downloading.

## Build Database

The following script will build a database (duckdb) that can be stored as a single file:

```
Rscript bin/build_database.R
```

The final database will be about 2.6 Gb. During the database building process, there will be an additional 10 Gb of disk space that is used temporarily - it will be cleaned-up as part of the script.

## Database Schema

The final database will contain the following tables and corresponding columns:

### gene_metadata

| **Column**   | **Type**    | **Example**        | **Notes**   |
|:-------------|:------------|:-------------------|:------------| 
| gene_id      | *USMALLINT* | 1                  | Internal ID |
| ensembl_id   | *VARCHAR*   | ENSG00000000003.15 |             |
| gene_name    | *VARCHAR*   | TSPAN6             |             |
| gene_type_id | *USMALLINT* | 1                  | Internal ID |

### gene_type_key

| **Column**   | **Type**    | **Example**    | **Notes**   |
|:-------------|:------------|:---------------|:------------| 
| gene_type_id | *USMALLINT* | 1              | Internal ID |
| gene_type    | *VARCHAR*   | protein_coding |             |

### sample_metadata
 
| **Column**                  | **Type**    | **Example**                          | **Notes**   |
|:----------------------------|:------------|:-------------------------------------|:------------| 
| tumor_type                  | *VARCHAR*   | ACC                                  |             |
| barcode_id                  | *USMALLINT* | 1                                    | Internal ID |
| barcode                     | *VARCHAR*   | TCGA-OR-A5J1-01A-11R-A29S-07         |             |
| sample                      | *VARCHAR*   | TCGA-OR-A5J1-01A                     |             |
| patient                     | *VARCHAR*   | TCGA-OR-A5J1                         |             |
| tcga_sample_id              | *VARCHAR*   | e4038ebb-6e6d-44b1-84ad-e35aafca7b70 |             |
| sample_type                 | *VARCHAR*   | Primary Tumor                        |             |
| overall_survival            | *INTEGER*   | 1355                                 |             |
| deceased                    | *BOOLEAN*   | TRUE                                 |             |
| preservation_method         | *VARCHAR*   | OCT                                  |             |
| age_at_diagnosis            | *INTEGER*   | 22250                                |             |
| site_of_resection_or_biopsy | *VARCHAR*   | Adrenal gland, NOS                   |             |
| race                        | *VARCHAR*   | white                                |             |
| gender                      | *VARCHAR*   | female                               |             |
| ethnicity                   | *VARCHAR*   | hispanic or latino                   |             |

### gene_data

| **Column**     | **Type**    | **Example** | **Notes**   |
|:---------------|:------------|:------------|:------------| 
| gene_id        | *USMALLINT* | 1           | Internal ID |
| barcode_id     | *USMALLINT* | 4           | Internal ID |
| expected_count | *INTEGER*   | 1396        |             |
| tpm            | *FLOAT*     | 17.1286     |             |
| copy_number    | *UTINYINT*  | 4           |             |

## Example Usage

```
library(duckdb)
con = dbConnect(duckdb::duckdb(), dbdir = file.path(DB_PATH, DB_NAME))
```

Get the expected counts for TP53 in all ACC samples:

```
df = 
  dbGetQuery(con, "
    SELECT 
      sm.barcode,
      sm.sample_type,
      gd.expected_count
    FROM gene_data gd
    JOIN gene_metadata gm 
      ON gd.gene_id = gm.gene_id
    JOIN sample_metadata sm 
      ON gd.barcode_id = sm.barcode_id
    WHERE gm.gene_name = 'TP53'
      AND sm.tumor_type = 'ACC'
  ")
```

Get the TPM  of all protein coding genes from all BRCA Primary Tumors:

```
df = 
  dbGetQuery(con, "
    SELECT 
      gm.gene_name,
      gm.ensembl_id,
      sm.barcode,
      gd.tpm
    FROM gene_data gd
    JOIN gene_metadata gm 
      ON gd.gene_id = gm.gene_id
    JOIN gene_type_key gtk 
      ON gm.gene_type_id = gtk.gene_type_id
    JOIN sample_metadata sm 
      ON gd.barcode_id = sm.barcode_id
    WHERE sm.tumor_type = 'BRCA'
      AND sm.sample_type = 'Primary Tumor'
      AND gtk.gene_type = 'protein_coding'
    ORDER BY gm.gene_name, sm.barcode;
  ")
```


## Notes / Limitations

* The database currently limits selected sample metadata. The columns can be modified in the the `preprocess_tumor_type` function in `build_database.R`
* Currently, copy number is limited to samples that have expression data.

## License

This project is licensed under the MIT License.

Note: This repository does not distribute TCGA data. Users are responsible for complying with GDC data usage policies.