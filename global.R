GDC_DATA = file.path("data", "raw", "GDC")
DB_PATH = file.path("data", "db")
DB_NAME_TEMP = 'tcga_temp.duckdb'
DB_NAME = 'tcga.duckdb'

# No LAML
TUMOR_TYPES = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")