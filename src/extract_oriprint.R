library(Seurat)
library(dplyr)

# 1. Load the data you downloaded
message("Reading Lo Riso & Villa RDS files...")
seu <- readRDS("CellOfOriginDataset.rds")
anno <- readRDS("AnnotationCellOfOrigin.rds")

# 2. Extract the 'OriPrint' Signatures
# These are the genes that define Ovarian vs Fallopian origin
message("Extracting Gene Signatures...")
signatures <- as.data.frame(as.matrix(GetAssayData(seu, slot = "data")))
# Keep only the top variable genes used for the OriPrint classifier
top_genes <- VariableFeatures(seu)
signature_subset <- signatures[top_genes, ]

# 3. Export to OvaLens data folder
dir.create("data/external", showWarnings = FALSE)
write.csv(signature_subset, "data/external/testa_oriprint_signatures.csv")
write.csv(seu@meta.data, "data/external/testa_origin_metadata.csv")

message("✅ Success! Files moved to data/external/")
