# ===============================================
# STEP 1: DATA EXPLORATION 
# ===============================================

rm(list =ls())

# Check the current working directory 
getwd()

# Set the project working directory 
setwd("C:/Users/margo/OneDrive/Escritorio/Master/Análisis_de_datos_y_Estadística/actividadUD2")

# Check if the count matrix file exists 
file.exists("data/GSE294583_RawCount_Matrix_ann.txt")

# Read the first eight lines of the count matrix 
data_preview <- readLines(
  "data/GSE294583_RawCount_Matrix_ann.txt",
  n = 8
)

# Display the file preview 
cat(data_preview, sep = "\n")

# Calculate the file size in megabytes 
file_size <- file.info(
  "data/GSE294583_RawCount_Matrix_ann.txt"
)$size / 1024^2   #$size selecciona el tamaño enn bytes y 1024^2 convierte los bytes en megabytes

# Display the file size
print(file_size)


# =================================================================
# STEP 2: OBTAIN THE COUNT MATRIX (COUNTS.TSV)
# =================================================================

# Read the expression matrix - 
raw_data <- read.delim(        
  "data/GSE294583_RawCount_Matrix_ann.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

#We have verified that the import has worked 

# Check the dimensions of the imported matrix 
dim(raw_data)

# Display the column names - 
colnames(raw_data)

# Check the data types of the expression columns 
count_columns <- c("T1", "T2", "T3", "Con1", "Con2", "Con3")

sapply(raw_data[count_columns], class)

# Display only the six expression columns 
raw_data[, count_columns]

# Count missing values in each expression column 
colSums(is.na(raw_data[count_columns]))

# We keep only a matrix containing the data we need plus the Gene ID
final_matrix<-raw_data[,c("GeneID",count_columns)]
dim(final_matrix)
colnames(final_matrix)

#Save the matrix
write.table(final_matrix, file = "data/counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#Check if the matrix exists
file.exists("data/counts.tsv")

# =================================================================
# STEP 3: OBTAIN THE METADATA FILE
# =================================================================

# We keep only the sample names (excluding GeneID)
samples <- colnames(final_matrix)[-1]

# We create the corresponding groups
groups<- c("Treated", "Treated", "Treated", "Control","Control","Control")

# We combine both into a table
metadata <- data.frame(sample=samples,group=groups)

# We check that everything worked correctly
print(metadata)

#Save it
write.table(metadata, file = "data/metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# We check that the file has been saved correctly
file.exists("data/metadata.tsv")
