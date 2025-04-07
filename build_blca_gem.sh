#!/bin/bash

# Define parameters
MANIFEST_FILE="gdc_manifest.2025-04-02.112803.txt"
OUTPUT_FILE="blca_gene_expression_matrix.tsv"
DOWNLOAD_DIR="gdc_downloads"

# Check if manifest file exists
if [ ! -f "$MANIFEST_FILE" ]; then
    echo "Error: Manifest file $MANIFEST_FILE not found."
    echo "Available manifest files:"
    find . -name "gdc_manifest*"
    exit 1
fi

# Create download directory if it doesn't exist
mkdir -p "$DOWNLOAD_DIR"

# Download files using GDC client (if you have gdc-client installed)
echo "Downloading files from GDC using manifest..."
which gdc-client > /dev/null
if [ $? -eq 0 ]; then
    gdc-client download -m "$MANIFEST_FILE" -d "$DOWNLOAD_DIR"
else
    echo "GDC client not found in PATH. Skipping download."
    echo "Please install GDC client or make sure it's in your PATH"
fi

# Find all downloaded RNA-seq files
echo "Locating RNA-seq files..."
find "$DOWNLOAD_DIR" -name "*.htseq.counts" -o -name "*.FPKM.txt" -o -name "*.rsem.genes.results" > temp_file_list.txt

# Check if any files were found
if [ ! -s temp_file_list.txt ]; then
    echo "No RNA-seq files found. Checking available files:"
    find "$DOWNLOAD_DIR" -type f | head -n 10
    exit 1
fi

echo "Found $(wc -l < temp_file_list.txt) RNA-seq files."

# Create a temporary directory for processing
mkdir -p temp_gem

# Process each RNA-seq file
echo "Processing RNA-seq files..."
while read expr_file; do
    # Extract sample ID from the file path
    # GDC files typically have TCGA barcodes in their directory structure
    
    # Extract TCGA barcode from file path
    tcga_barcode=$(basename "$(dirname "$expr_file")" | grep -o "TCGA-[A-Z0-9]\+-[A-Z0-9]\+-[A-Z0-9]\+-[A-Z0-9]\+-[A-Z0-9]\+-[A-Z0-9]\+")
    
    # If no barcode found, use directory name
    if [ -z "$tcga_barcode" ]; then
        tcga_barcode=$(basename "$(dirname "$expr_file")")
    fi
    
    echo "Processing $tcga_barcode: $expr_file"
    
    # Extract gene IDs and expression values
    # Adjust the awk command based on the file format
    if [[ "$expr_file" == *htseq.counts ]]; then
        # HTSeq format: gene_id count
        awk '{print $1"\t"$2}' "$expr_file" > "temp_gem/${tcga_barcode}.processed.tsv"
    elif [[ "$expr_file" == *FPKM.txt ]]; then
        # FPKM format: gene_id fpkm
        awk '{print $1"\t"$2}' "$expr_file" > "temp_gem/${tcga_barcode}.processed.tsv"
    elif [[ "$expr_file" == *rsem.genes.results ]]; then
        # RSEM format: gene_id expected_count TPM FPKM
        awk '{print $1"\t"$6}' "$expr_file" > "temp_gem/${tcga_barcode}.processed.tsv"
    fi
done < temp_file_list.txt

# Check if any files were processed
if [ ! "$(ls -A temp_gem)" ]; then
    echo "Error: No files were successfully processed. Check file formats and paths."
    rm -rf temp_gem
    exit 1
fi

# Create the gene list (first column of the matrix)
echo "Creating gene list..."
cat temp_gem/*.processed.tsv | cut -f1 | sort | uniq > temp_gem/gene_list.txt

# Create header row with sample IDs
echo -n "Gene_ID" > temp_gem/header.txt
for sample_file in temp_gem/*.processed.tsv; do
    sample_id=$(basename "$sample_file" .processed.tsv)
    echo -ne "\t$sample_id" >> temp_gem/header.txt
done
echo "" >> temp_gem/header.txt

# Build the matrix row by row
echo "Building expression matrix..."
while read gene_id; do
    echo -n "$gene_id" > temp_gem/current_row.txt
    
    for sample_file in temp_gem/*.processed.tsv; do
        # Find expression value for this gene in this sample
        expr_value=$(grep -w "^$gene_id" "$sample_file" | cut -f2)
        if [ -z "$expr_value" ]; then
            expr_value="0"  # Use 0 for missing values
        fi
        echo -ne "\t$expr_value" >> temp_gem/current_row.txt
    done
    
    echo "" >> temp_gem/current_row.txt
    cat temp_gem/current_row.txt >> temp_gem/matrix_body.txt
done < temp_gem/gene_list.txt

# Combine header and data rows
cat temp_gem/header.txt temp_gem/matrix_body.txt > "$OUTPUT_FILE"

# Clean up
rm -rf temp_gem
rm -f temp_file_list.txt

echo "Gene Expression Matrix created: $OUTPUT_FILE"
echo "Matrix dimensions: $(wc -l < "$OUTPUT_FILE") rows × $(head -1 "$OUTPUT_FILE" | tr '\t' '\n' | wc -l) columns"
