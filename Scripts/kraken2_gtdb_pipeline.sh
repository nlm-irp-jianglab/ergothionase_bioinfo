#!/bin/bash

# ----------- CONFIGURATION ------------
DB_PATH="kraken2"                # Path to GTDB Kraken2 database
READ1="sample_R1.fastq"          # Input FASTQ R1
READ2="sample_R2.fastq"          # Input FASTQ R2
THREADS=16                       # Number of CPU threads
OUTPUT_PREFIX="sample"           # Prefix for all output files
READ_LEN=150                     # Read length for Bracken
CONFIDENCE=0.1                   # Kraken2 confidence threshold
# --------------------------------------

# Step 1: Kraken2 classification
echo "[1/3] Running Kraken2..."
kraken2 \
  --db $DB_PATH \
  --threads $THREADS \
  --paired \
  --use-names \
  --confidence $CONFIDENCE \
  --report ${OUTPUT_PREFIX}.report \
  --output ${OUTPUT_PREFIX}.kraken \
  $READ1 $READ2

# Step 2: Bracken abundance estimation (species level)
echo "[2/3] Running Bracken..."
bracken \
  -d $DB_PATH \
  -i ${OUTPUT_PREFIX}.report \
  -o ${OUTPUT_PREFIX}_bracken_species.tsv \
  -r $READ_LEN \
  -l S

# Step 3: Convert Bracken output to table (cleaned TSV with only useful columns)
echo "[3/3] Extracting tabular data..."
cut -f1,2,5 ${OUTPUT_PREFIX}_bracken_species.tsv > ${OUTPUT_PREFIX}_species_table.tsv

echo "Pipeline complete. Output files:"
echo "- Kraken2 report: ${OUTPUT_PREFIX}.report"
echo "- Bracken TSV: ${OUTPUT_PREFIX}_bracken_species.tsv"
echo "- Table: ${OUTPUT_PREFIX}_species_table.tsv (columns: name, taxonomy_id, fraction_total_reads)"

