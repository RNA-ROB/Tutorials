#!/bin/bash

# Parse command line arguments
INFILE="$1"
VEP="$2"
THREADS="$3"
OUTDIR="$4"

# Extract file prefix from INFILE
PREFIX=$(basename "$INFILE" | rev | cut -d '.' -f2- | rev)

# Prepare full paths to output files
OUTFILE="$OUTDIR/$PREFIX.vep.vcf"
STATS="$OUTDIR/$PREFIX.stats.txt"
ERRLOG="$OUTDIR/$PREFIX.error.log"

# Annotate variants with VEP 
vep --everything \
    --species "homo_sapiens" \
    --assembly "GRCh37" \
    --input_file "$INFILE" \
    --format "vcf" \
    --output_file "$OUTFILE" \
    --stats_text "$STATS" \
    --warning_file "$ERRLOG" \
    --cache --refseq \
    --dir "$VEP" \
    --dir_cache "$VEP" \
    --vcf --exclude_predicted --offline --fork "$THREADS"
