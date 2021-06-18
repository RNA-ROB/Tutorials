#!/bin/bash

# This is a general script for calling variants from a WES/WGS read alignment file (BAM)
# using GATK4 Best Practices

# Parse command line arguments
INPUT="$1" # This is your input BAM file of WES/WGS reads
REFERENCE="$2" # This is your reference genome (make sure you have indexed this FASTA file and generated a sequence dictionary)
SAMPID="$3" # This is your sample ID
OUTDIR="$4" # This is your output directory

# Create a temporary directory in OUTDIR to store intermediate files
TMPDIR="$OUTDIR/tmp"
mkdir -p "$TMPDIR"

## GATK4 Best Practices Pipeline
# Variant calling (Round 1)
gatk HaplotypeCaller -R "$REFERENCE" -I "$INPUT" -O "$TMPDIR/$SAMPID.Round_1.vcf.gz"

for TYPE in "SNP" "INDEL"; do
    # Separate SNPs and INDELs
    gatk SelectVariants -R "$REFERENCE" -V "$TMPDIR/$SAMPID.Round_1.vcf.gz" -select-type "$TYPE" -O "$TMPDIR/$SAMPID.Round_1.$TYPE.vcf.gz"

    # Filter SNPs and INDELs based on thresholds recommended by the Broad
    if [ "$TYPE" == "SNP" ]; then
        # Recommended filters for SNPs
        gatk VariantFiltration -R "$REFERENCE" -V "$TMPDIR/$SAMPID.Round_1.$TYPE.vcf.gz" \
            -O "$TMPDIR/$SAMPID.Round_1.$TYPE.filtered.vcf.gz" \
            --filter-name "QD_filter" -filter "QD < 2.0" \
            --filter-name "FS_filter" -filter "FS > 60.0" \
            --filter-name "MQ_filter" -filter "MQ < 40.0" \
            --filter-name "SOR_filter" -filter "SOR > 4.0" \
            --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
            --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    else
        # Recommended filters for INDELs
        gatk VariantFiltration -R "$REFERENCE" -V "$TMPDIR/$SAMPID.Round_1.$TYPE.vcf.gz" \
            -O "$TMPDIR/$SAMPID.Round_1.$TYPE.filtered.vcf.gz" \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 200.0" \
            -filter-name "SOR_filter" -filter "SOR > 10.0"
    fi

    # Extract variants passing filter thresholds (needed for BQSR)
    gatk SelectVariants --exclude-filtered -V "$TMPDIR/$SAMPID.Round_1.$TYPE.filtered.vcf.gz" -O "$TMPDIR/$SAMPID.Round_1.$TYPE.BQSR.vcf.gz"

done

# Perform BQSR to generate a recalibration data table
gatk BaseRecalibrator -R "$REFERENCE" -I "$INPUT" --known-sites "$TMPDIR/$SAMPID.Round_1.SNP.BQSR.vcf.gz" \
    --known-sites "$TMPDIR/$SAMPID.Round_1.INDEL.BQSR.vcf.gz" -O "$TMPDIR/$SAMPID.BQSR.table"

# Generate a BAM file of recalibrated reads
gatk ApplyBQSR -R "$REFERENCE" -I "$INPUT" -bqsr "$TMPDIR/$SAMPID.BQSR.table" -O "$TMPDIR/$SAMPID.BQSR.bam"

# Variant calling (Round 2)
gatk HaplotypeCaller -R "$REFERENCE" -I "$TMPDIR/$SAMPID.BQSR.bam" -O "$TMPDIR/$SAMPID.Round_2.vcf.gz"

for TYPE in "SNP" "INDEL"; do
    # Separate SNPs and INDELs
    gatk SelectVariants -R "$REFERENCE" -V "$TMPDIR/$SAMPID.Round_2.vcf.gz" -select-type "$TYPE" -O "$TMPDIR/$SAMPID.Round_2.$TYPE.vcf.gz"

    # Filter SNPs and INDELs based on thresholds recommended by the Broad
    if [ "$TYPE" == "SNP" ]; then
        # Recommended filters for SNPs
        gatk VariantFiltration -R "$REFERENCE" -V "$TMPDIR/$SAMPID.Round_2.$TYPE.vcf.gz" \
            -O "$TMPDIR/$SAMPID.Round_2.$TYPE.filtered.vcf.gz" \
            --filter-name "QD_filter" -filter "QD < 2.0" \
            --filter-name "FS_filter" -filter "FS > 60.0" \
            --filter-name "MQ_filter" -filter "MQ < 40.0" \
            --filter-name "SOR_filter" -filter "SOR > 4.0" \
            --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
            --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    else
        # Recommended filters for INDELs
        gatk VariantFiltration -R "$REFERENCE" -V "$TMPDIR/$SAMPID.Round_2.$TYPE.vcf.gz" \
            -O "$TMPDIR/$SAMPID.Round_2.$TYPE.filtered.vcf.gz" \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 200.0" \
            -filter-name "SOR_filter" -filter "SOR > 10.0"
    fi

    # Extract variants passing filter thresholds
    gatk SelectVariants --exclude-filtered -V "$TMPDIR/$SAMPID.Round_2.$TYPE.filtered.vcf.gz" -O "$OUTDIR/$SAMPID.$TYPE.vcf.gz"

done

# Clean up TMPDIR
rm $TMPDIR/*
rmdir $TMPDIR
