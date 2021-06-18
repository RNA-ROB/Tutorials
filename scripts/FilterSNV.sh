#!/bin/bash

# This is a Bash script that will extract SNVs from an input ClinVar VCF file and segregates the SNVs
# into three user-defined categories of "variant confidence" based on ClinVar star ratings

# Bash function that filters for SNVs belonging to a particular
# clinicial significance category given a VCF
FilterSNV() {
    INFILE="$1"
    CLASS="$2"
    awk -v "CLASS=$CLASS" 'BEGIN{FS="\t";}{
        # Iterate over INFO field
        n=split($8,INFO,";");
        for(i=1;i<=n;i++){
			if(INFO[i]~"CLNSIG="){
				split(INFO[i],CLNSIG,"=");
				if(CLNSIG[2]==CLASS){
					# Check if (i) variant is SNV and (ii) alternate allele indeed exists
					if($5!="." && length($4)==1 && length($5)==1){
						printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,".");
					}
				}
			}
		}
    }' "$INFILE"
}
export -f FilterSNV

# Parse command-line arguments
INPUT="$1"
RATINGS="$2"
LOWER="$3"
UPPER="$4"
OUTDIR="$5"

# Report user-defined parameters used in running FilterSNV.sh
echo -e "Input VCF file: $INPUT" > "$OUTDIR/parameters.txt"
echo -e "Output directory: $OUTDIR" >> "$OUTDIR/parameters.txt"
echo -e "ClinVar Star Rating Threshold (Low-Confidence Variants): [0,$LOWER)" >> "$OUTDIR/parameters.txt"
echo -e "ClinVar Star Rating Threshold (Medium-Confidence Variants): [$LOWER,$UPPER]" >> "$OUTDIR/parameters.txt"
echo -e "ClinVar Star Rating Threshold (High-Confidence Variants): ($UPPER,4]" >> "$OUTDIR/parameters.txt"

# Create appropriate directories
mkdir -p "$OUTDIR/high_confidence" 
mkdir -p "$OUTDIR/medium_confidence" 
mkdir -p "$OUTDIR/low_confidence"

# Retrieve base name prefix of VCF file
PREFIX=$(basename "$INPUT" | rev | cut -d '.' -f2- | rev)

# Filter for variants with both CLNSIG and CLNREVSTAT annotations
grep -Pv "^#" "$INPUT" | awk 'BEGIN{FS="\t";}{
	# Define indicator variables for whether variant annotation comes with both CLNSIG and CLNREVSTAT
	CTR1=0; CTR2=0;
	
	# Process INFO field
	n=split($8,INFO,";");
	for(i=1;i<=n;i++){
		# Check if variant contains CLNREVSTAT field
		if(INFO[i]~"CLNREVSTAT="){
			CTR1=1;
		}
		# Check if variant contains CLNSIG field
		if(INFO[i]~"CLNSIG="){
			CTR2=1;
		}
	}

    if(CTR1==1 && CTR2==1){
        print;
    }
}' > "$OUTDIR/$PREFIX.filtered.vcf"

# Append a ClinVar star rating for each variant based on CLNREVSTAT field string
awk '{
	if(FNR==NR){
		# Build dictionary that maps CLNREVSTAT status to star rating based on ClinVar guidelines
		STAR[$1]=$2;
	}
	else{
		# Record number of ClinVar stars for each variant
		NSTAR="NA";
		
		# Process INFO field
		n=split($8,INFO,";");
		for(i=1;i<=n;i++){
			if(INFO[i]~"CLNREVSTAT="){
				split(INFO[i],CLNREVSTAT,"=");
				NSTAR=STAR[CLNREVSTAT[2]];
			}
		}
		
		# Print number of ClinVar stars for each variant
		printf("%s\t%s\n",NSTAR,$0);
	}
}' "$RATINGS" "$OUTDIR/$PREFIX.filtered.vcf" > "$OUTDIR/$PREFIX.filtered.star.vcf"

# Sort variants into each group based on defined thresholds
awk -v "LOWER=$LOWER" 'BEGIN{FS="\t";}{if($1<LOWER){print;}}' "$OUTDIR/$PREFIX.filtered.star.vcf" | cut -f2- > "$OUTDIR/low_confidence/$PREFIX.low_confidence.vcf"
awk -v "LOWER=$LOWER" -v "UPPER=$UPPER" 'BEGIN{FS="\t";}{if($1>=LOWER && $1<=UPPER){print;}}' "$OUTDIR/$PREFIX.filtered.star.vcf" | cut -f2- > "$OUTDIR/medium_confidence/$PREFIX.medium_confidence.vcf"
awk -v "UPPER=$UPPER" 'BEGIN{FS="\t";}{if($1>UPPER){print;}}' "$OUTDIR/$PREFIX.filtered.star.vcf" | cut -f2- > "$OUTDIR/high_confidence/$PREFIX.high_confidence.vcf"

## For each confidence level, isolate variants in each of the clinical significance categories:
# 	Benign	B
# 	Benign/Likely_benign	BLB
# 	Conflicting_interpretations_of_pathogenicity	CIP
# 	Likely_benign	LB
# 	Likely_pathogenic	LP
# 	Pathogenic	P
# 	Pathogenic/Likely_pathogenic	PLP
# 	Uncertain_significance	VUS

for GROUP in "low_confidence" "medium_confidence" "high_confidence"; do
	# Isolate variants from each class
	FilterSNV "$OUTDIR/$GROUP/$PREFIX.$GROUP.vcf" "Benign" > "$OUTDIR/$GROUP/$PREFIX.$GROUP.B.vcf"
	FilterSNV "$OUTDIR/$GROUP/$PREFIX.$GROUP.vcf" "Benign/Likely_benign" > "$OUTDIR/$GROUP/$PREFIX.$GROUP.BLB.vcf"
	FilterSNV "$OUTDIR/$GROUP/$PREFIX.$GROUP.vcf" "Conflicting_interpretations_of_pathogenicity" > "$OUTDIR/$GROUP/$PREFIX.$GROUP.CIP.vcf"
	FilterSNV "$OUTDIR/$GROUP/$PREFIX.$GROUP.vcf" "Likely_benign" > "$OUTDIR/$GROUP/$PREFIX.$GROUP.LB.vcf"
	FilterSNV "$OUTDIR/$GROUP/$PREFIX.$GROUP.vcf" "Likely_pathogenic" > "$OUTDIR/$GROUP/$PREFIX.$GROUP.LP.vcf"
	FilterSNV "$OUTDIR/$GROUP/$PREFIX.$GROUP.vcf" "Pathogenic" > "$OUTDIR/$GROUP/$PREFIX.$GROUP.P.vcf"
	FilterSNV "$OUTDIR/$GROUP/$PREFIX.$GROUP.vcf" "Pathogenic/Likely_pathogenic" > "$OUTDIR/$GROUP/$PREFIX.$GROUP.PLP.vcf"
	FilterSNV "$OUTDIR/$GROUP/$PREFIX.$GROUP.vcf" "Uncertain_significance" > "$OUTDIR/$GROUP/$PREFIX.$GROUP.VUS.vcf"
done

# Remove intermediate files
rm "$OUTDIR/$PREFIX.filtered.vcf"
rm "$OUTDIR/$PREFIX.filtered.star.vcf"
rm "$OUTDIR/low_confidence/$PREFIX.low_confidence.vcf"
rm "$OUTDIR/medium_confidence/$PREFIX.medium_confidence.vcf"
rm "$OUTDIR/high_confidence/$PREFIX.high_confidence.vcf"
