#!/usr/bin/env bash
set -e

# ==============================================================================
# 0. Setup & Directory Configuration
# ==============================================================================
WORKDIR="/home/bio244/bio244miguel/bio244data"
F3DIR="${WORKDIR}/bio244f3stats"
VCF_SRC="/home/bio244/244VCF/phIndoPim_filter5_pruned_forR.vcf.gz"

cd "$WORKDIR"
mkdir -p "$F3DIR"

echo "=== Starting f3 Statistic Workflow ==="

# ==============================================================================
# 1. Sample Selection & Metadata Subsetting (Eastern Samar-Southern Leyte)
# ==============================================================================
# Ensure samples.list and metadata_samples-population.tsv exist in $WORKDIR
if [ ! -f samples.list ] || [ ! -f metadata_samples-population.tsv ]; then
  echo "Error: Required 'samples.list' or 'metadata_samples-population.tsv' not found in $WORKDIR"
  exit 1
fi

# Filter out 'Red' individuals and generate PLINK keep file and metadata list
grep -f samples.list metadata_samples-population.tsv | grep -v "Red" | awk '{print $1, $1}' > "${F3DIR}/samples.keep"
grep -f samples.list metadata_samples-population.tsv | grep -v "Red" > "${F3DIR}/samples_metadata.list"

echo "Samples selected for f3 analysis:"
cat "${F3DIR}/samples_metadata.list"

# ==============================================================================
# 2. PLINK Binary File Generation
# ==============================================================================
cd "$F3DIR"

plink \
  --vcf "$VCF_SRC" \
  --allow-extra-chr \
  --chr host_JAEMWC010000002.1 host_JAEMWC010000009.1 host_JAEMWC010000016.1 host_JAEMWC010000023.1 \
  --keep samples.keep \
  --make-bed \
  --out plink_binary

# ==============================================================================
# 3. Convert PLINK Binary to EIGENSTRAT Format (using convertf)
# ==============================================================================
cat > convertf.params << EOF
genotypename:    plink_binary.bed
snpname:         plink_binary.bim
indivname:       plink_binary.fam
outputformat:    EIGENSTRAT
genotypeoutname: phIndoPim.geno
snpoutname:      phIndoPim.snp
indivoutname:    phIndoPim.ind
familynames:     NO
EOF

convertf -p convertf.params

# ==============================================================================
# 4. Label EIGENSTRAT Individuals with Population Colors (Blue, Green, Violet)
# ==============================================================================
awk 'NR==FNR {pop[$1]=$2; next} {print $1, $2, pop[$1]}' samples_metadata.list phIndoPim.ind > phIndoPim_labeled.ind

echo "Population breakdown in labeled .ind file:"
cut -d " " -f 3 phIndoPim_labeled.ind | sort | uniq -c

# ==============================================================================
# 5. Configure & Run qp3Pop
# ==============================================================================
# Define target/source triplets
cat > list.pop << EOF
Blue Green Violet
Green Violet Blue
Violet Blue Green
EOF

# Build parameter file for qp3Pop
cat > qp3pop.params << EOF
genotypename: phIndoPim.geno
snpname:      phIndoPim.snp
indivname:    phIndoPim_labeled.ind
popfilename:  list.pop
inbreed:      NO
outgroupmode: NO
EOF

# Execute qp3Pop
qp3Pop -p qp3pop.params > qp3pop_EasternSamarSouthernLeyte.results

echo "=== f3 Statistics Analysis Complete ==="
echo "Results written to: ${F3DIR}/qp3pop_EasternSamarSouthernLeyte.results"
cat qp3pop_EasternSamarSouthernLeyte.results
