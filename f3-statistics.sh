#!/usr/bin/env bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate admixtools
set -e

# ==============================================================================
# 0. Setup & Directory Configuration
# ==============================================================================
WORKDIR="/home/bio244/bio244miguel/bio244data"
F3DIR="${WORKDIR}/bio244f3stats"
VCF_SRC="/home/bio244/244VCF/phIndoPim_filter5_pruned_forR.vcf.gz"

cd "$WORKDIR"
mkdir -p "$F3DIR"

echo "=== Starting f3 Statistic Workflow for Eastern Samar-Southern Leyte ==="

# ==============================================================================
# 1. Sample List & Population Metadata Verification / Generation
# ==============================================================================
if [ ! -f samples.list ]; then
  echo "samples.list not found. Extracting sample IDs directly from VCF..."
  gzip -dc "$VCF_SRC" | grep -m1 "^#CHROM" | cut -f10- | tr '\t' '\n' > samples.list
fi

if [ ! -f metadata_samples-population.tsv ]; then
  echo "metadata_samples-population.tsv not found. Generating mapping..."
  awk '{
    if ($1 ~ /blue/i) pop="Blue";
    else if ($1 ~ /green/i) pop="Green";
    else if ($1 ~ /violet/i) pop="Violet";
    else pop="Other";
    if (pop != "Other") print $1, pop;
  }' samples.list > metadata_samples-population.tsv
fi

# ==============================================================================
# 2. Filter Samples (Exclude 'Red' Lineage) & Prepare Input Lists
# ==============================================================================
grep -f samples.list metadata_samples-population.tsv | grep -v "Red" | awk '{print $1, $1}' > "${F3DIR}/samples.keep"
grep -f samples.list metadata_samples-population.tsv | grep -v "Red" > "${F3DIR}/samples_metadata.list"

# ==============================================================================
# 3. PLINK Binary File Generation
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
# 4. Recode Scaffold Names to Integer IDs & Run convertf (EIGENSTRAT)
# ==============================================================================
awk 'BEGIN{chr_count=0} {
  if (!($1 in chr_map)) {
    chr_count++;
    chr_map[$1] = chr_count;
  }
  $1 = chr_map[$1];
  print $0;
}' OFS="\t" plink_binary.bim > phIndoPim_convertf.bim

cat > convertf.params << EOF
genotypename:    plink_binary.bed
snpname:         phIndoPim_convertf.bim
indivname:       plink_binary.fam
outputformat:    EIGENSTRAT
genotypeoutname: phIndoPim.geno
snpoutname:      phIndoPim.snp
indivoutname:    phIndoPim.ind
familynames:     NO
EOF

convertf -p convertf.params

# ==============================================================================
# 5. Label EIGENSTRAT Individuals with Population Colors & Run qp3Pop
# ==============================================================================
awk 'NR==FNR {pop[$1]=$2; next} {print $1, $2, pop[$1]}' samples_metadata.list phIndoPim.ind > phIndoPim_labeled.ind

cat > list.pop << EOF
Blue Green Violet
Green Violet Blue
Violet Blue Green
EOF

cat > qp3pop.params << EOF
genotypename: phIndoPim.geno
snpname:      phIndoPim.snp
indivname:    phIndoPim_labeled.ind
popfilename:  list.pop
inbreed:      NO
outgroupmode: NO
EOF

qp3Pop -p qp3pop.params > qp3pop_EasternSamarSouthernLeyte.results

echo "=== f3 Statistics Calculation Complete ==="
echo "Results saved to: ${F3DIR}/qp3pop_EasternSamarSouthernLeyte.results"
cat qp3pop_EasternSamarSouthernLeyte.results
