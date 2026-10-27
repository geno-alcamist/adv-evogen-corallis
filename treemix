#!/usr/bin/env bash
set -e

# ==============================================================================
# 0. Setup & Directory Configuration
# ==============================================================================
WORKDIR="/home/bio244/bio244miguel/bio244data"
TREEDIR="${WORKDIR}/bio244treemix"
F3DIR="${WORKDIR}/bio244f3stats"
VCF_SRC="/home/bio244/244VCF/phIndoPim_filter5_pruned_forR.vcf.gz"

cd "$WORKDIR"
mkdir -p "${TREEDIR}/input" "${TREEDIR}/sample_lists" "${TREEDIR}/results" "${TREEDIR}/scripts"

echo "=== Starting TreeMix Workflow for Eastern Samar-Southern Leyte ==="

# ==============================================================================
# 1. Prepare Population Lists & PLINK Clusters
# ==============================================================================
# Copy or build sample metadata list from f3 directory if available
if [ -f "${F3DIR}/samples_metadata.list" ]; then
  cp "${F3DIR}/samples_metadata.list" "${TREEDIR}/sample_lists/easternsamar_southernleyte_popmap.tsv"
else
  grep -f samples.list metadata_samples-population.tsv | grep -v "Red" > "${TREEDIR}/sample_lists/easternsamar_southernleyte_popmap.tsv"
fi

# Build PLINK cluster file: FID (0), IID ($1), CLUSTER ($2)
awk '{print 0, $1, $2}' "${TREEDIR}/sample_lists/easternsamar_southernleyte_popmap.tsv" > "${TREEDIR}/sample_lists/easternsamar_southernleyte_clusters.txt"

# Build population order list
awk '{print $2}' "${TREEDIR}/sample_lists/easternsamar_southernleyte_popmap.tsv" | sort -u > "${TREEDIR}/sample_lists/pop_order.txt"

# Build PLINK keep list
awk '{print 0, $1}' "${TREEDIR}/sample_lists/easternsamar_southernleyte_popmap.tsv" > "${TREEDIR}/sample_lists/samples.keep"

# ==============================================================================
# 2. VCF to PLINK Conversion (plink2)
# ==============================================================================
plink2 \
  --vcf "$VCF_SRC" \
  --keep "${TREEDIR}/sample_lists/samples.keep" \
  --make-bed \
  --allow-extra-chr \
  --out "${TREEDIR}/input/easternsamar_southernleyte"

# ==============================================================================
# 3. Stratified Allele-Frequency Calculation (PLINK 1.9)
# ==============================================================================
plink \
  --bfile "${TREEDIR}/input/easternsamar_southernleyte" \
  --within "${TREEDIR}/sample_lists/easternsamar_southernleyte_clusters.txt" \
  --freq gz \
  --allow-extra-chr \
  --out "${TREEDIR}/input/easternsamar_southernleyte_freq"

# ==============================================================================
# 4. Conversion to TreeMix Input Format
# ==============================================================================
# Create python conversion script in scripts/ convert_treemix_validated.py
cat > "${TREEDIR}/scripts/convert_treemix.py" << 'EOF'
import gzip, sys

freq_file = sys.argv[1]
pop_order_file = sys.argv[2]
out_file = sys.argv[3]

with open(pop_order_file, 'r') as f:
    pops = [line.strip() for line in f if line.strip()]

counts = {}
snps = []

with gzip.open(freq_file, 'rt') as infile:
    header = infile.readline()
    for line in infile:
        parts = line.strip().split()
        if len(parts) < 7: continue
        chr_id, snp, clst, a1, a2, mac, nchrobs = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
        
        key = f"{chr_id}:{snp}"
        if key not in counts:
            counts[key] = {}
            snps.append(key)
        
        ref_count = int(mac)
        alt_count = int(nchrobs) - int(mac)
        counts[key][clst] = f"{ref_count},{alt_count}"

with gzip.open(out_file, 'wt') as outfile:
    outfile.write(" ".join(pops) + "\n")
    for key in snps:
        row = [counts[key].get(p, "0,0") for p in pops]
        outfile.write(" ".join(row) + "\n")
EOF

python3 "${TREEDIR}/scripts/convert_treemix.py" \
  "${TREEDIR}/input/easternsamar_southernleyte_freq.frq.strat.gz" \
  "${TREEDIR}/sample_lists/pop_order.txt" \
  "${TREEDIR}/input/easternsamar_southernleyte.treemix.frq.gz"

# ==============================================================================
# 5. TreeMix Execution (m = 0, 1, 2 migration events)
# ==============================================================================
cd "$TREEDIR"

for m in 0 1 2; do
  echo "Running TreeMix for Eastern Samar-Southern Leyte with m=${m}..."
  treemix \
    -i input/easternsamar_southernleyte.treemix.frq.gz \
    -o "results/easternsamar_southernleyte_m${m}" \
    -m "$m" \
    -k 500 \
    -seed 12345 || echo "Warning: TreeMix run m=${m} returned non-zero exit code."
done

echo "=== TreeMix Analysis Complete ==="
echo "Results saved to ${TREEDIR}/results/"
