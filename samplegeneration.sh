cd /home/bio244/bio244miguel/bio244data

# 1. Extract all sample IDs directly from the VCF header (no bcftools required)
gzip -dc /home/bio244/244VCF/phIndoPim_filter5_pruned_forR.vcf.gz | grep -m1 "^#CHROM" | cut -f10- | tr '\t' '\n' > samples.list

# 2. Generate population metadata mapping (Blue, Green, Violet)
awk '{
  if ($1 ~ /blue/i) pop="Blue";
  else if ($1 ~ /green/i) pop="Green";
  else if ($1 ~ /violet/i) pop="Violet";
  else pop="Other";
  if (pop != "Other") print $1, pop;
}' samples.list > metadata_samples-population.tsv

# 3. Check the created sample list and population counts
cat metadata_samples-population.tsv
