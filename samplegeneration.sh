cd /home/bio244/bio244miguel/bio244data

# 1. Extract all sample IDs directly from the VCF
bcftools query -l /home/bio244/244VCF/phIndoPim_filter5_pruned_forR.vcf.gz > samples.list

# 2. Automatically generate metadata mapping based on sample names (Blue, Green, Violet)
awk '{
  if ($1 ~ /blue/i) pop="Blue";
  else if ($1 ~ /green/i) pop="Green";
  else if ($1 ~ /violet/i) pop="Violet";
  else pop="Other";
  if (pop != "Other") print $1, pop;
}' samples.list > metadata_samples-population.tsv

cat metadata_samples-population.tsv
