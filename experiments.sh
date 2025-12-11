python scripts/cojo_region_dependency.py \
  --gwas-file results/04_loci/hisp_euro_analysis_eur/LOC_100/LOC_100_gwas_ref_match.tsv \
  --bfile-prefix results/04_loci/hisp_euro_analysis_eur/LOC_100/reference_EUR \
  --output-dir results/05_cojo/hisp_euro_analysis_eur/LOC_100/experiments/region_dependency \
  --sample-size 25217 \
  --use-cojo-slct \
  --region "BLK_REGION:8:11000000:12000000" \
  --region "PPP1R3B_REGION:8:8500000:9000000"


Rscript scripts/plot_region_dependency.R   --summary results/05_cojo/hisp_euro_analysis_eur/LOC_100/experiments/region_dependency/region_dependency_summary.tsv   --output-prefix results/05_cojo/hisp_euro_analysis_eur/LOC_100/experiments/region_dependency/region_dependency_scree