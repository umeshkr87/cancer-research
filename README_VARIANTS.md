# 🧬 Variant Processing Pipeline

## Quick Setup
```bash
# 1. Upload COSMIC file
scp Cosmic_MutantCensus_v102_GRCh38.tsv aa107@cli-dtn.researchdata.illinois.edu:/u/aa107/uiuc-cancer-research/data/raw/variants/

# 2. Run automated processing
cd /u/aa107/uiuc-cancer-research/scripts/variants/
chmod +x run_cosmic_filter.sh
./run_cosmic_filter.sh