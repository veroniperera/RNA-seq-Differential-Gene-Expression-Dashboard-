# RNA-seq Differential Gene Expression Dashboard

An interactive dashboard for exploring differential gene expression, built on a full RNA-seq analysis pipeline: DESeq2 → GO/KEGG enrichment → Streamlit.

**Dataset:** [GSE110114](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110114) — Breast cancer tissue vs. adjacent normal tissue (13 samples: 10 tumor, 3 normal)

## Live Demo
🔗 [View the dashboard][(https://your-app-url.streamlit.app](https://gffy8kkjy9a5eamdafehcm.streamlit.app/)) <!-- replace with your actual deployed URL -->

## Features

- **Differential Expression**: DESeq2-based analysis (via `pydeseq2`) comparing tumor vs. adjacent normal tissue
- **Interactive Visualizations**:
  - PCA plot — sample clustering by condition
  - Volcano plot — significant genes with adjustable thresholds
  - MA plot — mean expression vs. fold change
  - Heatmap — top DEGs with hierarchical clustering
- **Pathway Enrichment**: GO (Biological Process, Cellular Component, Molecular Function) and KEGG pathway analysis via `gseapy`/Enrichr
- **Gene Explorer**: search any gene to see its expression, differential stats, and pathway membership
- **Adjustable filters**: live p-value and fold-change cutoffs update all plots in real time

## Tech Stack

| Component | Tool |
|---|---|
| Differential expression | `pydeseq2` |
| Gene ID mapping | `mygene` |
| Enrichment analysis | `gseapy` (Enrichr) |
| Dashboard | `streamlit` |
| Visualizations | `plotly` |
| Data source | `GEOparse` (NCBI GEO) |

## Project Structure
analysis.py # Data download, DESeq2, enrichment — run once to generate outputs
├── streamlit_app.py # Interactive dashboard — reads pre-computed results
├── requirements.txt # Dependencies for the dashboard
└── Data/ # Pre-computed analysis outputs (CSV files)
├── GSE110114_DESeq2_results.csv
├── gene_symbol_map.csv
├── GO_enrichment_results.csv
├── KEGG_enrichment_results.csv
├── metadata.csv
└── vst_counts.csv


## Running Locally

**1. Clone the repo**
```bash
git clone https://github.com/veroniperera/RNA-seq-Differential-Gene-Expression-Dashboard-.git
cd RNA-seq-Differential-Gene-Expression-Dashboard-
```

**2. Install dependencies**
```bash
pip install -r requirements.txt
```

**3. Run the dashboard** (uses pre-computed results in `Data/`)
```bash
streamlit run streamlit_app.py
```

**To re-run the full analysis pipeline from scratch** (optional — requires additional packages: `pydeseq2`, `gseapy`, `GEOparse`, `mygene`):
```bash
pip install pydeseq2 gseapy GEOparse mygene
python analysis.py
```

## Key Findings

- **8,577** genes significantly differentially expressed (padj < 0.05, |log2FC| > 1)
- Top enriched pathways are predominantly **immune-related** (MHC complexes, T cell proliferation, interferon-gamma signaling) — consistent with known immune cell infiltration in tumor tissue
- Full results and pathway tables browsable directly in the dashboard's Enrichment tab

## Limitations

- Small control group (n=3 adjacent normal samples) limits statistical power for lowly-expressed genes
- Enrichment results reflect Enrichr's 2021 gene set libraries (GO_2021, KEGG_2021_Human)

## Author

Veroni Perera
