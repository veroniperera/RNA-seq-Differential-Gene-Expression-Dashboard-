#=================================================
#VISUALIZATION
#=================================================
from ast import If
import os 
import pickle as pkl 
import pandas as pd 
import numpy as np
import streamlit as st 
import plotly.express as px 
import plotly.graph_objects as go 
from scipy.cluster.hierarchy  import linkage, dendrogram 

# ------------------------
# CONFIG
# ------------------------
OUTPUT_PATH = r"C:\Users\HP\Downloads\OneDrive\Portfolio\PROJECT 1"

st.set_page_config(
    page_title="RNA-seq DGE Dashboard | GSE110114",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ------------------------
# CUSTOM CSS
# ------------------------
st.markdown("""
<style>
    .main { background-color: #0e1117; }
    h1, h2, h3 { color: #f5f5f5; font-family: 'Helvetica Neue', sans-serif; }
    .stMetric { background-color: #1c1f26; border-radius: 12px; padding: 10px; border: 1px solid #2a2e37; }
    div[data-testid="stMetricValue"] { color: #00d4ff; font-weight: 700; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] {
        background-color: #1c1f26; border-radius: 8px 8px 0 0; padding: 10px 20px;
        color: #ccc; font-weight: 600;
    }
    .stTabs [aria-selected="true"] { background-color: #00d4ff; color: #0e1117; }
    section[data-testid="stSidebar"] { background-color: #12141a; }
</style>
""", unsafe_allow_html=True)

PLOTLY_TEMPLATE = "plotly_dark"
ACCENT = "#00d4ff"
UP_COLOR = "#ff4b4b"
DOWN_COLOR = "#4b9fff"
NEUTRAL_COLOR = "#4a4e57"

# ------------------------
# LOAD DATA
# ------------------------
@st.cache_data
def load_data():
    results_df = pd.read_csv(os.path.join(OUTPUT_PATH, "GSE110114_DESeq2_results.csv"), index_col=0)
    gene_map = pd.read_csv(os.path.join(OUTPUT_PATH, "gene_symbol_map.csv"), index_col=0)
    go_df = pd.read_csv(os.path.join(OUTPUT_PATH, "GO_enrichment_results.csv"))
    kegg_df = pd.read_csv(os.path.join(OUTPUT_PATH, "KEGG_enrichment_results.csv"))
    vst_counts = pd.read_csv(os.path.join(OUTPUT_PATH, "vst_counts.csv"), index_col=0)
    metadata = pd.read_csv(os.path.join(OUTPUT_PATH, "metadata.csv"), index_col=0)
    return results_df, gene_map, go_df, kegg_df, vst_counts, metadata

results_df, gene_map, go_df, kegg_df, vst_counts, metadata = load_data()

# Merge symbols into results
results_df = results_df.merge(
    gene_map[["symbol"]], left_index=True, right_index=True, how="left"
)

# ------------------------
# SIDEBAR CONTROLS
# ------------------------
st.sidebar.title("🧬 Controls")
padj_cutoff = st.sidebar.slider("Adjusted p-value cutoff", 0.001, 0.10, 0.05, 0.001)
lfc_cutoff = st.sidebar.slider("|log2 Fold Change| cutoff", 0.0, 5.0, 1.0, 0.1)
top_n_heatmap = st.sidebar.slider("Top N genes for heatmap", 10, 100, 50, 5)
gene_query = st.sidebar.text_input("🔍 Search a gene symbol", "")

results_df["significant"] = np.where(
    (results_df["padj"] < padj_cutoff) & (results_df["log2FoldChange"].abs() > lfc_cutoff),
    np.where(results_df["log2FoldChange"] > 0, "Up in Tumor", "Up in Normal"),
    "Not significant",
)

sig_df = results_df[results_df["significant"] != "Not significant"]

# ------------------------
# HEADER / KPIs
# ------------------------
st.title("RNA-seq Differential Expression Dashboard")
st.caption("GSE110114 — Breast cancer tissue vs. adjacent normal tissue")

k1, k2, k3, k4 = st.columns(4)
k1.metric("Total genes tested", f"{results_df.shape[0]:,}")
k2.metric("Significant DEGs", f"{sig_df.shape[0]:,}")
k3.metric("Up in Tumor", f"{(sig_df['significant']=='Up in Tumor').sum():,}")
k4.metric("Up in Normal", f"{(sig_df['significant']=='Up in Normal').sum():,}")

st.divider()

# ------------------------
# TABS
# ------------------------
tab_pca, tab_volcano, tab_ma, tab_heatmap, tab_enrich, tab_gene = st.tabs(
    ["📊 PCA", "🌋 Volcano", "📈 MA Plot", "🔥 Heatmap", "🧭 Enrichment", "🔎 Gene Explorer"]
)

# --- PCA ---
with tab_pca:
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    coords = pca.fit_transform(vst_counts.T)
    pca_df = pd.DataFrame(coords, columns=["PC1", "PC2"], index=vst_counts.columns)
    pca_df = pca_df.merge(metadata, left_index=True, right_index=True)

    fig_pca = px.scatter(
        pca_df, x="PC1", y="PC2", color="condition", text=pca_df.index,
        template=PLOTLY_TEMPLATE, height=550,
        color_discrete_sequence=[ACCENT, "#ff6b6b"],
        title=f"PCA — PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%) vs "
              f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)",
    )
    fig_pca.update_traces(marker=dict(size=14, line=dict(width=1, color="white")), textposition="top center")
    st.plotly_chart(fig_pca, use_container_width=True)

# --- VOLCANO ---
with tab_volcano:
    color_map = {"Up in Tumor": UP_COLOR, "Up in Normal": DOWN_COLOR, "Not significant": NEUTRAL_COLOR}
    fig_volcano = px.scatter(
        results_df, x="log2FoldChange", y=-np.log10(results_df["padj"].clip(lower=1e-300)),
        color="significant", color_discrete_map=color_map,
        hover_name="symbol", opacity=0.7, template=PLOTLY_TEMPLATE, height=600,
        labels={"y": "-log10(adjusted p-value)"},
        title="Volcano Plot",
    )
    fig_volcano.add_vline(x=lfc_cutoff, line_dash="dash", line_color="gray")
    fig_volcano.add_vline(x=-lfc_cutoff, line_dash="dash", line_color="gray")
    fig_volcano.add_hline(y=-np.log10(padj_cutoff), line_dash="dash", line_color="gray")

    top_labels = sig_df.reindex(sig_df["padj"].sort_values().index).head(15)
    for _, row in top_labels.iterrows():
        fig_volcano.add_annotation(
            x=row["log2FoldChange"], y=-np.log10(max(row["padj"], 1e-300)),
            text=row["symbol"], showarrow=True, arrowhead=1, font=dict(size=10, color="white"),
        )
    st.plotly_chart(fig_volcano, use_container_width=True)

# --- MA PLOT ---
with tab_ma:
    fig_ma = px.scatter(
        results_df, x=np.log10(results_df["baseMean"].clip(lower=1)), y="log2FoldChange",
        color="significant", color_discrete_map=color_map, hover_name="symbol",
        opacity=0.6, template=PLOTLY_TEMPLATE, height=600,
        labels={"x": "log10(mean expression)"}, title="MA Plot",
    )
    fig_ma.add_hline(y=0, line_color="white", line_width=1)
    st.plotly_chart(fig_ma, use_container_width=True)

# --- HEATMAP ---
with tab_heatmap:
    top_genes = sig_df.reindex(sig_df["padj"].sort_values().index).head(top_n_heatmap).index
    heat_data = vst_counts.loc[vst_counts.index.intersection(top_genes)]
    heat_z = heat_data.sub(heat_data.mean(axis=1), axis=0).div(heat_data.std(axis=1), axis=0)

    row_link = linkage(heat_z.fillna(0), method="ward")
    row_order = dendrogram(row_link, no_plot=True)["leaves"]
    heat_z = heat_z.iloc[row_order]

    labels = gene_map.reindex(heat_z.index)["symbol"].fillna(heat_z.index.to_series())

    fig_heat = go.Figure(data=go.Heatmap(
        z=heat_z.values, x=heat_z.columns, y=labels,
        colorscale="RdBu_r", zmid=0, colorbar=dict(title="z-score"),
    ))
    fig_heat.update_layout(template=PLOTLY_TEMPLATE, height=800, title=f"Top {top_n_heatmap} DEGs (z-scored, clustered)")
    st.plotly_chart(fig_heat, use_container_width=True)

# --- ENRICHMENT ---
with tab_enrich:
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("GO Biological Process")
        go_top = go_df.sort_values("Adjusted P-value").head(15)
        fig_go = px.bar(
            go_top, x=-np.log10(go_top["Adjusted P-value"]), y="Term", orientation="h",
            template=PLOTLY_TEMPLATE, color=-np.log10(go_top["Adjusted P-value"]),
            color_continuous_scale="Viridis", labels={"x": "-log10(adj p-value)"},
        )
        fig_go.update_layout(yaxis=dict(autorange="reversed"), height=550, showlegend=False)
        st.plotly_chart(fig_go, use_container_width=True)

    with col2:
        st.subheader("KEGG Pathways")
        kegg_top = kegg_df.sort_values("Adjusted P-value").head(15)
        fig_kegg = px.bar(
            kegg_top, x=-np.log10(kegg_top["Adjusted P-value"]), y="Term", orientation="h",
            template=PLOTLY_TEMPLATE, color=-np.log10(kegg_top["Adjusted P-value"]),
            color_continuous_scale="Plasma", labels={"x": "-log10(adj p-value)"},
        )
        fig_kegg.update_layout(yaxis=dict(autorange="reversed"), height=550, showlegend=False)
        st.plotly_chart(fig_kegg, use_container_width=True)

# --- GENE EXPLORER ---
with tab_gene:
    if gene_query:
        match = results_df[results_df["symbol"].str.upper() == gene_query.upper()]
        if match.empty:
            st.warning(f"No gene found matching '{gene_query}'")
        else:
            gene_id = match.index[0]
            row = match.iloc[0]

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("log2FC", f"{row['log2FoldChange']:.2f}")
            c2.metric("padj", f"{row['padj']:.2e}")
            c3.metric("baseMean", f"{row['baseMean']:.1f}")
            c4.metric("Status", row["significant"])

            expr = vst_counts.loc[gene_id].to_frame("expression").merge(metadata, left_index=True, right_index=True)
            fig_gene = px.box(
                expr, x="condition", y="expression", color="condition", points="all",
                template=PLOTLY_TEMPLATE, height=500,
                color_discrete_sequence=[ACCENT, "#ff6b6b"],
                title=f"{gene_query} expression (VST) by condition",
            )
            st.plotly_chart(fig_gene, use_container_width=True)

            st.subheader("Pathways containing this gene")
            go_hits = go_df[go_df["Genes"].str.contains(gene_query, case=False, na=False)]
            kegg_hits = kegg_df[kegg_df["Genes"].str.contains(gene_query, case=False, na=False)]
            st.dataframe(pd.concat([go_hits[["Term", "Adjusted P-value"]], kegg_hits[["Term", "Adjusted P-value"]]]))
    else:
        st.info("Enter a gene symbol in the sidebar to explore its expression and pathway membership.")

    

