#!/usr/bin/env python3
"""Regenerate the circuit base tables with CORRECT AGIs, matching original formats.
Outputs to work/clean_tables/ for distribution."""
import os, numpy as np, pandas as pd, anndata as ad

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
RESULTS = os.environ.get("ASO_ROOT", REPO)
ATLAS = os.environ.get("ASO_ATLAS", os.path.join(REPO, "atlas"))
OUT = os.environ.get("ASO_WORK", os.path.join(REPO, "work"))
OUT = os.path.join(OUT, "clean_tables"); os.makedirs(OUT, exist_ok=True)
CTS = ['Epidermal','Mesophyll','Stele','Meristematic','Guard']

gene_order = ['CBL1','CBL9','CBL2','CIPK23','CIPK1','CIPK9','CML38','CML42','CML24','CaM7','CaM3',
              'CDPK3','CDPK4','CAX3','NRT2.1','SOS1','ANNAT1','AKT1','AKT2','GORK','KAT1','KAT2','KC1',
              'HSP22.0','TRAF-like']
gmap = pd.read_csv(os.path.join(REPO,"tables","corrected_gene_agi_map.csv"))
name2agi = dict(zip(gmap.gene, gmap.correct_AGI))
name2agi.update({'HSP22.0':'AT4G10250','TRAF-like':'AT4G01390'})

adata = ad.read_h5ad(os.path.join(ATLAS,"seedling_6d_anndata.h5ad"))
var = set(map(str, adata.var_names))
genes = [g for g in gene_order if name2agi.get(g) in var]
agis = [name2agi[g] for g in genes]
X = adata[:, agis].X
X = np.asarray(X.todense()) if hasattr(X,"todense") else np.asarray(X)
expr = pd.DataFrame(X, columns=genes)
ct = adata.obs['CellType'].values

# 1) expr per celltype (gene=AGI, alias, celltype, mean_expr, pct_expr)
rows=[]
for g,a in zip(genes,agis):
    col=expr[g].values
    for c in CTS:
        m=ct==c
        rows.append({'gene':a,'alias':g,'celltype':c,'mean_expr':float(col[m].mean()),'pct_expr':float((col[m]>0).mean()*100)})
    rows.append({'gene':a,'alias':g,'celltype':'ALL','mean_expr':float(col.mean()),'pct_expr':float((col>0).mean()*100)})
pd.DataFrame(rows).to_csv(os.path.join(OUT,"ca2_k_expr_per_celltype.csv"),index=False)

# 2) global co-expression (gene-name index/cols)
g=expr.corr(method='pearson').reindex(index=genes,columns=genes)
g.to_csv(os.path.join(OUT,"ca2_k_coexpression_matrix.csv"))

# 3) per-cell-type co-expression, row index = "GENE_<ct>", cols = gene names (match original)
for c in CTS:
    sub=expr[ct==c]
    corr=sub[genes].corr(method='pearson').reindex(index=genes,columns=genes)
    corr.index=[f"{gn}_{c}" for gn in corr.index]
    corr.to_csv(os.path.join(OUT,f"ca2_k_coexpression_{c}.csv"))

print("wrote:", ", ".join(sorted(os.listdir(OUT))))
print(f"({len(genes)} genes, {adata.shape[0]} cells)")
