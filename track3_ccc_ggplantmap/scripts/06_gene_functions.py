#!/usr/bin/env python3
"""Gene function annotations for top 4 LASSO biomarker genes + Ca2+ signaling analysis."""
# --- portable paths (de-sandboxed; replaces /mnt/results and /workspace) ---
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
REPO_ROOT = _os.path.abspath(_os.path.join(_HERE, '..', '..'))
RESULTS = _os.environ.get("ASO_ROOT", REPO_ROOT)          # holds tables/ and figures/
ATLAS   = _os.environ.get("ASO_ATLAS", _os.path.join(REPO_ROOT, "atlas"))  # large intermediates (not shipped)
WORK    = _os.environ.get("ASO_WORK", _os.path.join(REPO_ROOT, "work"))    # scratch outputs
_os.makedirs(WORK, exist_ok=True)
# --- end portable paths ---

import csv
import json

# Gene function annotations based on TAIR, UniProt, KEGG, SUBA, and literature
genes = [
    {
        "gene_id": "AT4G10250",
        "alias": "ATHSP22.0",
        "coefficient": "+0.723",
        "direction": "Positive predictor (flight-enriched)",
        "selection_frequency": "100%",
        "protein_name": "HSP20-like chaperones superfamily protein (22.0 kDa heat shock protein)",
        "protein_family": "HSP20/alpha-crystallin family (Pfam PF00011)",
        "go_terms": "BP: response to heat, response to reactive oxygen species, protein folding; MF: unfolded protein binding; CC: endoplasmic reticulum, chloroplast",
        "subcellular_localization": "Endomembrane system (ER), chloroplast stroma",
        "kegg_pathway": "ath04141 (Protein processing in endoplasmic reticulum)",
        "molecular_weight": "22.0 kDa, 195 aa",
        "known_function": "Small heat shock protein (sHSP) that functions as a molecular chaperone, binding unfolded proteins to prevent aggregation. Localized to the endomembrane system and chloroplast. Involved in protein folding, response to heat stress and ROS. The ABI1-HSP22 interaction modulates ABA-auxin signaling crosstalk, providing a link between proteostress and hormone signaling.",
        "spaceflight_relevance": "As a chaperone upregulated in flight samples, HSP22.0 likely reflects proteostatic stress under microgravity. Spaceflight induces protein misfolding via mechanical stress, oxidative stress, and altered calcium homeostasis. sHSPs are first-line defenders preventing protein aggregation. Its positive coefficient indicates spaceflight conditions activate proteostress protective mechanisms. The ABA-auxin crosstalk link suggests hormone signaling redistribution under microgravity."
    },
    {
        "gene_id": "AT3G07365",
        "alias": "NAT (natural antisense transcript)",
        "coefficient": "+0.828",
        "direction": "Positive predictor (flight-enriched, strongest coefficient)",
        "selection_frequency": "100%",
        "protein_name": "Natural antisense transcript (non-coding RNA, NOT protein-coding)",
        "protein_family": "Long non-coding natural antisense transcript (lncNAT)",
        "go_terms": "N/A (non-coding RNA, no protein product)",
        "subcellular_localization": "Nuclear (transcriptional regulatory function)",
        "kegg_pathway": "N/A",
        "molecular_weight": "N/A (non-coding RNA)",
        "known_function": "AT3G07365 is a natural antisense transcript (NAT) that overlaps with AT3G46230 on the opposite DNA strand. NATs are endogenous non-coding RNAs that form double-stranded RNA structures with their sense transcripts. They regulate sense-strand gene expression through transcriptional interference, RNA-RNA pairing, or siRNA-mediated silencing (DCL1/DCL2/RDR6 pathway). Identified as heat-induced in proteomics studies (Babbar et al. 2021, PMC7835529). NATs are particularly associated with environmental stress responses.",
        "spaceflight_relevance": "As the strongest positive predictor of spaceflight (coefficient +0.828), this non-coding RNA suggests that regulatory RNA dynamics, not just protein-coding changes, are central to the spaceflight transcriptional signature. NATs can rapidly modulate gene expression without requiring translation, which may be advantageous for rapid environmental adaptation. The overlap with AT3G46230 (a protein-coding gene) suggests spaceflight may alter the sense-antisense regulatory balance. This finding highlights the importance of including non-coding RNAs in spaceflight biomarker panels."
    },
    {
        "gene_id": "AT2G14247",
        "alias": "Expressed protein",
        "coefficient": "-0.364",
        "direction": "Negative predictor (ground-enriched)",
        "selection_frequency": "100%",
        "protein_name": "Expressed protein (function unknown)",
        "protein_family": "No characterized Pfam domains; 78 aa small protein",
        "go_terms": "BP: biological_process unknown; MF: molecular_function unknown; CC: chloroplast",
        "subcellular_localization": "Chloroplast (SUBAcon consensus: plastid; TargetP: chloroplast Class 2; experimental: mitochondrion PMID 21311031)",
        "kegg_pathway": "Not assigned",
        "molecular_weight": "8.65 kDa, 78 aa",
        "known_function": "Small chloroplast-localized protein with unknown molecular function. Only 78 amino acids. Has broad BLAST hits across bacteria and eukaryotes (35,333 hits in 2,444 species), suggesting a conserved but uncharacterized function. The chloroplast localization suggests a role in photosynthesis or plastid metabolism.",
        "spaceflight_relevance": "As a negative predictor (ground-enriched), this chloroplast protein is suppressed under spaceflight. This is consistent with known spaceflight effects on photosynthesis: microgravity and altered light conditions in orbit disrupt chloroplast function, photosynthetic electron transport, and plastid gene expression. The suppression of a chloroplast protein in flight samples aligns with the broader pattern of photosynthetic apparatus downregulation observed in spaceflight transcriptomics studies."
    },
    {
        "gene_id": "AT4G01390",
        "alias": "TRAF-like / MATH domain protein",
        "coefficient": "-0.819",
        "direction": "Negative predictor (ground-enriched, strongest negative coefficient)",
        "selection_frequency": "100%",
        "protein_name": "Meprin and TRAF homology domain-containing protein / MATH domain-containing protein",
        "protein_family": "MATH/TRAF domain (Pfam PF00917); TNF receptor-associated factor-like",
        "go_terms": "BP: signal transduction, response to stress; CC: cellular_component unknown",
        "subcellular_localization": "Not definitively localized; H2O2-responsive (identified in H2O2-induced proteome)",
        "kegg_pathway": "Not assigned to specific KEGG pathway",
        "molecular_weight": "Predicted ~30-40 kDa",
        "known_function": "TRAF-like protein containing a MATH (meprin and TRAF homology) domain. In animals, TRAF proteins are key signal transducers in TNF receptor pathways, mediating stress and immune signaling. In plants, MATH-domain proteins are involved in abiotic stress responses (Kushwaha et al. 2016, Front Plant Sci). Listed in the Leaf Senescence Database (LSD_871) with signal transduction function. Identified as H2O2-responsive, linking it to oxidative stress signaling.",
        "spaceflight_relevance": "As the strongest negative predictor (coefficient -0.819), this TRAF-like signal transduction protein is strongly enriched in ground controls and suppressed in flight. This suggests that ground-condition stress signaling pathways (particularly oxidative stress and senescence-associated signaling) are attenuated under microgravity. The suppression of a TRAF-domain signal transducer in flight may reflect altered stress perception or dampened stress signaling cascades under microgravity conditions, where mechanical and gravitational stimuli are absent."
    }
]

# Write CSV
with open(RESULTS + '/tables/gene_function_annotations.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=genes[0].keys())
    writer.writeheader()
    writer.writerows(genes)

print(f"Gene function annotations saved: {len(genes)} genes")

# Ca2+ signaling analysis
ca2_analysis = {
    "pathway": "Ca2+ signaling",
    "total_communication_strength": 249.5,
    "rank": 1,
    "comparison": {
        "Ca2+": 249.5,
        "BR": 54.5,
        "flg22": 34.2,
        "ABA": 33.9,
        "chitin": 28.1
    },
    "biological_mechanism": "Ca2+ is a primary second messenger in plant cells, transducing extracellular stimuli into intracellular responses through cytosolic Ca2+ oscillations and gradients. In the context of spaceflight, Ca2+ signaling is particularly relevant because: (1) gravity perception in plants involves statolith sedimentation in specialized cells, which triggers Ca2+ channel activation; (2) microgravity alters cytosolic Ca2+ oscillation patterns and disrupts the normal Ca2+ gradients that guide directional growth; (3) Ca2+ signaling cascades activate downstream effectors including CBL-CIPK complexes (calcineurin B-like proteins interacting with CBL-interiting protein kinases), calmodulin (CaM), and calcium-dependent protein kinases (CDPKs/CPKs), which regulate stress responses, ion transport, and gene expression.",
    "ccc_relevance": "The dominance of Ca2+ signaling (249.5 total communication strength, 4.6x the next pathway BR) in the atlas cell-cell communication network indicates that Ca2+-mediated intercellular signaling is the primary communication channel in Arabidopsis seedling tissues. This is consistent with Ca2+ being a central hub in plant environmental signaling, coordinating responses across cell types. The Ca2+ pathway includes ligand-receptor pairs that trigger Ca2+ influx, which then propagates as intercellular Ca2+ waves through plasmodesmata and apoplastic routes.",
    "connection_to_lasso_genes": "Two of the four top LASSO genes connect to Ca2+ signaling: (1) AT4G10250 (HSP22.0) is induced by Ca2+-dependent stress signaling and participates in ABA-auxin crosstalk, which is downstream of Ca2+ signaling; (2) AT4G01390 (TRAF-like) is a signal transduction protein that may participate in Ca2+-activated stress signaling cascades. The convergence of the LASSO biomarker panel and the dominant CCC pathway on stress signaling is consistent with spaceflight being perceived as a multifaceted stressor.",
    "key_cbl_cipk_components": "The Ca2+ signaling pathway in the CCC database includes ligand-receptor pairs that activate Ca2+ channels and downstream CBL-CIPK, calmodulin, and CDPK cascades. These regulate ion transporters (including K+, NO3-, Mg2+ transport), stress-responsive gene expression, and hormone signaling crosstalk."
}

with open(RESULTS + '/tables/ca2_signaling_analysis.json', 'w') as f:
    json.dump(ca2_analysis, f, indent=2)

print("Ca2+ signaling analysis saved")
print("\n=== Summary ===")
for g in genes:
    print(f"\n{g['gene_id']} ({g['alias']}):")
    print(f"  Coefficient: {g['coefficient']} ({g['direction']})")
    print(f"  Function: {g['protein_name']}")
    print(f"  Spaceflight relevance: {g['spaceflight_relevance'][:120]}...")
