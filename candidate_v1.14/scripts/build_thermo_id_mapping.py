#!/usr/bin/env python3
"""
Build Thermodynamic ID Mapping for v1.14 TFA Layer
====================================================
Maps ModelSEED metabolite IDs in v1.13 to KEGG/ChEBI/InChI identifiers
needed for eQuilibrator deltaG estimation.

This is a preparatory step — no model edits.
"""

from pathlib import Path
import sys
import warnings
import re

import cobra
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
MODEL_PATH = ROOT / "baseline" / "mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_for_memote.xml"
OUT = ROOT / "candidate_v1.14" / "thermo_prep"

# Priority metabolites for TFA (central carbon, redox, CBB, TCA, PHB, aromatic)
PRIORITY_METABOLITES = {
    # Central carbon
    "cpd00020": "Pyruvate",
    "cpd00022": "Acetyl-CoA",
    "cpd00029": "Acetate",
    "cpd00036": "Succinate",
    "cpd00040": "Glyoxylate",
    "cpd00137": "Citrate",
    "cpd00142": "cis-Aconitate",
    "cpd00024": "2-Oxoglutarate",
    "cpd00106": "Fumarate",
    "cpd00130": "L-Malate",
    "cpd00032": "Oxaloacetate",
    "cpd00061": "Phosphoenolpyruvate",
    # Redox carriers
    "cpd00003": "NAD+",
    "cpd00004": "NADH",
    "cpd00005": "NADPH",
    "cpd00006": "NADP+",
    "cpd00015": "FAD",
    "cpd00982": "FADH2",
    # CBB cycle
    "cpd00101": "Ribose-5-phosphate",
    "cpd00171": "Ribulose-5-phosphate",
    "cpd00169": "Ribulose-1,5-bisphosphate",
    "cpd00102": "Glyceraldehyde-3-phosphate",
    "cpd00085": "Fructose-6-phosphate",
    "cpd00031": "Glycerate-3-phosphate",
    # Energy
    "cpd00002": "ATP",
    "cpd00008": "ADP",
    "cpd00018": "AMP",
    "cpd00009": "Phosphate",
    "cpd00012": "PPi",
    # PHB
    "cpd00211": "Butyrate",
    "cpd00196": "Acetoacetyl-CoA",
    "cpd00557": "3-Hydroxybutyryl-CoA",
    # CO2/HCO3
    "cpd00011": "CO2",
    "cpd00242": "HCO3-",
    # H2
    "cpd11640": "H2",
    # N2
    "cpd00528": "N2",
    "cpd00013": "NH3",
    # Ferredoxin
    "cpd11621": "Oxidized_ferredoxin",
    "cpd11620": "Reduced_ferredoxin",
    # Water/H+
    "cpd00001": "H2O",
    "cpd00067": "H+",
}

# Known ModelSEED -> KEGG mappings for common metabolites
MODELSEED_TO_KEGG = {
    "cpd00020": "C00022",  # Pyruvate
    "cpd00022": "C00024",  # Acetyl-CoA
    "cpd00029": "C00033",  # Acetate
    "cpd00036": "C00042",  # Succinate
    "cpd00040": "C00048",  # Glyoxylate
    "cpd00137": "C00158",  # Citrate
    "cpd00024": "C00026",  # 2-Oxoglutarate
    "cpd00106": "C00122",  # Fumarate
    "cpd00130": "C00149",  # L-Malate
    "cpd00032": "C00036",  # Oxaloacetate
    "cpd00061": "C00074",  # PEP
    "cpd00003": "C00003",  # NAD+
    "cpd00004": "C00004",  # NADH
    "cpd00005": "C00005",  # NADPH
    "cpd00006": "C00006",  # NADP+
    "cpd00015": "C00016",  # FAD
    "cpd00002": "C00002",  # ATP
    "cpd00008": "C00008",  # ADP
    "cpd00009": "C00009",  # Pi
    "cpd00011": "C00011",  # CO2
    "cpd00242": "C00288",  # HCO3-
    "cpd00085": "C00085",  # F6P
    "cpd00102": "C00118",  # G3P
    "cpd00031": "C00197",  # 3PG
    "cpd00211": "C00246",  # Butyrate
    "cpd00001": "C00001",  # H2O
    "cpd00067": "C00080",  # H+
    "cpd00013": "C00014",  # NH3
    "cpd00101": "C00117",  # R5P
    "cpd00171": "C00199",  # Ru5P
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    model = cobra.io.read_sbml_model(str(MODEL_PATH))

    # Build mapping table
    rows = []
    for met in model.metabolites:
        # Extract base cpd ID (without compartment)
        base_id = met.id.split("[")[0] if "[" in met.id else met.id.rstrip("_c0_e0")
        # Try to match to priority list
        match = None
        for cpd_id in PRIORITY_METABOLITES:
            if cpd_id in met.id:
                match = cpd_id
                break

        if match:
            kegg_id = MODELSEED_TO_KEGG.get(match, "NEEDS_LOOKUP")
            rows.append({
                "model_metabolite_id": met.id,
                "modelseed_base_id": match,
                "name": met.name or PRIORITY_METABOLITES[match],
                "formula": met.formula,
                "charge": met.charge,
                "compartment": met.compartment,
                "kegg_id": kegg_id,
                "chebi_id": "NEEDS_LOOKUP",
                "equilibrator_id": kegg_id if kegg_id != "NEEDS_LOOKUP" else "NEEDS_LOOKUP",
                "priority": "high" if match in list(PRIORITY_METABOLITES.keys())[:20] else "medium",
                "notes": "",
            })

    df = pd.DataFrame(rows)
    df = df.sort_values(["priority", "modelseed_base_id"])
    df.to_csv(OUT / "thermo_id_mapping.tsv", sep="\t", index=False)

    # Summary
    mapped = (df["kegg_id"] != "NEEDS_LOOKUP").sum()
    total = len(df)
    print(f"Thermodynamic ID mapping: {mapped}/{total} metabolites have KEGG IDs")
    print(f"Output: {OUT / 'thermo_id_mapping.tsv'}")

    # Also identify reactions in central carbon that need deltaG
    priority_rxns = [
        "rxn00336_c0", "rxn00330_c0", "rxn00256_c0", "rxn00974_c0",
        "rxn00199_c0", "rxn00285_c0", "rxn00288_c0", "rxn13974_c0",
        "rxn01116_c0", "rxn06874_c0",
    ]

    rxn_rows = []
    for rid in priority_rxns:
        if rid in model.reactions:
            rxn = model.reactions.get_by_id(rid)
            reactants = "; ".join([f"{v} {m.id}" for m, v in rxn.metabolites.items() if v < 0])
            products = "; ".join([f"{v} {m.id}" for m, v in rxn.metabolites.items() if v > 0])
            rxn_rows.append({
                "reaction_id": rid,
                "name": rxn.name,
                "equation_summary": rxn.reaction,
                "reversible": rxn.lower_bound < 0,
                "current_lb": rxn.lower_bound,
                "current_ub": rxn.upper_bound,
                "deltaG_source": "NEEDS_eQuilibrator",
                "tfa_priority": "high",
            })

    rxn_df = pd.DataFrame(rxn_rows)
    rxn_df.to_csv(OUT / "tfa_priority_reactions.tsv", sep="\t", index=False)
    print(f"TFA priority reactions: {len(rxn_df)} reactions identified")

    return 0


if __name__ == "__main__":
    sys.exit(main())
