"""
Shared Curated Reaction & Metabolite Mapping Module
=====================================================
Single source of truth for reaction/metabolite identity in CGA009 v1.13+.
All validation scripts MUST import from here instead of hardcoding IDs.

CORRECTED mappings (verified in Phase 1 baseline audit):
- rxn01116_c0 = Ru5P epimerase (NOT RuBisCO)
- rxn05040_c0 = DHBP synthase (NOT RuBisCO)
- rxn02507_c0 = indole-3-glycerol-P synthase (NOT RuBisCO)
- rxn00018_c0 = TRUE RuBisCO carboxylase
- rxn01111_c0 = phosphoribulokinase (PRK)
"""

# ═══════════════════════════════════════════════════════════════════════════════
# CALVIN-BENSON-BASSHAM (CBB) CYCLE
# ═══════════════════════════════════════════════════════════════════════════════
RUBISCO = "rxn00018_c0"          # EC 4.1.1.39: CO2 + RuBP + H2O -> 2x 3PG
PRK = "rxn01111_c0"              # EC 2.7.1.19: ATP + Ru5P -> RuBP + ADP
RU5P_EPIMERASE = "rxn01116_c0"   # EC 5.1.3.1: Ru5P <-> Xu5P (NOT RuBisCO!)

# These are NOT Calvin cycle / NOT RuBisCO:
NOT_RUBISCO = {
    "rxn01116_c0": "Ru5P_epimerase (pentose phosphate, NOT RuBisCO)",
    "rxn05040_c0": "DHBP_synthase (riboflavin biosynthesis, NOT RuBisCO)",
    "rxn02507_c0": "indole-3-glycerol-P_synthase (tryptophan, NOT RuBisCO)",
}

# ═══════════════════════════════════════════════════════════════════════════════
# TCA CYCLE
# ═══════════════════════════════════════════════════════════════════════════════
CITRATE_SYNTHASE = "rxn00256_c0"   # EC 2.3.3.1
ACONITASE = "rxn00974_c0"         # EC 4.2.1.3
IDH = "rxn00199_c0"               # EC 1.1.1.42 isocitrate dehydrogenase
SUCCINYL_COA_SYNTH = "rxn00285_c0" # EC 6.2.1.5
SDH = "rxn00288_c0"               # EC 1.3.5.1 succinate dehydrogenase

# ═══════════════════════════════════════════════════════════════════════════════
# GLYOXYLATE SHUNT
# ═══════════════════════════════════════════════════════════════════════════════
ICL = "rxn00336_c0"    # EC 4.1.3.1 isocitrate lyase
MS = "rxn00330_c0"     # EC 2.3.3.9 malate synthase

# ═══════════════════════════════════════════════════════════════════════════════
# PYRUVATE METABOLISM
# ═══════════════════════════════════════════════════════════════════════════════
POR = "rxn13974_c0"    # EC 1.2.7.1 pyruvate:ferredoxin oxidoreductase

# ═══════════════════════════════════════════════════════════════════════════════
# NITROGEN FIXATION & HYDROGEN
# ═══════════════════════════════════════════════════════════════════════════════
NITROGENASE = "rxn06874_c0"   # EC 1.18.6.1: N2 + 8Fd_red + 16ATP -> 2NH3 + H2
HYDROGENASE = "rxn05759_c0"   # EC 1.12.7.2: 2H+ + 2Fd_red <-> H2 + 2Fd_ox

# ═══════════════════════════════════════════════════════════════════════════════
# PHOTOSYNTHESIS & ENERGY
# ═══════════════════════════════════════════════════════════════════════════════
PHOTO_RC = "rxnTX73PHO001_c0"    # Type II photosynthetic reaction center
BC1_COMPLEX = "rxnTX73PHO002_c0" # Cytochrome bc1 Q-cycle
ATP_SYNTHASE = "rxn08173_c0"     # F-type ATP synthase

# ═══════════════════════════════════════════════════════════════════════════════
# PHB
# ═══════════════════════════════════════════════════════════════════════════════
PHB_SYNTHASE = "rxn15455_c0"

# ═══════════════════════════════════════════════════════════════════════════════
# EXCHANGE REACTIONS
# ═══════════════════════════════════════════════════════════════════════════════
EX_ACETATE = "EX_cpd00029_e0"
EX_SUCCINATE = "EX_cpd00036_e0"
EX_BUTYRATE = "EX_cpd00211_e0"
EX_CO2 = "EX_cpd00011_e0"
EX_H2 = "EX_cpd11640_e0"
EX_PHOTON = "EX_cpd11632_e0"
EX_NH3 = "EX_cpd00013_e0"

ALL_CARBON_EXCHANGES = [EX_ACETATE, EX_SUCCINATE, EX_BUTYRATE]

# ═══════════════════════════════════════════════════════════════════════════════
# BIOMASS
# ═══════════════════════════════════════════════════════════════════════════════
BIOMASS = "bio1"

# ═══════════════════════════════════════════════════════════════════════════════
# KEY METABOLITES (CORRECTED IDENTITIES)
# ═══════════════════════════════════════════════════════════════════════════════
MET_RUBP = "cpd00871"       # D-Ribulose-1,5-bisphosphate (RuBP) — RuBisCO substrate
MET_3PG = "cpd00169"        # 3-Phosphoglycerate — RuBisCO product (NOT RuBP!)
MET_RU5P = "cpd00171"       # D-Ribulose-5-phosphate
MET_XU5P = "cpd00198"       # D-Xylulose-5-phosphate
MET_G3P = "cpd00102"        # Glyceraldehyde-3-phosphate
MET_F6P = "cpd00085"        # D-Fructose-6-phosphate
MET_FBP = "cpd00072"        # D-Fructose-1,6-bisphosphate
MET_S7P = "cpd00238"        # Sedoheptulose-7-phosphate
MET_PYRUVATE = "cpd00020"   # Pyruvate
MET_ACETYL_COA = "cpd00022" # Acetyl-CoA
MET_FD_RED = "cpd11620"     # Reduced ferredoxin
MET_FD_OX = "cpd11621"      # Oxidized ferredoxin
MET_N2 = "cpd00528"         # Dinitrogen
MET_H2 = "cpd11640"         # Hydrogen gas
MET_PHOTON = "cpd11632"     # Photon

# ═══════════════════════════════════════════════════════════════════════════════
# SUBSTRATE SIMULATION CONFIGS
# ═══════════════════════════════════════════════════════════════════════════════
SUBSTRATE_CONFIGS = {
    "acetate": {"exchange": EX_ACETATE, "uptake": 10.0, "c_per_mol": 2},
    "succinate": {"exchange": EX_SUCCINATE, "uptake": 5.0, "c_per_mol": 4},
    "butyrate": {"exchange": EX_BUTYRATE, "uptake": 5.0, "c_per_mol": 4},
}

# Condition-specific enzyme constraints (evidence-supported)
CONDITION_CONSTRAINTS = {
    "succinate": {
        "block": [(ICL, 0), (MS, 0)],  # (rxn_id, upper_bound)
        "evidence": "McKinlay_Harwood_2011_Fig1: zero glyoxylate shunt on succinate",
    },
}
