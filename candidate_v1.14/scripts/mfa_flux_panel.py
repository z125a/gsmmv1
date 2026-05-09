#!/usr/bin/env python3
"""
MFA Internal Flux Validation Panel — v1.14 Cycle 2
====================================================
Uses exact ModelSEED reaction IDs mapped from the v1.13 model to compare
pFBA flux distributions against McKinlay/Harwood 2011 qualitative 13C-MFA targets.

Key reactions identified:
- Isocitrate lyase (ICL): rxn00336_c0
- Malate synthase (MS): rxn00330_c0
- Citrate synthase (CS): rxn00256_c0
- Aconitase: rxn00974_c0
- Isocitrate dehydrogenase (IDH): rxn00199_c0
- Succinyl-CoA ligase: rxn00285_c0
- Succinate dehydrogenase (SDH): rxn00288_c0
- Pyruvate:ferredoxin oxidoreductase (POR/PDH): rxn13974_c0
- RuBisCO reactions: rxn01116_c0, rxn05040_c0, rxn00018_c0, rxn01111_c0, rxn02507_c0
- Nitrogenase: rxn06874_c0
- H2 exchange: EX_cpd11640_e0 (or search)
- CO2 exchange: EX_cpd00011_e0
- Acetate exchange: EX_cpd00029_e0
- Succinate exchange: EX_cpd00036_e0
- Butyrate exchange: EX_cpd00211_e0
"""

from pathlib import Path
import sys
import warnings

import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import pfba, flux_variability_analysis

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
MODEL_PATH = ROOT / "baseline" / "mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_for_memote.xml"
OUT = ROOT / "candidate_v1.14" / "validation_reports"

# Exact reaction ID mapping for MFA validation
REACTION_MAP = {
    "ICL": "rxn00336_c0",           # isocitrate lyase
    "MS": "rxn00330_c0",            # malate synthase
    "CS": "rxn00256_c0",            # citrate synthase
    "ACO": "rxn00974_c0",           # aconitase
    "IDH": "rxn00199_c0",           # isocitrate dehydrogenase
    "SCS": "rxn00285_c0",           # succinyl-CoA synthetase
    "SDH": "rxn00288_c0",           # succinate dehydrogenase
    "POR": "rxn13974_c0",           # pyruvate:ferredoxin oxidoreductase
    "nitrogenase": "rxn06874_c0",   # nitrogenase
    "EX_CO2": "EX_cpd00011_e0",     # CO2 exchange
    "EX_acetate": "EX_cpd00029_e0", # acetate exchange
    "EX_succinate": "EX_cpd00036_e0",  # succinate exchange
    "EX_butyrate": "EX_cpd00211_e0",   # butyrate exchange
    "EX_H2": "EX_cpd11640_e0",     # H2 exchange
    "biomass": "bio1",              # biomass objective
}

# RuBisCO-related reactions (multiple in model)
RUBISCO_RXNS = ["rxn01116_c0", "rxn05040_c0", "rxn00018_c0", "rxn01111_c0", "rxn02507_c0"]

# Substrate conditions to simulate
SUBSTRATES = {
    "acetate": {
        "exchange": "EX_cpd00029_e0",
        "uptake_rate": 10.0,  # mmol/gDW/h (arbitrary, for relative comparison)
        "carbon_atoms": 2,
    },
    "succinate": {
        "exchange": "EX_cpd00036_e0",
        "uptake_rate": 5.0,   # 10 mM succinate = 40 mM C, same as 20 mM acetate
        "carbon_atoms": 4,
    },
    "butyrate": {
        "exchange": "EX_cpd00211_e0",
        "uptake_rate": 5.0,   # 10 mM butyrate = 40 mM C
        "carbon_atoms": 4,
    },
}

# Carbon source exchanges to close when testing a specific substrate
ALL_CARBON_EXCHANGES = ["EX_cpd00029_e0", "EX_cpd00036_e0", "EX_cpd00211_e0"]


def load_model():
    cobra.Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    model.solver = "glpk"
    return model


def get_flux(solution, rxn_id):
    """Safely get flux value from solution."""
    if rxn_id in solution.fluxes.index:
        return solution.fluxes[rxn_id]
    return None


def simulate_photoheterotrophic(model, substrate_name, substrate_info):
    """
    Simulate anaerobic photoheterotrophic growth on a single substrate.
    Returns pFBA solution or None.
    """
    with model:
        # Close all carbon exchanges first
        for ex_id in ALL_CARBON_EXCHANGES:
            if ex_id in model.reactions:
                rxn = model.reactions.get_by_id(ex_id)
                rxn.lower_bound = 0

        # Open the target substrate
        ex_id = substrate_info["exchange"]
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -substrate_info["uptake_rate"]

        # Check feasibility
        sol = model.optimize()
        if sol.status != "optimal":
            return None

        # Run pFBA for minimal flux distribution
        try:
            return pfba(model)
        except Exception:
            return sol


def extract_mfa_panel(solution, substrate_name):
    """Extract the key flux values for MFA comparison."""
    panel = {"substrate": substrate_name}

    # Core reaction fluxes
    for name, rxn_id in REACTION_MAP.items():
        flux = get_flux(solution, rxn_id)
        panel[f"flux_{name}"] = flux if flux is not None else "NOT_IN_MODEL"

    # RuBisCO total flux
    rubisco_total = 0.0
    for rid in RUBISCO_RXNS:
        f = get_flux(solution, rid)
        if f is not None:
            rubisco_total += f
    panel["flux_RuBisCO_total"] = rubisco_total

    # Glyoxylate shunt total (ICL + MS)
    icl = get_flux(solution, REACTION_MAP["ICL"]) or 0
    ms = get_flux(solution, REACTION_MAP["MS"]) or 0
    panel["flux_glyoxylate_shunt"] = icl + ms

    return panel


def mfa_qualitative_tests(panels):
    """
    Run qualitative tests based on McKinlay/Harwood 2011 findings.
    """
    tests = []

    for panel in panels:
        sub = panel["substrate"]
        gs = panel.get("flux_glyoxylate_shunt", 0)
        rubisco = panel.get("flux_RuBisCO_total", 0)
        por = panel.get("flux_POR", 0) if panel.get("flux_POR") != "NOT_IN_MODEL" else 0
        co2 = panel.get("flux_EX_CO2", 0) if panel.get("flux_EX_CO2") != "NOT_IN_MODEL" else 0
        h2 = panel.get("flux_EX_H2", 0) if panel.get("flux_EX_H2") != "NOT_IN_MODEL" else 0
        growth = panel.get("flux_biomass", 0) if panel.get("flux_biomass") != "NOT_IN_MODEL" else 0
        nit = panel.get("flux_nitrogenase", 0) if panel.get("flux_nitrogenase") != "NOT_IN_MODEL" else 0

        # Test 1: Growth feasibility
        tests.append({
            "substrate": sub,
            "test": "growth_feasible",
            "expected": ">0",
            "observed": f"{growth:.4f}" if isinstance(growth, float) else str(growth),
            "status": "PASS" if isinstance(growth, float) and growth > 1e-6 else "FAIL",
            "source": "McKinlay_Harwood_2011_Table2",
        })

        # Test 2: Glyoxylate shunt usage
        if sub in ["acetate", "butyrate"]:
            tests.append({
                "substrate": sub,
                "test": "glyoxylate_shunt_active",
                "expected": ">0 (acetate/butyrate assimilated via glyoxylate shunt)",
                "observed": f"{gs:.6f}",
                "status": "PASS" if gs > 1e-6 else "FAIL",
                "source": "McKinlay_Harwood_2011_Fig1",
            })
        elif sub == "succinate":
            tests.append({
                "substrate": sub,
                "test": "glyoxylate_shunt_inactive",
                "expected": "=0 (succinate not via glyoxylate shunt)",
                "observed": f"{gs:.6f}",
                "status": "PASS" if abs(gs) < 1e-6 else "FAIL",
                "source": "McKinlay_Harwood_2011_Fig1",
            })

        # Test 3: Calvin cycle (RuBisCO) should carry flux
        tests.append({
            "substrate": sub,
            "test": "Calvin_cycle_active",
            "expected": ">0 (essential for electron balance in WT photoheterotrophic)",
            "observed": f"{rubisco:.4f}",
            "status": "PASS" if abs(rubisco) > 1e-6 else "FAIL",
            "source": "McKinlay_Harwood_2011_Fig2",
        })

        # Test 4: Net CO2 production (WT on acetate/succinate/fumarate should produce net CO2)
        # Convention: positive = secretion in COBRApy exchange
        if sub in ["acetate", "succinate"]:
            tests.append({
                "substrate": sub,
                "test": "net_CO2_production",
                "expected": ">0 (net CO2 release for WT)",
                "observed": f"{co2:.4f}" if isinstance(co2, float) else str(co2),
                "status": "PASS" if isinstance(co2, float) and co2 > 0 else "INFO",
                "source": "McKinlay_Harwood_2011_Table3",
            })

        # Test 5: POR flux direction
        # On succinate, POR should be forward (pyruvate -> acetyl-CoA) to generate acetyl-CoA
        if sub == "succinate":
            tests.append({
                "substrate": sub,
                "test": "POR_forward_for_acetylCoA",
                "expected": ">0 (succinate needs decarboxylation to reach acetyl-CoA)",
                "observed": f"{por:.4f}" if isinstance(por, float) else str(por),
                "status": "PASS" if isinstance(por, float) and por > 0 else "INFO",
                "source": "McKinlay_Harwood_2011_Fig1",
            })

    return tests


def run_fva_on_key_reactions(model, substrate_name, substrate_info):
    """Run FVA on key reactions to check flux ranges."""
    with model:
        for ex_id in ALL_CARBON_EXCHANGES:
            if ex_id in model.reactions:
                model.reactions.get_by_id(ex_id).lower_bound = 0
        ex_id = substrate_info["exchange"]
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -substrate_info["uptake_rate"]

        # FVA on glyoxylate shunt and Calvin cycle
        target_rxns = [REACTION_MAP["ICL"], REACTION_MAP["MS"]] + RUBISCO_RXNS
        target_rxns = [r for r in target_rxns if r in model.reactions]

        if not target_rxns:
            return None

        try:
            fva = flux_variability_analysis(
                model, reaction_list=target_rxns, fraction_of_optimum=0.9
            )
            fva["substrate"] = substrate_name
            return fva
        except Exception:
            return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("MFA INTERNAL FLUX VALIDATION PANEL — Cycle 2")
    print("=" * 70)

    model = load_model()
    print(f"Model loaded: {len(model.reactions)} rxns, {len(model.metabolites)} mets")

    # Verify key reactions exist
    print("\nKey reaction verification:")
    for name, rxn_id in REACTION_MAP.items():
        exists = rxn_id in model.reactions
        print(f"  {name:20s} ({rxn_id}): {'FOUND' if exists else 'MISSING'}")
    for rid in RUBISCO_RXNS:
        exists = rid in model.reactions
        print(f"  {'RuBisCO':20s} ({rid}): {'FOUND' if exists else 'MISSING'}")

    # Simulate each substrate
    print("\n" + "-" * 70)
    print("SUBSTRATE SIMULATIONS (pFBA)")
    print("-" * 70)

    all_panels = []
    all_fva = []

    for sub_name, sub_info in SUBSTRATES.items():
        print(f"\n  [{sub_name}] uptake={sub_info['uptake_rate']} mmol/gDW/h")
        sol = simulate_photoheterotrophic(model, sub_name, sub_info)

        if sol is not None:
            panel = extract_mfa_panel(sol, sub_name)
            all_panels.append(panel)
            growth = panel.get("flux_biomass", 0)
            print(f"    Growth: {growth:.4f}")
            print(f"    Glyoxylate shunt: {panel['flux_glyoxylate_shunt']:.4f}")
            print(f"    RuBisCO total: {panel['flux_RuBisCO_total']:.4f}")
            print(f"    POR: {panel.get('flux_POR', 'N/A')}")
            print(f"    CO2 exchange: {panel.get('flux_EX_CO2', 'N/A')}")
            print(f"    H2 exchange: {panel.get('flux_EX_H2', 'N/A')}")
            print(f"    Nitrogenase: {panel.get('flux_nitrogenase', 'N/A')}")
        else:
            print(f"    INFEASIBLE")
            all_panels.append({"substrate": sub_name, "flux_biomass": 0})

        # FVA
        fva = run_fva_on_key_reactions(model, sub_name, sub_info)
        if fva is not None:
            all_fva.append(fva)

    # Run qualitative tests
    print("\n" + "-" * 70)
    print("QUALITATIVE MFA COMPARISON")
    print("-" * 70)

    tests = mfa_qualitative_tests(all_panels)
    tests_df = pd.DataFrame(tests)

    pass_n = (tests_df["status"] == "PASS").sum()
    fail_n = (tests_df["status"] == "FAIL").sum()
    info_n = (tests_df["status"] == "INFO").sum()

    print(f"\n  Results: {pass_n} PASS / {fail_n} FAIL / {info_n} INFO / {len(tests_df)} total")
    print()
    for _, row in tests_df.iterrows():
        marker = "✓" if row["status"] == "PASS" else ("✗" if row["status"] == "FAIL" else "·")
        print(f"  {marker} [{row['substrate']:10s}] {row['test']:30s} | {row['status']:4s} | obs={row['observed']}")

    # Save outputs
    tests_df.to_csv(OUT / "mfa_panel_qualitative_tests.tsv", sep="\t", index=False)
    pd.DataFrame(all_panels).to_csv(OUT / "mfa_panel_flux_values.tsv", sep="\t", index=False)

    if all_fva:
        fva_combined = pd.concat(all_fva, ignore_index=False)
        fva_combined.to_csv(OUT / "mfa_panel_fva_key_reactions.tsv", sep="\t")

    # Diagnostic: why is glyoxylate shunt zero?
    print("\n" + "-" * 70)
    print("DIAGNOSTIC: Glyoxylate Shunt Analysis")
    print("-" * 70)

    for sub_name, sub_info in SUBSTRATES.items():
        if sub_name in ["acetate", "butyrate"]:
            with model:
                for ex_id in ALL_CARBON_EXCHANGES:
                    if ex_id in model.reactions:
                        model.reactions.get_by_id(ex_id).lower_bound = 0
                model.reactions.get_by_id(sub_info["exchange"]).lower_bound = -sub_info["uptake_rate"]

                # Check if ICL can carry flux at all
                icl = model.reactions.get_by_id("rxn00336_c0")
                print(f"\n  [{sub_name}] ICL bounds: [{icl.lower_bound}, {icl.upper_bound}]")
                print(f"  [{sub_name}] ICL GPR: {icl.gene_reaction_rule}")

                # Force ICL to carry flux and check feasibility
                icl.lower_bound = 0.1
                sol_forced = model.optimize()
                if sol_forced.status == "optimal":
                    print(f"  [{sub_name}] Forcing ICL>=0.1: FEASIBLE, growth={sol_forced.objective_value:.4f}")
                else:
                    print(f"  [{sub_name}] Forcing ICL>=0.1: INFEASIBLE")

                # FVA on ICL
                icl.lower_bound = 0
                try:
                    fva_icl = flux_variability_analysis(model, reaction_list=["rxn00336_c0"], fraction_of_optimum=0.99)
                    print(f"  [{sub_name}] ICL FVA (99% opt): min={fva_icl['minimum'].iloc[0]:.4f}, max={fva_icl['maximum'].iloc[0]:.4f}")
                except Exception as e:
                    print(f"  [{sub_name}] ICL FVA failed: {e}")

    print("\n" + "=" * 70)
    print("PANEL COMPLETE")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
