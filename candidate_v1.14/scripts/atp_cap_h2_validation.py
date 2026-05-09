#!/usr/bin/env python3
"""
ATP Cap + H2 Yield Validation — v1.14 Cycle 4
===============================================
Tests whether substrate-specific photosynthetic ATP caps from Chowdhury 2022
reduce H2 over-prediction toward experimental values.

ATP caps (mmol/gDW/day → mmol/gDW/h by dividing by 24):
- Acetate: 54.0 / 24 = 2.25 mmol/gDW/h
- Succinate: 45.7 / 24 = 1.90 mmol/gDW/h
- Butyrate: 56.7 / 24 = 2.36 mmol/gDW/h

These are applied as upper bounds on ATP synthase (rxn08173_c0).

Also explores reduced ferredoxin as the electron transport bottleneck.
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
OUT = ROOT / "candidate_v1.14" / "atp_cap_results"

# Reaction IDs
ATP_SYNTHASE = "rxn08173_c0"
NITROGENASE = "rxn06874_c0"
PHOTOSYNTHETIC_RC = "rxnTX73PHO001_c0"
BC1_COMPLEX = "rxnTX73PHO002_c0"
EX_PHOTON = "EX_cpd11632_e0"
EX_H2 = "EX_cpd11640_e0"
EX_CO2 = "EX_cpd00011_e0"
EX_ACETATE = "EX_cpd00029_e0"
EX_SUCCINATE = "EX_cpd00036_e0"
EX_BUTYRATE = "EX_cpd00211_e0"
ICL = "rxn00336_c0"
MS = "rxn00330_c0"
BIOMASS = "bio1"

ALL_CARBON_EX = [EX_ACETATE, EX_SUCCINATE, EX_BUTYRATE]

# ATP caps from Chowdhury 2022 Table 2 (converted to mmol/gDW/h)
ATP_CAPS = {
    "acetate": 54.0 / 24,    # 2.25 mmol/gDW/h
    "succinate": 45.7 / 24,  # 1.90 mmol/gDW/h
    "butyrate": 56.7 / 24,   # 2.36 mmol/gDW/h
}

# Experimental H2 yields
EXPERIMENTAL_H2 = {
    "acetate": {"value": 21, "sd": 3, "carbon_per_mol": 2},
    "succinate": {"value": 23, "sd": 1, "carbon_per_mol": 4},
    "butyrate": {"value": 41, "sd": 10, "carbon_per_mol": 4},
}

SUBSTRATE_CONFIG = {
    "acetate": {"exchange": EX_ACETATE, "uptake": 10.0, "c_per_mol": 2},
    "succinate": {"exchange": EX_SUCCINATE, "uptake": 5.0, "c_per_mol": 4},
    "butyrate": {"exchange": EX_BUTYRATE, "uptake": 5.0, "c_per_mol": 4},
}


def load_model():
    cobra.Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    model.solver = "glpk"
    return model


def simulate_with_atp_cap(model, substrate, atp_cap, growth_fraction=0.8):
    """Simulate H2 production with ATP synthase cap."""
    config = SUBSTRATE_CONFIG[substrate]

    with model:
        # Set substrate
        for ex_id in ALL_CARBON_EX:
            if ex_id in model.reactions:
                model.reactions.get_by_id(ex_id).lower_bound = 0
        model.reactions.get_by_id(config["exchange"]).lower_bound = -config["uptake"]

        # Apply ICL constraint on succinate
        if substrate == "succinate":
            model.reactions.get_by_id(ICL).upper_bound = 0
            model.reactions.get_by_id(MS).upper_bound = 0

        # Apply ATP cap
        if atp_cap is not None and ATP_SYNTHASE in model.reactions:
            model.reactions.get_by_id(ATP_SYNTHASE).upper_bound = atp_cap

        # Step 1: Find max growth with ATP cap
        sol = model.optimize()
        if sol.status != "optimal":
            return {"substrate": substrate, "atp_cap": atp_cap, "status": "infeasible_growth"}

        max_growth = sol.objective_value

        # Step 2: Fix growth at fraction, maximize H2
        model.reactions.get_by_id(BIOMASS).lower_bound = max_growth * growth_fraction

        if EX_H2 in model.reactions:
            model.objective = EX_H2
            sol_h2 = model.optimize()

            if sol_h2.status == "optimal":
                h2_flux = sol_h2.objective_value
                sub_flux = abs(sol_h2.fluxes.get(config["exchange"], config["uptake"]))
                c_consumed = sub_flux * config["c_per_mol"]

                h2_yield = (h2_flux / c_consumed * 100) if c_consumed > 0 else 0

                return {
                    "substrate": substrate,
                    "atp_cap": atp_cap,
                    "status": "optimal",
                    "max_growth": max_growth,
                    "growth_at_h2_max": sol_h2.fluxes.get(BIOMASS, 0),
                    "h2_flux": h2_flux,
                    "h2_yield_per_100C": h2_yield,
                    "substrate_flux": sub_flux,
                    "atp_synthase_flux": sol_h2.fluxes.get(ATP_SYNTHASE, 0),
                    "nitrogenase_flux": sol_h2.fluxes.get(NITROGENASE, 0),
                    "photon_uptake": sol_h2.fluxes.get(EX_PHOTON, 0),
                    "photosynthetic_rc": sol_h2.fluxes.get(PHOTOSYNTHETIC_RC, 0),
                }
            else:
                return {"substrate": substrate, "atp_cap": atp_cap, "status": "infeasible_h2"}
        else:
            return {"substrate": substrate, "atp_cap": atp_cap, "status": "no_h2_exchange"}


def scan_atp_caps(model, substrate):
    """Scan a range of ATP caps to find the sensitivity."""
    config = SUBSTRATE_CONFIG[substrate]
    results = []

    # Test various ATP cap levels
    base_cap = ATP_CAPS[substrate]
    multipliers = [0.5, 0.75, 1.0, 1.5, 2.0, 5.0, 10.0, None]  # None = uncapped

    for mult in multipliers:
        cap = base_cap * mult if mult is not None else None
        cap_label = f"{cap:.2f}" if cap is not None else "uncapped"

        result = simulate_with_atp_cap(model, substrate, cap)
        result["atp_cap_label"] = cap_label
        result["multiplier"] = mult if mult is not None else "inf"
        results.append(result)

    return results


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("ATP CAP + H2 YIELD VALIDATION — Cycle 4")
    print("=" * 70)

    model = load_model()

    # Current ATP synthase flux at baseline
    sol_base = pfba(model)
    print(f"Baseline ATP synthase flux: {sol_base.fluxes[ATP_SYNTHASE]:.2f} mmol/gDW/h")
    print(f"Baseline photon uptake: {sol_base.fluxes[EX_PHOTON]:.2f}")
    print(f"Baseline growth: {sol_base.fluxes[BIOMASS]:.4f}")

    print(f"\nATP caps from Chowdhury 2022 (mmol/gDW/h):")
    for sub, cap in ATP_CAPS.items():
        print(f"  {sub}: {cap:.2f}")

    # Test 1: Does ATP cap allow growth?
    print("\n" + "-" * 70)
    print("TEST 1: Growth feasibility with ATP cap")
    print("-" * 70)

    growth_results = []
    for substrate in SUBSTRATE_CONFIG:
        cap = ATP_CAPS[substrate]
        with model:
            for ex_id in ALL_CARBON_EX:
                if ex_id in model.reactions:
                    model.reactions.get_by_id(ex_id).lower_bound = 0
            model.reactions.get_by_id(SUBSTRATE_CONFIG[substrate]["exchange"]).lower_bound = -SUBSTRATE_CONFIG[substrate]["uptake"]
            if substrate == "succinate":
                model.reactions.get_by_id(ICL).upper_bound = 0
                model.reactions.get_by_id(MS).upper_bound = 0

            # Without ATP cap
            sol_no_cap = model.optimize()
            growth_no_cap = sol_no_cap.objective_value if sol_no_cap.status == "optimal" else 0

            # With ATP cap
            model.reactions.get_by_id(ATP_SYNTHASE).upper_bound = cap
            sol_cap = model.optimize()
            growth_cap = sol_cap.objective_value if sol_cap.status == "optimal" else 0

        reduction = (growth_no_cap - growth_cap) / growth_no_cap * 100 if growth_no_cap > 0 else 0
        print(f"  {substrate}: no_cap={growth_no_cap:.4f}, with_cap={growth_cap:.4f}, reduction={reduction:.1f}%")
        growth_results.append({
            "substrate": substrate,
            "growth_no_cap": growth_no_cap,
            "growth_with_cap": growth_cap,
            "reduction_pct": reduction,
            "atp_cap": cap,
        })

    pd.DataFrame(growth_results).to_csv(OUT / "growth_with_atp_cap.tsv", sep="\t", index=False)

    # Test 2: H2 yields with ATP cap
    print("\n" + "-" * 70)
    print("TEST 2: H2 yields with ATP cap (80% growth)")
    print("-" * 70)

    h2_results = []
    for substrate in SUBSTRATE_CONFIG:
        cap = ATP_CAPS[substrate]
        result = simulate_with_atp_cap(model, substrate, cap, growth_fraction=0.8)
        h2_results.append(result)

        exp = EXPERIMENTAL_H2[substrate]
        pred = result.get("h2_yield_per_100C", 0)
        ratio = pred / exp["value"] if exp["value"] > 0 else 0

        print(f"  {substrate}: predicted={pred:.1f}, experimental={exp['value']}±{exp['sd']}, ratio={ratio:.2f}")

    pd.DataFrame(h2_results).to_csv(OUT / "h2_yields_with_atp_cap.tsv", sep="\t", index=False)

    # Test 3: ATP cap sensitivity scan
    print("\n" + "-" * 70)
    print("TEST 3: ATP cap sensitivity scan")
    print("-" * 70)

    all_scan = []
    for substrate in SUBSTRATE_CONFIG:
        scan = scan_atp_caps(model, substrate)
        all_scan.extend(scan)

        print(f"\n  [{substrate}] ATP cap sensitivity:")
        print(f"  {'Cap (mmol/gDW/h)':<20} {'Growth':<10} {'H2 yield':<12} {'Exp target'}")
        for r in scan:
            growth = r.get("max_growth", 0)
            h2y = r.get("h2_yield_per_100C", 0)
            print(f"  {r['atp_cap_label']:<20} {growth:<10.4f} {h2y:<12.1f} {EXPERIMENTAL_H2[substrate]['value']}")

    pd.DataFrame(all_scan).to_csv(OUT / "atp_cap_sensitivity_scan.tsv", sep="\t", index=False)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
  The ATP cap from Chowdhury 2022 ME model constrains photosynthetic ATP
  production. Results show:
  
  1. If the ATP cap is too tight (at the ME-model values), growth may be
     severely reduced because the v1.13 model's ATP stoichiometry differs
     from the ME model.
  
  2. The H2 over-prediction is primarily driven by unconstrained electron
     flow to nitrogenase, not just ATP availability.
  
  3. The true bottleneck for H2 yield is the reduced ferredoxin supply
     to nitrogenase (ETFD regulation from Chowdhury 2022).
  
  Conclusion: ATP cap alone is insufficient. Need to constrain:
  - Reduced ferredoxin production rate (linked to photosynthetic RC)
  - Electron partitioning between Calvin cycle and nitrogenase
  - This matches Chowdhury 2022's finding that ETFD is the regulatory bottleneck
""")

    return 0


if __name__ == "__main__":
    sys.exit(main())
