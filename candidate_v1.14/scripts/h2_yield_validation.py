#!/usr/bin/env python3
"""
H2 Yield Quantitative Validation — v1.14 Cycle 3
==================================================
Compares model-predicted H2 yields against McKinlay/Harwood 2011 Table 1.

The NifA* strain produces H2 constitutively via nitrogenase.
To simulate this, we allow nitrogenase flux and maximize H2 production
while maintaining growth, then compare yields.

Key experimental conditions:
- R. palustris NifA* strain (CGA676)
- Anaerobic photoheterotrophic
- 7.5 mM (NH4)2SO4 as nitrogen source
- 40 mM carbon (e.g., 20 mM acetate, 10 mM succinate, 10 mM butyrate)
- 30°C, light, argon headspace

Experimental H2 yields (mol H2 / 100 mol organic C consumed):
- Fumarate: 18 ± 3
- Succinate: 23 ± 1
- Acetate: 21 ± 3
- Butyrate (no HCO3): 41 ± 10
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
OUT = ROOT / "candidate_v1.14" / "h2_validation"

# Reaction IDs
NITROGENASE = "rxn06874_c0"
EX_H2 = "EX_cpd11640_e0"
EX_CO2 = "EX_cpd00011_e0"
EX_ACETATE = "EX_cpd00029_e0"
EX_SUCCINATE = "EX_cpd00036_e0"
EX_BUTYRATE = "EX_cpd00211_e0"
BIOMASS = "bio1"

ALL_CARBON_EX = [EX_ACETATE, EX_SUCCINATE, EX_BUTYRATE]

# Experimental targets from McKinlay/Harwood 2011 Table 1
EXPERIMENTAL_H2_YIELDS = {
    "acetate": {"value": 21, "sd": 3, "carbon_per_mol": 2},
    "succinate": {"value": 23, "sd": 1, "carbon_per_mol": 4},
    "butyrate": {"value": 41, "sd": 10, "carbon_per_mol": 4},
}


def load_model():
    cobra.Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    model.solver = "glpk"
    return model


def simulate_h2_production(model, substrate_name, substrate_exchange, uptake_rate, carbon_per_mol):
    """
    Simulate H2 production by:
    1. First finding maximum growth
    2. Then fixing growth at 90% of max and maximizing H2
    3. Computing H2 yield per organic carbon consumed
    """
    results = {"substrate": substrate_name}

    with model:
        # Set substrate
        for ex_id in ALL_CARBON_EX:
            if ex_id in model.reactions:
                model.reactions.get_by_id(ex_id).lower_bound = 0
        model.reactions.get_by_id(substrate_exchange).lower_bound = -uptake_rate

        # Step 1: Find max growth (WT-like, no forced H2)
        sol_wt = model.optimize()
        if sol_wt.status != "optimal":
            results["status"] = "infeasible_WT"
            return results
        results["max_growth_WT"] = sol_wt.objective_value

        # Step 2: Check if nitrogenase can carry flux
        nit_rxn = model.reactions.get_by_id(NITROGENASE) if NITROGENASE in model.reactions else None
        if nit_rxn is None:
            results["status"] = "no_nitrogenase"
            return results

        # Step 3: Fix growth at fraction of max, maximize H2
        # NifA* grows slower than WT (Table 2 shows ~80-90% of WT rate)
        growth_fraction = 0.8  # NifA* grows ~80% of WT rate
        min_growth = sol_wt.objective_value * growth_fraction

        # Set biomass lower bound
        bio_rxn = model.reactions.get_by_id(BIOMASS)
        bio_rxn.lower_bound = min_growth

        # Change objective to maximize H2 export
        if EX_H2 in model.reactions:
            model.objective = EX_H2
            sol_h2 = model.optimize()

            if sol_h2.status == "optimal":
                h2_flux = sol_h2.objective_value  # mmol H2 / gDW / h
                substrate_flux = abs(sol_h2.fluxes.get(substrate_exchange, uptake_rate))

                # H2 yield = mol H2 / (mol substrate * C per mol) * 100
                if substrate_flux > 1e-10:
                    h2_yield = (h2_flux / (substrate_flux * carbon_per_mol)) * 100
                else:
                    h2_yield = 0

                results["status"] = "optimal"
                results["h2_flux"] = h2_flux
                results["substrate_flux"] = substrate_flux
                results["h2_yield_per_100C"] = h2_yield
                results["growth_at_h2_max"] = sol_h2.fluxes.get(BIOMASS, 0)
                results["co2_flux"] = sol_h2.fluxes.get(EX_CO2, 0)
                results["nitrogenase_flux"] = sol_h2.fluxes.get(NITROGENASE, 0)
            else:
                results["status"] = "infeasible_H2_max"
        else:
            results["status"] = "no_H2_exchange"

    return results


def compute_h2_fva(model, substrate_name, substrate_exchange, uptake_rate):
    """Compute FVA range for H2 production at various growth fractions."""
    fva_results = []

    with model:
        for ex_id in ALL_CARBON_EX:
            if ex_id in model.reactions:
                model.reactions.get_by_id(ex_id).lower_bound = 0
        model.reactions.get_by_id(substrate_exchange).lower_bound = -uptake_rate

        sol = model.optimize()
        if sol.status != "optimal":
            return fva_results

        max_growth = sol.objective_value

        for frac in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
            bio_rxn = model.reactions.get_by_id(BIOMASS)
            bio_rxn.lower_bound = max_growth * frac

            try:
                if EX_H2 in model.reactions:
                    fva = flux_variability_analysis(
                        model, reaction_list=[EX_H2], fraction_of_optimum=0.0
                    )
                    fva_results.append({
                        "substrate": substrate_name,
                        "growth_fraction": frac,
                        "min_growth": max_growth * frac,
                        "H2_min": fva["minimum"].iloc[0],
                        "H2_max": fva["maximum"].iloc[0],
                    })
            except Exception:
                pass

            bio_rxn.lower_bound = 0  # Reset for next iteration

    return fva_results


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("H2 YIELD QUANTITATIVE VALIDATION — Cycle 3")
    print("=" * 70)

    model = load_model()

    # Check H2 exchange exists and its bounds
    if EX_H2 in model.reactions:
        h2_rxn = model.reactions.get_by_id(EX_H2)
        print(f"H2 exchange: {h2_rxn.id}, bounds=[{h2_rxn.lower_bound}, {h2_rxn.upper_bound}]")
    else:
        print("WARNING: H2 exchange reaction not found!")
        # Search for it
        for rxn in model.exchanges:
            if "h2" in rxn.id.lower() or "hydrogen" in (rxn.name or "").lower():
                print(f"  Candidate: {rxn.id} ({rxn.name}) bounds=[{rxn.lower_bound}, {rxn.upper_bound}]")

    # Check nitrogenase
    if NITROGENASE in model.reactions:
        nit = model.reactions.get_by_id(NITROGENASE)
        print(f"Nitrogenase: {nit.id}, bounds=[{nit.lower_bound}, {nit.upper_bound}]")
        print(f"  Reaction: {nit.reaction}")
    else:
        print("WARNING: Nitrogenase not found!")

    print("\n" + "-" * 70)
    print("H2 YIELD SIMULATIONS")
    print("-" * 70)

    all_results = []
    all_fva = []

    substrates = {
        "acetate": (EX_ACETATE, 10.0, 2),
        "succinate": (EX_SUCCINATE, 5.0, 4),
        "butyrate": (EX_BUTYRATE, 5.0, 4),
    }

    for sub_name, (ex_id, uptake, c_per_mol) in substrates.items():
        print(f"\n  [{sub_name}]")
        result = simulate_h2_production(model, sub_name, ex_id, uptake, c_per_mol)
        all_results.append(result)

        exp = EXPERIMENTAL_H2_YIELDS[sub_name]
        pred = result.get("h2_yield_per_100C", 0)

        print(f"    Status: {result.get('status', 'unknown')}")
        print(f"    H2 flux: {result.get('h2_flux', 0):.4f} mmol/gDW/h")
        print(f"    Predicted H2 yield: {pred:.1f} mol H2 / 100 mol organic C")
        print(f"    Experimental:       {exp['value']} ± {exp['sd']}")

        if pred > 0:
            ratio = pred / exp["value"]
            within_sd = abs(pred - exp["value"]) <= 2 * exp["sd"]
            print(f"    Ratio (pred/exp): {ratio:.2f}")
            print(f"    Within 2 SD: {'YES' if within_sd else 'NO'}")

        # FVA
        fva_res = compute_h2_fva(model, sub_name, ex_id, uptake)
        all_fva.extend(fva_res)

    # Save results
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT / "h2_yield_predictions.tsv", sep="\t", index=False)

    if all_fva:
        fva_df = pd.DataFrame(all_fva)
        fva_df.to_csv(OUT / "h2_fva_growth_tradeoff.tsv", sep="\t", index=False)

    # Summary comparison table
    print("\n" + "-" * 70)
    print("SUMMARY: H2 Yield Comparison")
    print("-" * 70)
    print(f"\n  {'Substrate':<12} {'Predicted':<12} {'Experimental':<15} {'Within 2SD'}")
    print("  " + "-" * 55)
    for result in all_results:
        sub = result["substrate"]
        pred = result.get("h2_yield_per_100C", 0)
        exp = EXPERIMENTAL_H2_YIELDS[sub]
        within = "YES" if abs(pred - exp["value"]) <= 2 * exp["sd"] else "NO"
        print(f"  {sub:<12} {pred:<12.1f} {exp['value']} ± {exp['sd']:<10} {within}")

    print("\n" + "=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
