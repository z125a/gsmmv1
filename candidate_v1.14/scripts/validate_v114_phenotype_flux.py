#!/usr/bin/env python3
"""
Candidate v1.14 Phenotype & Internal Flux Validation Script
============================================================
Reproduces v1.13 baseline and builds the first internal-flux validation panel
using McKinlay/Harwood 2011 and Chowdhury 2022 data.

This script:
1. Loads the v1.13 SBML model
2. Reproduces baseline FBA/pFBA/FVA
3. Simulates photoheterotrophic growth on acetate, succinate, butyrate
4. Compares predicted flux distributions against 13C-MFA qualitative targets
5. Outputs a structured validation report

No model edits are made here — this is validation-only.
"""

from pathlib import Path
import sys
import warnings

import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import flux_variability_analysis, pfba

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
MODEL_PATH = ROOT / "baseline" / "mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_for_memote.xml"
OUT = ROOT / "candidate_v1.14" / "validation_reports"


def load_model():
    """Load v1.13 SBML and configure solver."""
    cobra.Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    model.solver = "glpk"
    return model


def find_exchange_reactions(model):
    """Identify exchange reactions and map metabolite IDs."""
    exchanges = {}
    for rxn in model.exchanges:
        for met in rxn.metabolites:
            met_name = met.name.lower() if met.name else met.id.lower()
            exchanges[rxn.id] = {
                "reaction_id": rxn.id,
                "metabolite_id": met.id,
                "metabolite_name": met.name,
                "lower_bound": rxn.lower_bound,
                "upper_bound": rxn.upper_bound,
            }
    return exchanges


def find_reaction_by_substring(model, substrings):
    """Find reactions whose ID or name contains any of the given substrings."""
    results = []
    for rxn in model.reactions:
        rxn_text = (rxn.id + " " + (rxn.name or "")).lower()
        for sub in substrings:
            if sub.lower() in rxn_text:
                results.append(rxn)
                break
    return results


def baseline_reproduction(model):
    """Reproduce baseline FBA metrics."""
    sol = model.optimize()
    results = {
        "status": sol.status,
        "objective_value": sol.objective_value,
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "genes": len(model.genes),
    }
    return results, sol


def run_pfba_and_fva(model):
    """Run pFBA and FVA on the objective."""
    pfba_sol = pfba(model)

    obj_rxns = [r for r in model.reactions if r.objective_coefficient != 0]
    if not obj_rxns:
        obj_rxns = [model.reactions[0]]

    fva_result = flux_variability_analysis(
        model,
        reaction_list=[r.id for r in obj_rxns],
        fraction_of_optimum=0.9,
    )
    return pfba_sol, fva_result


def identify_key_reactions(model):
    """Identify reactions relevant to 13C-MFA validation axes."""
    key_pathways = {
        "Calvin_cycle_RuBisCO": ["rubisco", "cbb", "ribulose"],
        "Calvin_cycle_GAPDH": ["gapdh", "glyceraldehyde-3-phosphate dehydrogenase"],
        "TCA_citrate_synthase": ["citrate synthase", "cs_"],
        "TCA_isocitrate_dehydrogenase": ["isocitrate dehydrogenase"],
        "TCA_alpha_ketoglutarate_dehydrogenase": ["alpha-ketoglutarate", "2-oxoglutarate dehydrogenase"],
        "TCA_succinate_dehydrogenase": ["succinate dehydrogenase", "sdh"],
        "glyoxylate_shunt_ICL": ["isocitrate lyase", "icl"],
        "glyoxylate_shunt_MS": ["malate synthase"],
        "pyruvate_dehydrogenase": ["pyruvate dehydrogenase", "pdh"],
        "pyruvate_ferredoxin_oxidoreductase": ["pyruvate:ferredoxin", "por"],
        "nitrogenase": ["nitrogenase", "nif"],
        "H2_exchange": ["h2", "hydrogen"],
        "CO2_exchange": ["co2", "carbon dioxide"],
        "acetate_exchange": ["acetate", "acet"],
        "succinate_exchange": ["succinate", "succ"],
        "butyrate_exchange": ["butyrate", "butyr", "butanoate"],
        "PHB_synthesis": ["phb", "polyhydroxybutyrate", "poly-beta-hydroxy"],
        "phosphoribulokinase": ["phosphoribulokinase", "prk"],
    }

    found = {}
    for pathway_name, search_terms in key_pathways.items():
        matches = find_reaction_by_substring(model, search_terms)
        found[pathway_name] = [(r.id, r.name, r.lower_bound, r.upper_bound) for r in matches]

    return found


def simulate_substrate_condition(model, substrate_exchange_id, substrate_uptake_rate,
                                  other_exchanges_to_close=None):
    """
    Simulate growth on a specific substrate by setting exchange bounds.
    Returns FBA solution or None if infeasible.
    """
    with model:
        # Close other carbon sources
        if other_exchanges_to_close:
            for ex_id in other_exchanges_to_close:
                if ex_id in model.reactions:
                    model.reactions.get_by_id(ex_id).lower_bound = 0

        # Set substrate uptake
        if substrate_exchange_id in model.reactions:
            model.reactions.get_by_id(substrate_exchange_id).lower_bound = -abs(substrate_uptake_rate)

        sol = model.optimize()
        if sol.status == "optimal":
            # Also get pFBA
            try:
                pfba_sol = pfba(model)
                return pfba_sol
            except Exception:
                return sol
        return None


def analyze_flux_distribution(solution, key_reactions):
    """Extract flux values for key pathway reactions from a solution."""
    flux_summary = {}
    if solution is None:
        return flux_summary

    fluxes = solution.fluxes if hasattr(solution, 'fluxes') else solution.fluxes

    for pathway, rxn_list in key_reactions.items():
        pathway_flux = 0.0
        rxn_ids_found = []
        for rxn_id, rxn_name, lb, ub in rxn_list:
            if rxn_id in fluxes.index:
                f = fluxes[rxn_id]
                pathway_flux += f
                rxn_ids_found.append((rxn_id, f))
        flux_summary[pathway] = {
            "total_flux": pathway_flux,
            "reactions": rxn_ids_found,
        }

    return flux_summary


def qualitative_mfa_comparison(flux_summary, substrate):
    """
    Compare predicted flux patterns against McKinlay/Harwood qualitative targets.
    Returns list of (test_name, expected, observed, pass/fail/inconclusive).
    """
    results = []

    # Test 1: Glyoxylate shunt should be active on acetate/butyrate, inactive on succinate
    gs_flux = flux_summary.get("glyoxylate_shunt_ICL", {}).get("total_flux", 0)
    gs_ms_flux = flux_summary.get("glyoxylate_shunt_MS", {}).get("total_flux", 0)
    gs_total = abs(gs_flux) + abs(gs_ms_flux)

    if substrate in ["acetate", "butyrate"]:
        expected = "active (>0)"
        status = "PASS" if gs_total > 1e-6 else "FAIL"
    elif substrate == "succinate":
        expected = "inactive (=0)"
        status = "PASS" if gs_total < 1e-6 else "FAIL"
    else:
        expected = "unknown"
        status = "INCONCLUSIVE"

    results.append({
        "test": f"glyoxylate_shunt_on_{substrate}",
        "expected": expected,
        "observed": f"{gs_total:.4f}",
        "status": status,
    })

    # Test 2: Calvin cycle (RuBisCO) should carry flux in WT photoheterotrophic growth
    cbb_flux = flux_summary.get("Calvin_cycle_RuBisCO", {}).get("total_flux", 0)
    expected_cbb = "active (>0) for WT photoheterotrophic"
    status_cbb = "PASS" if abs(cbb_flux) > 1e-6 else "FAIL"
    results.append({
        "test": f"Calvin_cycle_active_on_{substrate}",
        "expected": expected_cbb,
        "observed": f"{cbb_flux:.4f}",
        "status": status_cbb,
    })

    # Test 3: H2 exchange should be possible (not blocked)
    h2_flux = flux_summary.get("H2_exchange", {}).get("total_flux", 0)
    results.append({
        "test": f"H2_exchange_feasible_on_{substrate}",
        "expected": "non-zero if nitrogenase active",
        "observed": f"{h2_flux:.4f}",
        "status": "INFO",
    })

    # Test 4: CO2 should be produced (substrate oxidation)
    co2_flux = flux_summary.get("CO2_exchange", {}).get("total_flux", 0)
    results.append({
        "test": f"CO2_production_on_{substrate}",
        "expected": "positive (net CO2 release for WT on most substrates)",
        "observed": f"{co2_flux:.4f}",
        "status": "PASS" if co2_flux > 0 else "INFO",
    })

    return results


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("CANDIDATE v1.14 — PHENOTYPE & INTERNAL FLUX VALIDATION")
    print("=" * 70)

    # Step 1: Load model
    print("\n[1] Loading v1.13 SBML model...")
    model = load_model()
    print(f"    Loaded: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, {len(model.genes)} genes")

    # Step 2: Baseline reproduction
    print("\n[2] Baseline FBA reproduction...")
    baseline, sol = baseline_reproduction(model)
    print(f"    Status: {baseline['status']}")
    print(f"    Objective (bio1): {baseline['objective_value']:.6f}")

    # Step 3: pFBA and FVA
    print("\n[3] Running pFBA and FVA...")
    pfba_sol, fva_result = run_pfba_and_fva(model)
    print(f"    pFBA total flux: {pfba_sol.fluxes.abs().sum():.2f}")
    print(f"    FVA on objective: min={fva_result['minimum'].iloc[0]:.4f}, max={fva_result['maximum'].iloc[0]:.4f}")

    # Step 4: Identify key reactions
    print("\n[4] Identifying key pathway reactions for MFA comparison...")
    key_rxns = identify_key_reactions(model)
    key_rxn_report = []
    for pathway, rxns in key_rxns.items():
        key_rxn_report.append({
            "pathway": pathway,
            "n_reactions_found": len(rxns),
            "reaction_ids": "; ".join([r[0] for r in rxns]) if rxns else "NONE_FOUND",
        })
    key_rxn_df = pd.DataFrame(key_rxn_report)
    key_rxn_df.to_csv(OUT / "key_pathway_reactions.tsv", sep="\t", index=False)
    print(f"    Found reactions for {sum(1 for r in key_rxn_report if r['n_reactions_found'] > 0)}/{len(key_rxn_report)} pathways")

    # Step 5: Identify exchange reactions
    print("\n[5] Cataloging exchange reactions...")
    exchanges = find_exchange_reactions(model)
    ex_df = pd.DataFrame(exchanges.values())
    ex_df.to_csv(OUT / "exchange_reactions_catalog.tsv", sep="\t", index=False)
    print(f"    Found {len(exchanges)} exchange reactions")

    # Step 6: Substrate-specific simulations
    print("\n[6] Substrate-specific growth simulations...")

    # Find exchange reaction IDs for key substrates
    substrate_exchanges = {}
    for ex_id, ex_info in exchanges.items():
        name = (ex_info["metabolite_name"] or "").lower()
        mid = ex_info["metabolite_id"].lower()
        if "acetat" in name or "acetat" in mid:
            substrate_exchanges.setdefault("acetate", []).append(ex_id)
        if "succinat" in name or "succinat" in mid:
            substrate_exchanges.setdefault("succinate", []).append(ex_id)
        if "butyrat" in name or "butanoat" in name or "butyrat" in mid or "butanoat" in mid:
            substrate_exchanges.setdefault("butyrate", []).append(ex_id)

    print(f"    Substrate exchange candidates: {dict((k, v) for k, v in substrate_exchanges.items())}")

    all_mfa_results = []
    all_flux_summaries = []

    for substrate_name, ex_ids in substrate_exchanges.items():
        if not ex_ids:
            continue
        # Use the first matching exchange
        ex_id = ex_ids[0]
        print(f"\n    Simulating growth on {substrate_name} (exchange: {ex_id})...")

        # Close other carbon exchanges
        other_carbon = []
        for other_sub, other_ids in substrate_exchanges.items():
            if other_sub != substrate_name:
                other_carbon.extend(other_ids)

        sol_sub = simulate_substrate_condition(
            model, ex_id, substrate_uptake_rate=10.0,
            other_exchanges_to_close=other_carbon
        )

        if sol_sub is not None:
            obj_val = sol_sub.fluxes.get("bio1", sol_sub.objective_value if hasattr(sol_sub, 'objective_value') else 0)
            print(f"      Growth: {obj_val:.4f}")

            flux_summary = analyze_flux_distribution(sol_sub, key_rxns)
            mfa_tests = qualitative_mfa_comparison(flux_summary, substrate_name)
            all_mfa_results.extend(mfa_tests)

            for pathway, info in flux_summary.items():
                all_flux_summaries.append({
                    "substrate": substrate_name,
                    "pathway": pathway,
                    "total_flux": info["total_flux"],
                    "n_reactions": len(info["reactions"]),
                })
        else:
            print(f"      INFEASIBLE — no growth on {substrate_name}")
            all_mfa_results.append({
                "test": f"growth_on_{substrate_name}",
                "expected": "feasible",
                "observed": "infeasible",
                "status": "FAIL",
            })

    # Step 7: Write reports
    print("\n[7] Writing validation reports...")

    # Baseline report
    pd.DataFrame([baseline]).to_csv(OUT / "baseline_metrics.tsv", sep="\t", index=False)

    # MFA comparison
    mfa_df = pd.DataFrame(all_mfa_results)
    mfa_df.to_csv(OUT / "mfa_qualitative_comparison.tsv", sep="\t", index=False)

    # Flux summaries
    flux_df = pd.DataFrame(all_flux_summaries)
    flux_df.to_csv(OUT / "substrate_flux_summaries.tsv", sep="\t", index=False)

    # Summary statistics
    if len(mfa_df) > 0:
        pass_count = (mfa_df["status"] == "PASS").sum()
        fail_count = (mfa_df["status"] == "FAIL").sum()
        info_count = (mfa_df["status"] == "INFO").sum()
        total = len(mfa_df)
        print(f"\n    MFA Qualitative Tests: {pass_count} PASS / {fail_count} FAIL / {info_count} INFO / {total} total")
    else:
        print("\n    No MFA tests could be run.")

    print(f"\n    Reports written to: {OUT}")
    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE — No model edits made.")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
