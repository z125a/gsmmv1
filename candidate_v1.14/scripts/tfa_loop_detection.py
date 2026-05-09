#!/usr/bin/env python3
"""
TFA Loop Detection and Thermodynamic Feasibility — v1.14 Cycle 6
==================================================================
Identifies thermodynamically infeasible internal loops in the v1.13 model
using FVA-based loop detection (loopless FBA approach).

Also applies known standard deltaG values for central carbon reactions
to check if any predicted flux directions violate thermodynamics.

Standard deltaG'° values (pH 7.0, 25°C, I=0.25M) from literature:
- Citrate synthase: -31.4 kJ/mol (strongly irreversible forward)
- Isocitrate dehydrogenase: -8.4 kJ/mol
- Isocitrate lyase: -1.1 kJ/mol (near equilibrium)
- Malate synthase: -32.5 kJ/mol (strongly irreversible forward)
- Pyruvate:Fd oxidoreductase: -19.2 kJ/mol (forward)
- RuBisCO: -35.1 kJ/mol (strongly irreversible forward)
- Succinate dehydrogenase: +0.4 kJ/mol (near equilibrium, reversible)
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
OUT = ROOT / "candidate_v1.14" / "tfa_results"

# Standard deltaG'° values (kJ/mol) at pH 7.0, 25°C, I=0.25M
# Sources: eQuilibrator, Alberty 2003, Flamholz et al. 2012
THERMO_DATA = {
    "rxn00256_c0": {"name": "Citrate synthase", "deltaG0": -31.4, "direction": "forward_only", "source": "eQuilibrator/Alberty"},
    "rxn00974_c0": {"name": "Aconitase", "deltaG0": 7.6, "direction": "reversible", "source": "eQuilibrator"},
    "rxn00199_c0": {"name": "Isocitrate dehydrogenase", "deltaG0": -8.4, "direction": "forward_preferred", "source": "eQuilibrator"},
    "rxn00336_c0": {"name": "Isocitrate lyase", "deltaG0": -1.1, "direction": "near_equilibrium", "source": "eQuilibrator"},
    "rxn00330_c0": {"name": "Malate synthase", "deltaG0": -32.5, "direction": "forward_only", "source": "eQuilibrator"},
    "rxn00288_c0": {"name": "Succinate dehydrogenase", "deltaG0": 0.4, "direction": "reversible", "source": "eQuilibrator"},
    "rxn00285_c0": {"name": "Succinyl-CoA synthetase", "deltaG0": -2.9, "direction": "near_equilibrium", "source": "eQuilibrator"},
    "rxn13974_c0": {"name": "Pyruvate:Fd oxidoreductase", "deltaG0": -19.2, "direction": "forward_preferred", "source": "Thauer 1977/eQuilibrator"},
    "rxn01116_c0": {"name": "RuBisCO (form I)", "deltaG0": -35.1, "direction": "forward_only", "source": "eQuilibrator"},
    "rxn08173_c0": {"name": "ATP synthase", "deltaG0": -36.0, "direction": "forward_only", "source": "Alberty"},
}


def load_model():
    cobra.Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    model.solver = "glpk"
    return model


def detect_type_III_loops(model):
    """
    Detect Type III (thermodynamically infeasible) loops using FVA.
    A Type III loop exists when a set of internal reactions can carry flux
    without any net exchange with the environment.
    
    Method: Set all exchange reactions to zero, then check which internal
    reactions can still carry flux (these form loops).
    """
    with model:
        # Block all exchanges
        for rxn in model.exchanges:
            rxn.bounds = (0, 0)

        # Also block demand/sink reactions
        for rxn in model.demands:
            rxn.bounds = (0, 0)

        # Set objective to zero (feasibility only)
        model.objective = model.reactions[0]
        model.objective.direction = "max"

        # Check which reactions can carry flux
        loop_reactions = []
        
        # Use FVA on all internal reactions
        internal_rxns = [r.id for r in model.reactions 
                        if not r.id.startswith("EX_") and not r.id.startswith("DM_")
                        and r.id != "bio1"]
        
        # FVA with fraction_of_optimum=0 (just feasibility)
        try:
            # Do FVA in batches to avoid memory issues
            batch_size = 100
            all_fva = []
            for i in range(0, len(internal_rxns), batch_size):
                batch = internal_rxns[i:i+batch_size]
                fva = flux_variability_analysis(
                    model, reaction_list=batch, fraction_of_optimum=0.0
                )
                all_fva.append(fva)
            
            fva_combined = pd.concat(all_fva)
            
            # Reactions that can carry flux with all exchanges closed = loop members
            for rxn_id in fva_combined.index:
                min_flux = fva_combined.loc[rxn_id, "minimum"]
                max_flux = fva_combined.loc[rxn_id, "maximum"]
                if abs(min_flux) > 1e-6 or abs(max_flux) > 1e-6:
                    loop_reactions.append({
                        "reaction_id": rxn_id,
                        "fva_min": min_flux,
                        "fva_max": max_flux,
                        "can_carry_flux": True,
                    })
            
            return loop_reactions, fva_combined
            
        except Exception as e:
            print(f"  FVA loop detection failed: {e}")
            return [], None


def check_thermo_consistency(model):
    """
    Check if pFBA flux directions are consistent with known deltaG values.
    """
    sol = pfba(model)
    
    results = []
    for rxn_id, thermo in THERMO_DATA.items():
        if rxn_id not in model.reactions:
            continue
        
        flux = sol.fluxes.get(rxn_id, 0)
        deltaG = thermo["deltaG0"]
        direction = thermo["direction"]
        
        # Check consistency
        # If deltaG < -10 kJ/mol: reaction should be forward (positive flux)
        # If deltaG > 10 kJ/mol: reaction should be reverse (negative flux)
        # If -10 < deltaG < 10: either direction is thermodynamically feasible
        
        if direction == "forward_only" and flux < -1e-6:
            status = "VIOLATION"
            note = f"Flux is reverse ({flux:.4f}) but deltaG strongly favors forward"
        elif direction == "forward_preferred" and flux < -1e-6:
            status = "WARNING"
            note = f"Flux is reverse ({flux:.4f}); possible but deltaG favors forward"
        elif direction == "reversible" or direction == "near_equilibrium":
            status = "OK"
            note = f"Flux={flux:.4f}; reaction is thermodynamically reversible"
        else:
            status = "OK"
            note = f"Flux={flux:.4f}; consistent with deltaG"
        
        results.append({
            "reaction_id": rxn_id,
            "name": thermo["name"],
            "deltaG0_kJ_mol": deltaG,
            "thermo_direction": direction,
            "pfba_flux": flux,
            "status": status,
            "note": note,
        })
    
    return results


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("TFA LOOP DETECTION & THERMODYNAMIC FEASIBILITY — Cycle 6")
    print("=" * 70)

    model = load_model()

    # Step 1: Check thermodynamic consistency of pFBA fluxes
    print("\n[1] Thermodynamic consistency check (pFBA vs deltaG)")
    print("-" * 50)

    thermo_results = check_thermo_consistency(model)
    thermo_df = pd.DataFrame(thermo_results)
    thermo_df.to_csv(OUT / "thermo_consistency_check.tsv", sep="\t", index=False)

    violations = [r for r in thermo_results if r["status"] == "VIOLATION"]
    warnings_list = [r for r in thermo_results if r["status"] == "WARNING"]

    for r in thermo_results:
        marker = "✗" if r["status"] == "VIOLATION" else ("!" if r["status"] == "WARNING" else "✓")
        print(f"  {marker} {r['name']:35s} | deltaG={r['deltaG0_kJ_mol']:+6.1f} | flux={r['pfba_flux']:+8.4f} | {r['status']}")

    print(f"\n  Violations: {len(violations)}, Warnings: {len(warnings_list)}")

    # Step 2: Detect internal loops
    print("\n[2] Type III loop detection (FVA with exchanges closed)")
    print("-" * 50)
    print("  Running FVA on internal reactions with all exchanges blocked...")

    loop_rxns, fva_data = detect_type_III_loops(model)

    if loop_rxns:
        print(f"  Found {len(loop_rxns)} reactions that can carry flux in loops")
        loop_df = pd.DataFrame(loop_rxns)
        loop_df.to_csv(OUT / "type_III_loop_reactions.tsv", sep="\t", index=False)

        # Check if any priority reactions are in loops
        priority_in_loops = [r for r in loop_rxns if r["reaction_id"] in THERMO_DATA]
        if priority_in_loops:
            print(f"\n  Priority reactions in loops:")
            for r in priority_in_loops:
                name = THERMO_DATA.get(r["reaction_id"], {}).get("name", "unknown")
                print(f"    {r['reaction_id']}: {name} (FVA: [{r['fva_min']:.2f}, {r['fva_max']:.2f}])")
    else:
        print("  No loops detected (or detection failed)")

    # Step 3: Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
  Thermodynamic consistency:
  - {len(violations)} violations (flux direction contradicts strong deltaG)
  - {len(warnings_list)} warnings (flux direction possible but unfavored)
  - {len(thermo_results) - len(violations) - len(warnings_list)} OK
  
  Internal loops: {len(loop_rxns)} reactions can carry flux without exchange
  
  TFA implementation status:
  - deltaG values assigned to {len(THERMO_DATA)} priority reactions
  - eQuilibrator API not available; using literature values
  - Full TFA (pyTFA) requires metabolite concentration data (not available)
  - Current approach: direction constraints only (no concentration-dependent bounds)
""")

    return 0


if __name__ == "__main__":
    sys.exit(main())
