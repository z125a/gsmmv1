#!/usr/bin/env python3
"""
Ferredoxin Electron Transport Constraint — v1.14 Cycle 5
==========================================================
Implements the ETFD (electron transport through ferredoxin) constraint
identified by Chowdhury 2022 as the regulatory bottleneck for H2 production.

Strategy:
1. In WT photoheterotrophic growth, the model shows Fd_red production
   balanced by Fd:NADP+ oxidoreductase and 2-OG:Fd oxidoreductase.
2. The total Fd_red production rate at baseline (WT, no H2) represents
   the maximum electron flow through ferredoxin under normal conditions.
3. For NifA* simulation, we cap total Fd_red production at the WT level
   and allow nitrogenase to compete for those electrons.

This implements the biological reality that:
- Photosynthetic electron transport has a fixed capacity
- Nitrogenase competes with Calvin cycle for electrons (McKinlay/Harwood)
- ETFD is the bottleneck (Chowdhury 2022)

Evidence:
- McKinlay/Harwood 2011 Fig 2: electron balance within 10% for all substrates
- Chowdhury 2022: ETFD as regulatory bottleneck
- McKinlay/Harwood 2011: Calvin cycle decrease accounts for H2 electrons
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
OUT = ROOT / "candidate_v1.14" / "etfd_results"

# Reaction IDs
POR = "rxn13974_c0"           # pyruvate:Fd oxidoreductase (produces 2 Fd_red in forward)
OG_FD = "rxn14048_c0"         # 2-oxoglutarate:Fd oxidoreductase (produces 2 Fd_red in reverse)
FD_NADP = "rxn14159_c0"       # Fd:NADP+ oxidoreductase (produces 2 Fd_red in reverse)
NITROGENASE = "rxn06874_c0"   # consumes 8 Fd_red
HYDROGENASE = "rxn05759_c0"   # consumes 2 Fd_red
HMB_PP = "rxn14070_c0"        # consumes 2 Fd_red
CO_FD = "rxn07189_c0"         # CO:Fd oxidoreductase
NH3_FD = "rxn05893_c0"        # NH3:Fd oxidoreductase
CRESOL_FD = "rxn14227_c0"     # 4-cresol:Fd oxidoreductase
CHLORO_FD = "rxn15930_c0"     # chlorophyll:Fd oxidoreductase

ICL = "rxn00336_c0"
MS = "rxn00330_c0"
EX_H2 = "EX_cpd11640_e0"
EX_ACETATE = "EX_cpd00029_e0"
EX_SUCCINATE = "EX_cpd00036_e0"
EX_BUTYRATE = "EX_cpd00211_e0"
BIOMASS = "bio1"

ALL_CARBON_EX = [EX_ACETATE, EX_SUCCINATE, EX_BUTYRATE]

SUBSTRATE_CONFIG = {
    "acetate": {"exchange": EX_ACETATE, "uptake": 10.0, "c_per_mol": 2},
    "succinate": {"exchange": EX_SUCCINATE, "uptake": 5.0, "c_per_mol": 4},
    "butyrate": {"exchange": EX_BUTYRATE, "uptake": 5.0, "c_per_mol": 4},
}

EXPERIMENTAL_H2 = {
    "acetate": {"value": 21, "sd": 3},
    "succinate": {"value": 23, "sd": 1},
    "butyrate": {"value": 41, "sd": 10},
}


def load_model():
    cobra.Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    model.solver = "glpk"
    return model


def get_wt_fd_red_production(model, substrate):
    """
    Get the total Fd_red production rate in WT (no H2) growth.
    This represents the maximum electron flow through ferredoxin.
    """
    config = SUBSTRATE_CONFIG[substrate]

    with model:
        for ex_id in ALL_CARBON_EX:
            if ex_id in model.reactions:
                model.reactions.get_by_id(ex_id).lower_bound = 0
        model.reactions.get_by_id(config["exchange"]).lower_bound = -config["uptake"]

        if substrate == "succinate":
            model.reactions.get_by_id(ICL).upper_bound = 0
            model.reactions.get_by_id(MS).upper_bound = 0

        # Block nitrogenase and hydrogenase (WT conditions)
        model.reactions.get_by_id(NITROGENASE).upper_bound = 0
        model.reactions.get_by_id(HYDROGENASE).upper_bound = 0

        sol = pfba(model)

        # Calculate total Fd_red production
        # Fd_red is produced when:
        # - POR runs forward (pyruvate → acetyl-CoA): +2 Fd_red per flux unit
        # - OG_FD runs reverse (negative flux): +2 Fd_red per |flux| unit
        # - FD_NADP runs reverse (negative flux): +2 Fd_red per |flux| unit
        # - CO_FD runs forward: +1 Fd_red
        # - NH3_FD runs forward: +6 Fd_red
        # - CRESOL_FD runs forward: +4 Fd_red

        fd_production = 0
        por_flux = sol.fluxes.get(POR, 0)
        if por_flux > 0:  # Forward = produces Fd_red
            fd_production += por_flux * 2  # But wait - POR consumes Fd_red in forward!

        # Actually let me re-check the stoichiometry
        # rxn13974: pyruvate + CoA + 2 Fd_ox → acetyl-CoA + CO2 + 2 Fd_red
        # So FORWARD produces Fd_red (positive flux = Fd_red production)
        # But the model has it as: consumes 2 Fd_red in the metabolite dict
        # Let me check the actual stoichiometry sign

        # From the earlier analysis:
        # rxn13974_c0: pyruvate:ferredoxin 2-oxidoreductase (CoA-acetylating) (consumes 2.0 Fd_red)
        # This means in the model's convention, Fd_red is a REACTANT (consumed)
        # So when POR flux is positive (forward), it CONSUMES Fd_red
        # When POR flux is negative (reverse), it PRODUCES Fd_red

        # Let me recalculate based on actual stoichiometry:
        # Fd_red CONSUMERS (positive flux = consuming Fd_red):
        #   POR (rxn13974): 2 Fd_red consumed per unit forward flux
        #   Nitrogenase (rxn06874): 8 Fd_red consumed per unit
        #   Hydrogenase (rxn05759): 2 Fd_red consumed per unit
        #   HMB-PP (rxn14070): 2 Fd_red consumed per unit
        #   Acetylene reductase (rxn06926): 1 Fd_red consumed per unit

        # Fd_red PRODUCERS (negative flux = producing Fd_red for these):
        #   OG_FD (rxn14048): 2 Fd_red produced per unit forward flux
        #     Wait - let me re-read: "produces 2 Fd_red (reverse)"
        #     So rxn14048 forward CONSUMES Fd_red, reverse PRODUCES it
        #     Baseline flux = -55.27, so it's running reverse = PRODUCING Fd_red

        #   FD_NADP (rxn14159): "produces 2 Fd_red (reverse)"
        #     Baseline flux = -34.28, running reverse = PRODUCING Fd_red

        # Actually from the metabolite analysis:
        # cpd11620[c0] (Fd_red):
        #   rxn13974: consumes 2 Fd_red (positive flux = consuming)
        #   rxn14048: consumes 2 Fd_red (positive flux = consuming)
        #   rxn14159: consumes 2 Fd_red (positive flux = consuming)
        #   rxn06874: consumes 8 Fd_red
        #   rxn05759: consumes 2 Fd_red
        #   rxn14070: consumes 2 Fd_red

        # cpd11620[c0] (Fd_red):
        #   rxn07189: produces 1 Fd_red
        #   rxn05893: produces 6 Fd_red
        #   rxn14227: produces 4 Fd_red

        # So the Fd_red PRODUCERS are:
        #   rxn07189 (CO:Fd): +1 per unit forward flux
        #   rxn05893 (NH3:Fd): +6 per unit forward flux
        #   rxn14227 (cresol:Fd): +4 per unit forward flux
        #   rxn13974 (POR): +2 per unit REVERSE flux (negative flux)
        #   rxn14048 (OG:Fd): +2 per unit REVERSE flux (negative flux)
        #   rxn14159 (Fd:NADP): +2 per unit REVERSE flux (negative flux)

        # And Fd_red CONSUMERS are:
        #   rxn13974 (POR): 2 per unit FORWARD flux (positive flux)
        #   rxn14048 (OG:Fd): 2 per unit FORWARD flux
        #   rxn14159 (Fd:NADP): 2 per unit FORWARD flux
        #   rxn06874 (nitrogenase): 8 per unit
        #   rxn05759 (hydrogenase): 2 per unit
        #   rxn14070 (HMB-PP): 2 per unit
        #   rxn06926 (acetylene): 1 per unit

        # At baseline (WT, pFBA):
        # POR: +89.05 (forward, CONSUMING 2*89.05 = 178.1 Fd_red)
        # OG_FD: -55.27 (reverse, PRODUCING 2*55.27 = 110.5 Fd_red)
        # FD_NADP: -34.28 (reverse, PRODUCING 2*34.28 = 68.6 Fd_red)
        # HMB-PP: +0.51 (forward, CONSUMING 2*0.51 = 1.0 Fd_red)
        # CO_FD: +0.007 (forward, PRODUCING 1*0.007 = 0.007 Fd_red)

        # Net Fd_red balance: production - consumption
        # Production: 110.5 + 68.6 + 0.007 = 179.1
        # Consumption: 178.1 + 1.0 = 179.1
        # Balanced! (as expected for steady state)

        # The KEY insight: total Fd_red PRODUCTION = total Fd_red CONSUMPTION
        # In WT, all Fd_red goes to POR (reverse direction would produce it,
        # but actually POR is forward and CONSUMING it)

        # Wait, I'm confusing myself. Let me just compute net Fd_red production
        # from the actual fluxes:

        fluxes = sol.fluxes

        # Net Fd_red production (positive = producing Fd_red)
        fd_red_net = {}
        # Reactions where Fd_red is a PRODUCT (positive stoichiometry)
        fd_red_net["CO_FD"] = fluxes.get(CO_FD, 0) * 1      # produces 1 when forward
        fd_red_net["NH3_FD"] = fluxes.get(NH3_FD, 0) * 6    # produces 6 when forward
        fd_red_net["CRESOL_FD"] = fluxes.get(CRESOL_FD, 0) * 4  # produces 4 when forward

        # Reactions where Fd_red is a REACTANT (negative stoichiometry)
        # When these run REVERSE (negative flux), they produce Fd_red
        fd_red_net["POR"] = -fluxes.get(POR, 0) * 2         # produces 2 when reverse
        fd_red_net["OG_FD"] = -fluxes.get(OG_FD, 0) * 2     # produces 2 when reverse
        fd_red_net["FD_NADP"] = -fluxes.get(FD_NADP, 0) * 2 # produces 2 when reverse

        # Consumption
        fd_red_net["NITROGENASE"] = -fluxes.get(NITROGENASE, 0) * 8  # consumes 8
        fd_red_net["HYDROGENASE"] = -fluxes.get(HYDROGENASE, 0) * 2  # consumes 2
        fd_red_net["HMB_PP"] = -fluxes.get(HMB_PP, 0) * 2           # consumes 2

        total_production = sum(v for v in fd_red_net.values() if v > 0)
        total_consumption = sum(-v for v in fd_red_net.values() if v < 0)

        return {
            "substrate": substrate,
            "growth": sol.fluxes.get(BIOMASS, 0),
            "fd_red_production": total_production,
            "fd_red_consumption": total_consumption,
            "fd_red_details": fd_red_net,
            "por_flux": fluxes.get(POR, 0),
            "og_fd_flux": fluxes.get(OG_FD, 0),
            "fd_nadp_flux": fluxes.get(FD_NADP, 0),
        }


def simulate_h2_with_fd_constraint(model, substrate, fd_cap_fraction=1.0):
    """
    Simulate H2 production with a ferredoxin electron transport constraint.

    The constraint limits total Fd_red availability by capping the reverse
    flux of Fd:NADP+ oxidoreductase (the main Fd_red source from NADPH).
    """
    config = SUBSTRATE_CONFIG[substrate]

    # First get WT Fd_red production
    wt_info = get_wt_fd_red_production(model, substrate)
    wt_fd_nadp_flux = wt_info["fd_nadp_flux"]  # This is negative (reverse)

    with model:
        for ex_id in ALL_CARBON_EX:
            if ex_id in model.reactions:
                model.reactions.get_by_id(ex_id).lower_bound = 0
        model.reactions.get_by_id(config["exchange"]).lower_bound = -config["uptake"]

        if substrate == "succinate":
            model.reactions.get_by_id(ICL).upper_bound = 0
            model.reactions.get_by_id(MS).upper_bound = 0

        # Apply Fd constraint: limit Fd:NADP+ reverse flux to WT level * fraction
        # wt_fd_nadp_flux is negative (reverse direction)
        # We want to allow it to go MORE negative (more Fd_red production) only up to WT level
        fd_nadp_rxn = model.reactions.get_by_id(FD_NADP)
        fd_nadp_rxn.lower_bound = wt_fd_nadp_flux * fd_cap_fraction

        # Find max growth with constraint
        sol = model.optimize()
        if sol.status != "optimal":
            return {"substrate": substrate, "fd_cap_fraction": fd_cap_fraction, "status": "infeasible_growth"}

        max_growth = sol.objective_value

        # Fix growth at 80% and maximize H2
        model.reactions.get_by_id(BIOMASS).lower_bound = max_growth * 0.8

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
                    "fd_cap_fraction": fd_cap_fraction,
                    "status": "optimal",
                    "max_growth": max_growth,
                    "h2_flux": h2_flux,
                    "h2_yield_per_100C": h2_yield,
                    "nitrogenase_flux": sol_h2.fluxes.get(NITROGENASE, 0),
                    "fd_nadp_flux": sol_h2.fluxes.get(FD_NADP, 0),
                    "por_flux": sol_h2.fluxes.get(POR, 0),
                }

    return {"substrate": substrate, "fd_cap_fraction": fd_cap_fraction, "status": "failed"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("FERREDOXIN ELECTRON TRANSPORT CONSTRAINT — Cycle 5")
    print("=" * 70)

    model = load_model()

    # Step 1: Characterize WT Fd_red production
    print("\n[1] WT Ferredoxin balance (no H2 production)")
    print("-" * 50)

    wt_results = []
    for substrate in SUBSTRATE_CONFIG:
        info = get_wt_fd_red_production(model, substrate)
        wt_results.append(info)
        print(f"\n  [{substrate}]")
        print(f"    Growth: {info['growth']:.4f}")
        print(f"    Fd_red production: {info['fd_red_production']:.2f}")
        print(f"    Fd_red consumption: {info['fd_red_consumption']:.2f}")
        print(f"    POR flux: {info['por_flux']:.2f}")
        print(f"    OG:Fd flux: {info['og_fd_flux']:.2f} (negative=producing Fd_red)")
        print(f"    Fd:NADP flux: {info['fd_nadp_flux']:.2f} (negative=producing Fd_red)")

    # Step 2: H2 with Fd constraint at various levels
    print("\n\n[2] H2 yields with Fd:NADP+ constraint")
    print("-" * 50)

    all_h2_results = []
    for substrate in SUBSTRATE_CONFIG:
        print(f"\n  [{substrate}]")
        print(f"  {'Fd cap fraction':<18} {'H2 yield':<12} {'Experimental':<15} {'Ratio'}")

        for frac in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 5.0]:
            result = simulate_h2_with_fd_constraint(model, substrate, frac)
            all_h2_results.append(result)

            h2y = result.get("h2_yield_per_100C", 0)
            exp = EXPERIMENTAL_H2[substrate]["value"]
            ratio = h2y / exp if exp > 0 else 0
            print(f"  {frac:<18.2f} {h2y:<12.1f} {exp:<15} {ratio:.2f}")

    pd.DataFrame(all_h2_results).to_csv(OUT / "h2_with_fd_constraint.tsv", sep="\t", index=False)

    # Step 3: Find the Fd constraint level that matches experimental H2
    print("\n\n[3] Finding Fd constraint that matches experimental H2 yields")
    print("-" * 50)

    for substrate in SUBSTRATE_CONFIG:
        exp_h2 = EXPERIMENTAL_H2[substrate]["value"]
        best_frac = None
        best_diff = float("inf")

        for frac in np.arange(0.01, 3.0, 0.05):
            result = simulate_h2_with_fd_constraint(model, substrate, frac)
            h2y = result.get("h2_yield_per_100C", 0)
            diff = abs(h2y - exp_h2)
            if diff < best_diff:
                best_diff = diff
                best_frac = frac

        # Run at best fraction
        best_result = simulate_h2_with_fd_constraint(model, substrate, best_frac)
        h2y = best_result.get("h2_yield_per_100C", 0)
        print(f"  {substrate}: best_frac={best_frac:.2f}, H2_yield={h2y:.1f} (target={exp_h2})")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
