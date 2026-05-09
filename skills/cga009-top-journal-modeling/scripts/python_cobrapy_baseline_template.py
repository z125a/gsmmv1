from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis, pfba


ROOT = Path(__file__).resolve().parents[3]
MODEL = ROOT / "baseline" / "mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_for_memote.xml"
OUT = ROOT / "candidate_v1.14" / "python_cobrapy_baseline"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cobra.Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(str(MODEL))
    model.solver = "glpk"

    summary = {
        "model_id": model.id,
        "reaction_count": len(model.reactions),
        "metabolite_count": len(model.metabolites),
        "gene_count": len(model.genes),
        "solver": str(model.solver.interface),
        "objective": str(model.objective.expression),
    }
    pd.DataFrame([summary]).to_csv(OUT / "model_summary.tsv", sep="\t", index=False)

    solution = model.optimize()
    pd.DataFrame(
        [
            {
                "status": solution.status,
                "objective_value": solution.objective_value,
            }
        ]
    ).to_csv(OUT / "fba_solution_summary.tsv", sep="\t", index=False)

    try:
        pfba_solution = pfba(model)
        pfba_solution.fluxes.rename("flux").to_csv(OUT / "pfba_fluxes.tsv", sep="\t")
    except Exception as exc:
        pd.DataFrame([{"step": "pfba", "error": repr(exc)}]).to_csv(
            OUT / "pfba_error.tsv", sep="\t", index=False
        )

    try:
        objective_rxns = [rxn.id for rxn in model.reactions if rxn.objective_coefficient]
        reaction_list = objective_rxns or [rxn.id for rxn in model.reactions[:25]]
        fva = flux_variability_analysis(model, reaction_list=reaction_list, fraction_of_optimum=0.9)
        fva.to_csv(OUT / "fva_objective_or_first25.tsv", sep="\t")
    except Exception as exc:
        pd.DataFrame([{"step": "fva", "error": repr(exc)}]).to_csv(
            OUT / "fva_error.tsv", sep="\t", index=False
        )

    print(f"Wrote Python/COBRApy baseline outputs to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
