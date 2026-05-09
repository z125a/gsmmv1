# Prompt For The Next Agent

You are taking over a Rhodopseudomonas palustris CGA009/TX73 genome-scale metabolic model project. Your job is to improve the model from a trustworthy frozen baseline to top-journal-level predictive credibility.

Start from this repository root. Use the bundled skill:

`skills/cga009-top-journal-modeling/SKILL.md`

Baseline:

`baseline/mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_20260508.xlsx`

SBML for MEMOTE/COBRA:

`baseline/mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_for_memote.xml`

Read first:

1. `baseline/CLAUDE_HANDOFF_V1.13_20260509.md`
2. `evidence/literature_evidence.tsv`
3. `evidence/LITERATURE_READING_PROTOCOL.md`
4. `evidence/phenotype_data_requirements.tsv`
5. `skills/cga009-top-journal-modeling/references/top_journal_strategy.md`

Mission:

Build the next publishable version by integrating three constraint layers: enzyme-capacity constraints (ecGEM/ME-model), thermodynamic flux analysis (TFA), and 13C-MFA validation. Do not chase a single reaction or single phenotype. Improve the whole model's reliability.

Execution environment:

Use Python/COBRApy as the primary path. Do not require MATLAB. If MATLAB/COBRA Toolbox is installed, it can be used as an additional cross-check, but lack of MATLAB is not a blocker. Implement reproducible Python scripts for SBML loading, FBA, pFBA, FVA, loopless checks where available, phenotype panels, evidence-table processing, and report generation. Use memote from Python/CLI when available. Prefer open solvers such as GLPK/HiGHS-compatible setups unless a licensed solver is already configured.

Persistence requirement:

Do not stop after one analysis pass, one report, or one candidate patch. Work in repeated build-test-criticize-rebuild cycles until the model reaches publication-ready quality or an explicit external blocker is reached. If a cycle fails, record why, revise the hypothesis, and continue with the next most defensible route. A normal stopping point is not "I found the next steps"; a normal stopping point is "the next version is reproducibly stronger than v1.13 and satisfies the acceptance gates below."

Autonomy requirement:

Do not repeatedly ask the user for yes/no approval during normal work. Proceed autonomously through reading literature, extracting evidence, writing scripts, running validation, creating candidate folders, committing reproducible artifacts, and iterating. Ask the user only for true blockers: missing proprietary data, unavailable credentials, paid database access, destructive overwrite, ambiguous biological convention that would permanently alter model interpretation, or a wet-lab decision that cannot be inferred from evidence.

Rules:

1. Keep v1.13 frozen until a candidate release passes stronger validation.
2. Every model edit must have a source, audit record, before/after validation, and counterargument.
3. Run adversarial analysis at every step: write the best argument against your change, then test it.
4. Separate measured constraints, literature-derived constraints, and sensitivity-only assumptions.
5. Never add exchanges/transports or change bounds just to force growth.
6. Do not make one-off edits in porphyrin/chlorophyll/cobalamin, tRNA pseudo-metabolites, or carrier-convention modules until pathway-wide conventions are solved.
7. Keep all scripts deterministic so another person can rerun from a clean clone.
8. Continue iterating across candidate versions (`candidate_v1.14`, `candidate_v1.15`, etc.) until the acceptance gates are met; do not hand back only a plan unless blocked by missing experimental data, missing credentials, or impossible software setup.
9. Do not stop because MATLAB is absent; use Python/COBRApy and document any MATLAB-only comparison as optional.
10. Do not pause for routine confirmations; keep working until a real stop condition is reached.

First deliverables:

1. A `candidate_v1.14/` folder with scripts, changed files, and audit TSVs.
2. A baseline reproduction report showing that v1.13 metrics can be regenerated.
3. A data-gap report mapping each desired ecGEM/TFA/13C-MFA constraint to available or missing data.
4. A no-hallucination evidence ledger with URLs/DOIs and exact fields extracted.
5. A rejection log of tempting but unsupported shortcuts.

Mandatory literature gate:

Before editing the model, create `candidate_v1.14/evidence_extraction.tsv` from the papers/databases listed in `evidence/literature_evidence.tsv` using the schema in `evidence/LITERATURE_READING_PROTOCOL.md`. The model can only be changed when the proposed edit or constraint has a corresponding evidence row. If literature is missing or inaccessible, record it and move to another evidence-supported task.

Success criterion:

The model is only better when it predicts growth phenotypes, exchange fluxes, and internal flux distributions more accurately without losing stoichiometric, thermodynamic, and biological credibility.

Publication-ready acceptance gates:

1. Baseline v1.13 can be reproduced from scripts.
2. Candidate version has no regression in the validated phenotype panel unless new experimental evidence justifies the change.
3. Exchange-flux predictions are compared against measured uptake/secretion data where available.
4. Internal flux predictions are tested against 13C-MFA or isotope-informed literature constraints.
5. ecGEM/ME constraints use documented proteomics and kcat confidence tiers, with sensitivity analysis for uncertain kcat.
6. TFA removes infeasible cycles without relying on unsupported metabolite concentration assumptions.
7. MEMOTE/COBRA reports, model diffs, evidence ledgers, and rejected-hypothesis logs are committed.
8. A manuscript-style validation narrative can defend the model against reviewer questions about overfitting, thermodynamic infeasibility, enzyme-capacity realism, and strain-specific evidence.

If these gates are not met, keep iterating.
