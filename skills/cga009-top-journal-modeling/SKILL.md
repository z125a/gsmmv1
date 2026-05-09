---
name: cga009-top-journal-modeling
description: Use when advancing the Rhodopseudomonas palustris CGA009/TX73 genome-scale metabolic model from the frozen v1.13 baseline toward publication-grade predictive accuracy using adversarial model curation, ecGEM/ME-model constraints, thermodynamic flux analysis, 13C-MFA validation, MEMOTE/COBRA checks, and evidence-tracked GitHub handoff artifacts.
---

# CGA009 Top-Journal Modeling

## Operating Rule

Treat `mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_20260508.xlsx` as the frozen trustworthy baseline until a candidate version beats it on documented evidence, reproducible validation, and adversarial review.

Keep iterating. Do not end with a plan, a single report, or a partially tested candidate unless an external blocker makes further work impossible. A candidate version can replace v1.13 only after repeated build-test-criticize-rebuild cycles satisfy the acceptance gates in `references/top_journal_strategy.md`.

Do not tune one reaction to satisfy one phenotype. Improve the whole model only when the change survives these gates:

1. Literature or experimental source is cited.
2. Stoichiometry, formula, charge, GPR, and bounds are audited.
3. The existing phenotype panel is rerun.
4. MEMOTE/COBRA checks are rerun.
5. A counter-hypothesis is tested and recorded.
6. The change is reversible and versioned.

## Required Workflow

1. Confirm baseline files in `baseline/` and read `baseline/CLAUDE_HANDOFF_V1.13_20260509.md`.
2. Load `references/top_journal_strategy.md` for the full experimental/modeling roadmap.
3. Load `evidence/LITERATURE_READING_PROTOCOL.md` and `evidence/literature_evidence.tsv` before adding constraints or manuscript claims.
4. Load `evidence/phenotype_data_requirements.tsv` before requesting new wet-lab data.
5. Run `scripts/check_release_files.py` after edits to confirm the handoff package still has the required files.

Before editing the model, create a candidate-specific `evidence_extraction.tsv` using the extraction schema in `evidence/LITERATURE_READING_PROTOCOL.md`.

## Modeling Priorities

Priority 1: add enzyme-capacity constraints.

- Use absolute or calibrated protein abundances for condition-specific enzyme upper bounds.
- Use BRENDA/SABIO-RK/primary literature kcat values; keep organism, enzyme, temperature, pH, substrate, and confidence tier.
- Apply first to CBB cycle, nitrogenases, photosynthetic reaction center/electron transport, acetate assimilation, aromatic degradation, PHB, TCA/glyoxylate, ATP maintenance.

Priority 2: add thermodynamic constraints.

- Use metabolite concentration ranges only when measured for R. palustris or justified by close context.
- Use component-contribution/eQuilibrator-style standard transformed Gibbs energies.
- Avoid one-off direction edits when a pathway-wide redox/proton convention is unresolved.

Priority 3: validate against 13C-MFA.

- Extract flux ratios and absolute exchange rates from McKinlay/Harwood and related R. palustris studies.
- Compare internal flux distributions, not only growth/no-growth.
- Treat CBB flux, TCA/glyoxylate branch, acetate assimilation route, aromatic funneling, nitrogenase/H2, and PHB as the primary validation axes.

## Adversarial Review Checklist

For every candidate change, write a short record:

- Claim: what improves?
- Evidence: source and exact data used.
- Alternative: what else could explain the phenotype?
- Failure test: what result would falsify the change?
- Scope: which conditions are affected?
- Regression: did any existing validated phenotype or QC metric degrade?

Reject changes that improve a single headline metric but reduce biological credibility.
