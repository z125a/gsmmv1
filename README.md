# gsmmv1

# CGA009/TX73 GEM v1.13 GitHub Release Handoff

This package freezes the current trustworthy baseline and gives the next agent a rigorous path toward a stronger model.

## Trusted Baseline

Use `baseline/mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_20260508.xlsx` as the current trusted workbook.

Why v1.13:

- It is the latest validated workbook found in the project.
- It supersedes v1.12.
- MATLAB/COBRA QC reports 1021 reactions, 979 metabolites, 945 internal reactions, 47 imbalanced internal reactions, and 0 metabolites missing formula.
- Current phenotype panel remains intact with expected-design pass rate 100%.
- MEMOTE annotated total score is 67.8836%.

## What This Package Adds

- `baseline/`: copied v1.13 model, SBML, QC summaries, phenotype validation, and handoff notes.
- `evidence/literature_evidence.tsv`: prioritized papers/databases for ecGEM, TFA, and 13C-MFA.
- `evidence/LITERATURE_READING_PROTOCOL.md`: mandatory literature-reading and evidence-extraction protocol before model editing.
- `evidence/phenotype_data_requirements.tsv`: wet-lab and literature data needed to raise model credibility.
- `next_agent/NEXT_AGENT_PROMPT.md`: direct prompt for the next person/agent.
- `skills/cga009-top-journal-modeling/`: a bundled Codex skill for adversarial, publication-grade model improvement.

## Run A Package Check

```powershell
python .\skills\cga009-top-journal-modeling\scripts\check_release_files.py
```

## Execution Environment

Use Python/COBRApy as the default execution path. MATLAB/COBRA Toolbox is optional and should only be used when available. A next agent without MATLAB must still continue using Python tools such as COBRApy, memote, pandas, scipy, optlang, GLPK/HiGHS-compatible solvers, and notebook or script-based reports.

The agent should work autonomously: do not ask for repeated yes/no confirmation during normal literature extraction, script writing, validation, candidate-version creation, or Git commits. Ask the user only when there is a genuine external blocker, destructive action, credential need, paid/proprietary data access issue, or an irreversible scientific decision.

## Next Scientific Direction

Do not continue one-reaction patching as the main strategy. The next credible improvement is a whole-model validation program:

1. ecGEM/ME-model layer using protein abundance, kcat, and enzyme molecular weight.
2. TFA layer using metabolite concentration bounds and Gibbs energy constraints.
3. 13C-MFA layer using isotope-derived internal fluxes as the strict validation target.

Each candidate version must be evidence-tracked, reproducible, and adversarially reviewed before replacing v1.13.
