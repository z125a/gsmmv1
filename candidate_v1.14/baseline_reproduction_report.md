# Candidate v1.14 Baseline Reproduction Report

Date: 2026-05-09

## Purpose

Reproduce the frozen v1.13 baseline using the Python/COBRApy path before making any model edits.

## Input Model

`baseline/mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_for_memote.xml`

## Execution

Command:

```powershell
E:\model\memote_py312\Scripts\python.exe skills\cga009-top-journal-modeling\scripts\python_cobrapy_baseline_template.py
```

The script forces COBRApy to use the open `glpk` solver so the workflow does not depend on MATLAB or Gurobi.

## Outputs

`candidate_v1.14/python_cobrapy_baseline/`

## Reproduced Python/COBRApy Baseline

| Metric | Value |
| --- | --- |
| model_id | COBRAModel |
| reactions | 1021 |
| metabolites | 979 |
| genes | 685 |
| solver | optlang GLPK interface |
| objective | `bio1` |
| FBA status | optimal |
| FBA objective value | 10.545509652968153 |
| FVA reaction checked | `bio1` |
| FVA minimum at 90% optimum | 9.490958687688606 |
| FVA maximum | 10.545509652968134 |

## Interpretation

The Python-only route can load and optimize the v1.13 SBML baseline without MATLAB. This satisfies the first portability requirement for the next iteration.

This report does not claim v1.14 is improved. No model edit has been made yet because the literature evidence gate must be satisfied first.
