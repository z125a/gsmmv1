# Candidate v1.14 Iteration Log

## Cycle 1

Hypothesis: Before editing reactions, v1.13 should be reproducible through a Python-only path and the first evidence ledger should identify safe next constraints.

Actions:

- Cloned `https://github.com/yms-xjtu/gsmmv1.git`.
- Read the next-agent prompt, bundled skill, literature protocol, evidence ledger, phenotype data requirements, baseline handoff, and top-journal strategy.
- Ran the package integrity check successfully.
- Ran Python/COBRApy baseline reproduction with GLPK.
- Created first-pass evidence extraction, data-gap report, and rejected-shortcuts log.

Decision:

- Keep the cycle outputs.
- Do not edit the model yet. The evidence gate is now started, but hard ecGEM/TFA/13C-MFA constraints require more extraction from full papers, supplements, BRENDA/SABIO-RK, and metabolomics sources.

Next falsifiable hypothesis:

The McKinlay/Harwood 13C flux paper plus Navid/iRpa940 conditions can define a first internal-flux validation panel for acetate, succinate, and butyrate without changing v1.13 stoichiometry.
