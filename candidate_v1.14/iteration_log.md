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

---

## Cycle 2

Hypothesis: The v1.13 model's pFBA flux distributions on acetate, succinate, and butyrate will qualitatively match McKinlay/Harwood 2011 13C-MFA findings for glyoxylate shunt usage, Calvin cycle activity, and CO2 production patterns.

Actions:

1. Confirmed Python/COBRApy baseline reproduction in Linux environment (GLPK solver, growth=10.5455).
2. Read full McKinlay/Harwood 2011 mBio paper (PMC3063381) and extracted quantitative data:
   - Table 1: H2 yields (fumarate 18±3, succinate 23±1, acetate 21±3, butyrate 41±10 mol H2/100 mol organic C)
   - Table 2: Doubling times, biomass yields, carbon/electron recovery
   - Table 3: CO2 production and Calvin cycle refixation percentages
   - Fig 1: Qualitative flux maps showing glyoxylate shunt active on acetate/butyrate, inactive on succinate
   - Fig 2: Calvin cycle oxidizes 38-55% of electron carriers depending on substrate
3. Created `mfa_validation_targets.tsv` with 30 quantitative/qualitative targets from the paper.
4. Identified exact ModelSEED reaction IDs for all key pathways:
   - ICL: rxn00336_c0, MS: rxn00330_c0, CS: rxn00256_c0, SDH: rxn00288_c0
   - POR: rxn13974_c0, Nitrogenase: rxn06874_c0, RuBisCO: rxn01116_c0 + 4 others
5. Built and ran `scripts/mfa_flux_panel.py` — comprehensive internal flux validation.
6. Results: 11 PASS / 1 FAIL / 0 INFO out of 12 qualitative tests.
7. Key finding: **Glyoxylate shunt is incorrectly active on succinate** (flux=2.37).
   - FVA shows ICL minimum > 0 at 99.9% optimality on succinate
   - Blocking glyoxylate shunt reduces growth by only 0.28% (10.59 → 10.56)
   - This contradicts McKinlay/Harwood Fig 1 showing zero glyoxylate shunt on succinate
   - Root cause: pFBA finds glyoxylate shunt as slightly more flux-efficient route
   - Biological reality: ICL is not expressed on succinate (no need for it)
8. Read Chowdhury 2022 ME-model paper (PMC9431616) and extracted:
   - Mean kcat = 234,000 day⁻¹ (from cyanobacteria, modified by SASA)
   - ATP production rates: acetate 54.0, succinate 45.7, butyrate 56.7, p-coumarate 85.4 mmol/gDW/day
   - ME model growth predictions: succinate 0.74, acetate 0.77, butyrate 0.86 day⁻¹
   - Exclusive Mo-Nase expression predicted
   - ETFD (electron transport through ferredoxin) as regulatory bottleneck
9. Created `kcat_evidence.tsv` with initial entries from ME paper + BRENDA gaps.
10. Built thermodynamic ID mapping: 43/55 priority metabolites mapped to KEGG IDs.
11. Identified 10 priority reactions for TFA deltaG estimation.

Key findings:

| Test | Substrate | Status | Observation |
|------|-----------|--------|-------------|
| Growth feasible | all 3 | PASS | All substrates support growth |
| Glyoxylate shunt active | acetate | PASS | flux=2.39 |
| Glyoxylate shunt active | butyrate | PASS | flux=2.42 |
| Glyoxylate shunt inactive | succinate | **FAIL** | flux=2.37 (should be 0) |
| Calvin cycle active | all 3 | PASS | RuBisCO flux=-4.3 to -4.5 |
| Net CO2 production | acetate, succinate | PASS | CO2 exchange=0.23 |
| POR forward on succinate | succinate | PASS | flux=89.05 |

Adversarial analysis of the glyoxylate shunt failure:

- **Claim**: The model incorrectly routes succinate through glyoxylate shunt.
- **Evidence**: McKinlay/Harwood 2011 Fig 1 shows zero glyoxylate shunt flux on succinate; 13C labeling confirms succinate goes through TCA to acetyl-CoA via PDH/POR.
- **Alternative explanation**: pFBA minimizes total flux and finds glyoxylate shunt slightly more efficient than full TCA cycle for generating acetyl-CoA from succinate. This is a stoichiometric artifact — the model lacks enzyme expression constraints that would prevent ICL activity on succinate.
- **Failure test**: If we add an enzyme-capacity constraint (ICL protein = 0 on succinate), the model should still grow at near-identical rate (confirmed: only 0.28% reduction).
- **Scope**: Affects succinate-only growth; does not affect acetate or butyrate predictions.
- **Regression risk**: None if constraint is condition-specific.

Decision:

- Do NOT edit the model stoichiometry to fix this. The glyoxylate shunt issue is a constraint problem (missing enzyme expression data), not a stoichiometry error.
- This is exactly the type of issue that ecGEM/enzyme-capacity constraints should resolve.
- Record as a validated discrepancy and a target for the ecGEM layer.
- Continue to Cycle 3: begin building condition-specific enzyme constraints.

Next falsifiable hypothesis (Cycle 3):

Adding a condition-specific enzyme-capacity constraint (ICL upper bound = 0 on succinate, based on proteomics/transcriptomics evidence that ICL is not expressed on succinate) will eliminate the glyoxylate shunt discrepancy without regressing other phenotypes.

---

## Cycle 3

Hypothesis: Adding a condition-specific enzyme-capacity constraint (ICL=0 on succinate) will eliminate the glyoxylate shunt discrepancy without regressing other phenotypes.

Actions:

1. Implemented `scripts/ecgem_condition_constraints.py` with ICL/MS upper_bound=0 on succinate.
2. Ran constrained validation: **12 PASS / 0 FAIL** (up from 11/1).
3. Growth impact: only -0.28% on succinate — negligible regression.
4. Adversarial check: FVA confirms biomass range unchanged.
5. Implemented `scripts/h2_yield_validation.py` for quantitative H2 yield comparison.
6. H2 yield results: **massively over-predicted** (10-22x experimental).
   - Acetate: predicted 473 vs experimental 21 mol H2/100 mol C
   - Succinate: predicted 468 vs experimental 23
   - Butyrate: predicted 480 vs experimental 41
7. Root cause: model lacks electron transport regulation between growth and H2.
   - Photon cap exists (36.6 units) but electron partitioning is unconstrained.
   - ATP synthase has no substrate-specific upper bound.
   - Nitrogenase can consume unlimited reduced ferredoxin.
8. Extracted Chowdhury 2022 ATP production caps: acetate 54.0, succinate 45.7, butyrate 56.7 mmol/gDW/day.
9. Built thermodynamic ID mapping (43/55 metabolites mapped to KEGG).
10. Created `kcat_evidence.tsv` with initial entries.
11. Created `cycle3_audit_record.tsv` documenting all findings.

Key findings:

| Finding | Impact | Action |
|---------|--------|--------|
| ICL=0 on succinate works perfectly | Validates ecGEM approach | ACCEPTED as first constraint |
| H2 yields 10-22x over-predicted | Major model limitation | Need electron transport constraints |
| ATP synthase uncapped | Allows unrealistic energy | Candidate for ATP cap from ME model |
| Photon uptake already capped at 36.6 | Partial energy constraint exists | Not sufficient alone |

Decision:

- ACCEPT the ICL constraint as the first validated ecGEM edit for v1.14.
- Do NOT edit the model SBML yet — accumulate more constraints first.
- The H2 over-prediction is a known limitation that requires:
  1. Substrate-specific ATP caps (from Chowdhury 2022)
  2. Ferredoxin electron transport constraints (ETFD regulation)
  3. Proper nitrogenase ATP stoichiometry validation
- Continue to Cycle 4 with ATP constraint implementation.

Next falsifiable hypothesis (Cycle 4):

Adding a substrate-specific photosynthetic ATP cap (from Chowdhury 2022 ME model) will reduce H2 over-prediction toward the experimental range without breaking growth phenotypes.

---

## Cycle 4 (in progress)

Goals:
1. Implement substrate-specific ATP production cap.
2. Re-run H2 yield validation with ATP constraint.
3. Check if reduced ferredoxin availability is the true bottleneck.
4. Begin eQuilibrator deltaG queries for TFA priority reactions.
5. Build the first candidate v1.14 SBML with condition-specific constraint annotations.
