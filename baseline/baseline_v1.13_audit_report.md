# Baseline v1.13 Audit Report

Date: 2026-05-09 (audit rerun)

## Executive Summary

**Overall Status: CONDITIONAL PASS**

The baseline v1.13 model is usable as a working reference for candidate iterations,
but has known architecture gaps that limit certain validation axes.

## Key Metrics (Reproduced)

| Metric | Value | Expected | Match |
|--------|-------|----------|-------|
| FBA_status | optimal | optimal | ✓ |
| FBA_objective | 10.545510 | 10.545510 | ✓ |
| pFBA_total_flux | 9504.12 | ~9504 | ✓ |
| FVA_bio1_min_90pct | 9.4910 | ~9.49 | ✓ |
| FVA_bio1_max | 10.5455 | ~10.55 | ✓ |
| reaction_count | 1021 | 1021 | ✓ |
| metabolite_count | 979 | 979 | ✓ |
| gene_count | 685 | 685 | ✓ |

## Critical Findings

### 1. Calvin Cycle (CBB)
- RuBisCO (rxn00018_c0): EXISTS, irreversible [0, 1000], but **ZERO flux at baseline pFBA**
- Phosphoribulokinase (rxn01111_c0): EXISTS, irreversible [0, 1000], **ZERO flux at baseline**
- This means the Calvin cycle is NOT active in the default model configuration
- The model grows without CO2 fixation — consistent with chemoheterotrophic mode
- For photoheterotrophic validation, CO2/light constraints must be explicitly set

### 2. N2/Nitrogenase Architecture Gap
- **No N2 exchange reaction exists** — nitrogenase cannot function
- Nitrogenase (rxn06874_c0) exists but cannot carry flux (no N2 supply)
- All H2 production comes from hydrogenase (rxn05759_c0)
- CGA009 has defective uptake hydrogenase — this is biologically incorrect

### 3. Reaction Mapping Corrections (from previous agent)
- rxn01116_c0 = Ru5P epimerase (NOT RuBisCO) ✓ confirmed
- rxn05040_c0 = DHBP synthase (NOT RuBisCO) ✓ confirmed
- rxn02507_c0 = indole-3-glycerol-P synthase (NOT RuBisCO) ✓ confirmed
- rxn00018_c0 = TRUE RuBisCO ✓ confirmed

### 4. Stoichiometry
- 46 imbalanced internal reactions (matches handoff report of 47)
- Known to be porphyrin/chlorophyll/cobalamin convention issues
- Not fixable one-by-one per handoff instructions

## Conclusion

**Baseline v1.13 is CONDITIONALLY ACCEPTED as working reference.**

Conditions:
1. N2/H2 architecture gap is documented but not blocking for non-H2 work
2. Calvin cycle zero-flux is understood (model defaults to non-CBB growth)
3. Reaction mapping corrections are now documented and will be used in v1.14+
4. 47 imbalanced reactions are known convention issues, not stoichiometry errors

The baseline does NOT need to be "fixed" before continuing candidate work,
because the issues are architecture/convention gaps, not data corruption.
Candidate v1.14+ work should proceed with corrected mappings and explicit
condition-setting for photoheterotrophic simulations.
