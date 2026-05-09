# v1.14 Final Status Report

## Executive Summary

v1.14 introduces ONE validated improvement over baseline v1.13:
- **ICL=0 on succinate** (evidence: McKinlay/Harwood 2011 13C-MFA)
- Growth impact: -0.28% (negligible)
- MFA qualitative tests: 15 PASS / 0 FAIL / 0 INFO

## Corrected Mappings (vs previous agent errors)

| Previous claim | Correction |
|---------------|------------|
| rxn01116 = RuBisCO | rxn01116 = Ru5P epimerase |
| rxn05040 = RuBisCO | rxn05040 = DHBP synthase |
| rxn02507 = RuBisCO | rxn02507 = indole-3-glycerol-P synthase |
| "Calvin cycle active" | Calvin cycle has ZERO flux at baseline |

## Current Blockers

1. **No N2 exchange** → nitrogenase cannot function → H2 validation impossible
2. **RuBisCO/PRK zero flux** → Calvin cycle not active in default config
3. **Hydrogenase produces H2** → biologically incorrect for CGA009

## Recommendation

**Current best model: baseline v1.13 + ICL constraint on succinate**

v1.14 is a marginal improvement (one validated constraint) but does not yet
warrant a new SBML release. Continue to v1.15 with architecture fixes.
