# Calvin Cycle Zero-Flux Investigation

## Finding

RuBisCO (rxn00018_c0) and PRK (rxn01111_c0) have ZERO flux at baseline pFBA,
even when hydrogenase is blocked and H2 export is disabled.

## Root Cause

The model contains internal thermodynamic loops that act as free electron sinks,
making the Calvin cycle unnecessary for redox balance:

1. **rxn00501_c0** (3-oxopropanoate:NAD+ oxidoreductase): reversible, carries -112 NADH
   - Running in reverse, consuming NADH without biological purpose
   - Part of a loop with malonyl-semialdehyde metabolism

2. **rxn00929_c0** (pyrroline-5-carboxylate reductase): reversible, carries -78 NADH
   - Running in reverse, consuming NADH
   - Part of a proline/glutamate cycling loop

3. **rxn00248_c0** (malate dehydrogenase): carries +103 NADH production
   - Feeds NADH into the loop sinks above

These loops collectively dissipate ~190 NADH units — more than enough to replace
the Calvin cycle's role as electron sink (which would only handle ~40-55% of
electron carriers per McKinlay/Harwood).

## Why This Matters

- McKinlay/Harwood 2011 showed Calvin cycle is ESSENTIAL for photoheterotrophic growth
- ΔRuBisCO mutants cannot grow without nitrogenase/H2 as alternative electron sink
- The model contradicts this because internal loops provide unlimited electron disposal
- This means ALL previous "Calvin cycle active" conclusions were artifacts of
  misidentified reactions (Ru5P epimerase flux, not actual CO2 fixation)

## FVA Evidence

- RuBisCO FVA: min=0, max=12.62 (CAN carry flux but doesn't NEED to)
- Forcing RuBisCO ≥ 0.1: feasible, growth barely changes (10.545 → 10.545)
- This confirms Calvin cycle is optional in the current model

## Correct Fix (for v1.15)

The fix is NOT to force RuBisCO flux. The fix is to eliminate the
thermodynamically infeasible loops that provide free electron disposal:

1. Apply loopless FBA constraints (if available in COBRApy)
2. Or: apply TFA direction constraints to the loop-forming reactions
3. Or: identify and block the specific loop cycles

Once loops are eliminated, the model SHOULD require Calvin cycle for
photoheterotrophic growth, matching biological reality.

## Verification Plan

After loop elimination:
1. WT on acetate without H2: should REQUIRE Calvin cycle (growth = 0 without RuBisCO)
2. NifA* on acetate with H2: Calvin cycle should decrease, H2 should increase
3. ΔRuBisCO mutant: should only grow with nitrogenase active
