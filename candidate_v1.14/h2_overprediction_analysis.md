# H2 Over-Prediction Root Cause Analysis

## Problem

Model predicts H2 yields 10-22x higher than experimental values from McKinlay/Harwood 2011.

## Root Cause

The model has no constraint on the rate of reduced ferredoxin (Fd_red) production. When H2 is maximized:

1. POR (rxn13974_c0) runs at high flux, producing Fd_red from pyruvate oxidation
2. 2-oxoglutarate:Fd oxidoreductase (rxn14048_c0) runs in reverse, producing Fd_red
3. Fd:NADP+ oxidoreductase (rxn14159_c0) runs in reverse, converting NADPH → Fd_red
4. All Fd_red is consumed by nitrogenase (rxn06874_c0, 8 Fd_red per H2)

In reality, Fd_red production is limited by:
- Photosynthetic electron transport (cyclic photophosphorylation)
- The rate of electron flow through PSI → ferredoxin
- Chowdhury 2022 identifies this as the ETFD (electron transport through ferredoxin) bottleneck

## Why ATP Cap Doesn't Help

ATP cap only limits ATP synthase. But H2 production via nitrogenase requires both ATP AND Fd_red. The model can generate Fd_red through substrate oxidation (POR, TCA) independently of photosynthetic ATP, so capping ATP barely affects H2 yield.

## Evidence-Supported Solutions (for future cycles)

1. **Ferredoxin production rate constraint** (Priority 1):
   - Limit total Fd_red production to match photosynthetic electron transport capacity
   - Evidence: Chowdhury 2022 ETFD analysis; McKinlay/Harwood electron balance (Fig 2)
   - Implementation: Add upper bound on net Fd_red production reactions

2. **Electron balance constraint** (Priority 2):
   - McKinlay/Harwood Fig 2 shows electron balance within 10% for all substrates
   - Total electron carrier reduction ≈ total electron carrier oxidation
   - Calvin cycle oxidizes 38-55% of carriers; biosynthesis uses the rest

3. **Nitrogenase capacity constraint** (Priority 3):
   - From ME model: nitrogenase activity is coupled to growth via μ/kcat
   - Need kcat for nitrogenase (Mo-Nase) from literature

## What NOT To Do

- Do NOT simply cap nitrogenase flux to match experimental H2 — this is circular
- Do NOT remove nitrogenase from the model
- Do NOT add artificial electron sinks

## Status

This is a **known limitation** of stoichiometric models without enzyme-capacity constraints.
The ecGEM layer (with proper Fd_red production limits) is the correct solution path.
Recording as a data gap requiring either:
- Measured Fd_red production rates under photoheterotrophic conditions
- Or: calibrated photosynthetic electron transport rate from the ME model

## Ferredoxin Reactions in v1.13

| Reaction | Name | Fd_red role | Baseline flux |
|----------|------|-------------|---------------|
| rxn13974_c0 | POR | produces 2 Fd_red | 89.05 |
| rxn14048_c0 | 2-OG:Fd oxidoreductase | produces 2 Fd_red (reverse) | -55.27 |
| rxn14159_c0 | Fd:NADP+ oxidoreductase | produces 2 Fd_red (reverse) | -34.28 |
| rxn06874_c0 | Nitrogenase | consumes 8 Fd_red | 0.00 (WT) |
| rxn05759_c0 | Hydrogenase | consumes 2 Fd_red | 0.00 |
| rxn14070_c0 | HMB-PP reductase | consumes 2 Fd_red | 0.51 |
