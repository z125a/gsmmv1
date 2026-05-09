# Top-Journal Strategy For CGA009/TX73 GEM

## Baseline

Use v1.13 as the current credible version. It has 1021 reactions, 979 metabolites, no metabolites missing formula in MATLAB/COBRA QC, 47 imbalanced internal reactions, 100% expected-design phenotype pass rate in the current panel, and MEMOTE total score 67.8836%.

The next version must not merely increase MEMOTE. It must increase predictive trust by reducing under-constrained FBA behavior and validating internal fluxes.

## Execution Policy

Use Python/COBRApy as the primary, portable execution path. MATLAB/COBRA Toolbox results from v1.13 are trusted baseline evidence, but future agents must not require MATLAB to continue work. If MATLAB is unavailable, reproduce and extend validation using Python scripts, COBRApy, memote CLI/Python, pandas, scipy, optlang, and available open solvers.

The agent should proceed without repeated user confirmations. Continue literature extraction, evidence-table creation, candidate model construction, validation, rejection logging, and Git commits autonomously. Stop to ask only for credentials, paid/proprietary data, destructive overwrite approval, or a curator-level biological convention that cannot be resolved from evidence.

## Three Evidence Layers

### 1. ecGEM / ME-model constraints

Goal: limit flux by enzyme capacity, not only stoichiometric feasibility.

Required data:

- Condition-resolved absolute or calibrated protein abundance for CGA009/TX73.
- Carbon sources: acetate, succinate, malate, butyrate, p-coumarate, benzoate or other aromatic substrates.
- Regimes: light anaerobic, dark aerobic, nitrogen-fixing versus ammonium, varied light intensity.
- kcat and molecular weight table for core enzymes.

Minimum implementation:

- Map genes to proteins to reactions through GPRs.
- Convert protein abundance to enzyme capacity.
- Add v <= kcat * E constraints for high-confidence enzymes.
- Run sensitivity bands for uncertain kcat values instead of using a single false-precise value.

### 2. TFA constraints

Goal: eliminate thermodynamically impossible loops and unrealistic reaction directions.

Required data:

- Intracellular concentration ranges for central metabolites under matched conditions.
- pH, ionic strength, temperature, and compartment assumptions.
- Standard transformed Gibbs energies and uncertainty.

Minimum implementation:

- Apply TFA first to central carbon, redox carrier, CBB, TCA/glyoxylate, PHB, aromatic degradation, and nitrogenase-adjacent electron transfer modules.
- Flag reactions whose feasibility depends entirely on broad concentration bounds.
- Keep unresolved pseudo-metabolite and carrier conventions out of hard TFA until the whole module is curated.

### 3. 13C-MFA validation

Goal: test whether internal flux predictions match isotope-derived flux distributions.

Required data:

- Labeling substrate, enrichment, sampling time, growth regime, exchange fluxes, and reported flux ratios.
- Core fluxes around CBB, TCA/glyoxylate, acetyl-CoA assimilation, pyruvate/PEP/OAA nodes, PHB, CO2 fixation/release, and H2.

Minimum implementation:

- Reproduce the experimental medium and exchange bounds before comparing fluxes.
- Compare flux ratios and branch usage, not only objective value.
- Use pFBA/FVA/loopless FBA as model-side distributions, then test whether MFA fluxes lie inside constrained ranges.

## Acceptance Criteria

A candidate release can supersede v1.13 only if it has:

- All baseline phenotype tests passing or justified with new experimental evidence.
- No unexplained increase in internal imbalanced reactions.
- MEMOTE and COBRA reports saved.
- A table of all changed reactions/metabolites with source evidence.
- ecGEM/TFA/MFA constraints separated into measured, literature-derived, and sensitivity-only tiers.
- An adversarial review note listing failed hypotheses and rejected shortcuts.

If these criteria are not met, continue iterating. Do not stop at "next steps" or "recommendations." Create the next candidate, run the checks again, and keep the rejected-hypothesis log current until the model can be defended as publication-ready.

## Iteration Loop

Each cycle must produce artifacts:

1. Reproduce the previous baseline metrics.
2. Form one falsifiable model-improvement hypothesis.
3. Implement the smallest evidence-supported candidate change set.
4. Run phenotype, COBRA, MEMOTE, pFBA/FVA, loopless/TFA-relevant checks, and any available MFA comparison.
5. Write the strongest counterargument and test it.
6. Keep or reject the candidate.
7. If rejected, preserve the failure record and start the next hypothesis.
8. If kept, version the candidate and continue until the publication-ready gates are satisfied.

## Stop Conditions

Do not edit the porphyrin/chlorophyll/cobalamin module one reaction at a time. Solve formula/proton/charge conventions pathway-wide.

Do not edit tRNA(Gln) pseudo-metabolites without a curator decision.

Do not add exchange or transport reactions just to force growth.

Do not report optimized PHB or H2 as experimental truth without measured uptake/secretion or 13C-MFA support.
