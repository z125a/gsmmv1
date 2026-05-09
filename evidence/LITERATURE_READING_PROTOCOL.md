# Mandatory Literature Reading Protocol

The next agent must read and extract evidence before editing the model.

## What Is Included

This repository includes a machine-readable literature ledger:

- `evidence/literature_evidence.tsv`

The ledger gives prioritized papers, databases, URLs/DOIs, expected data to extract, and how each source should be used for CGA009/TX73 model improvement.

Full article PDFs are not bundled by default. Use the URLs/DOIs in the ledger to retrieve the latest accessible versions, supplementary files, and database records. Do not rely on memory or uncited claims.

## Required Reading Order

1. CGA009 ME-model paper:
   - Use for Rubisco/nitrogenase/enzyme-capacity precedent.
   - Extract model assumptions, enzyme constraints, growth conditions, kcat sensitivity, and GPR/protein mappings.

2. McKinlay/Harwood Calvin cycle and 13C/redox papers:
   - Use for internal flux validation, CBB cycle requirements, H2 yield, redox balancing, and isotope-informed constraints.

3. Navid et al. iRpa940 paper and local supplementary model:
   - Use for reference reaction content, photoheterotrophic tradeoff modeling, and comparison to the current v1.13 model.

4. BisA53 GEM paper:
   - Use as a benchmark for breadth of phenotype validation and reviewer-grade model evaluation.
   - Do not transfer strain-specific claims directly to CGA009/TX73.

5. CGA009 quantitative proteomics / aromatic growth paper:
   - Use for condition-specific protein abundance priors and growth phenotypes.
   - Treat LFQ or relative abundance as priors unless absolute calibration is available.

6. TFA, GECKO, BRENDA, SABIO-RK, eQuilibrator, MEMOTE:
   - Use for method implementation and data provenance.

## Extraction Template

For every source, create or append to `candidate_v*/evidence_extraction.tsv` with these columns:

- source_id
- url_or_doi
- exact_claim
- organism_or_strain
- growth_condition
- substrate
- oxygen_light_nitrogen_regime
- measured_quantity
- value
- unit
- uncertainty_or_range
- model_reaction_or_gene
- use_type: measured_constraint | literature_constraint | validation_target | sensitivity_prior | background_only
- confidence_tier: high | medium | low
- reason_for_tier
- model_action_allowed: yes | no
- notes

## Gate Before Model Editing

Do not change the model until at least one evidence row justifies the proposed change.

For ecGEM/ME constraints, every enzyme constraint must cite abundance and kcat provenance separately.

For TFA, every hard thermodynamic direction/concentration constraint must cite Gibbs energy and concentration assumptions.

For 13C-MFA, every internal flux validation target must cite the isotope experiment, growth condition, and reported flux or flux ratio.

If a source is unavailable, record it in `candidate_v*/missing_literature_or_data.tsv` and continue with another evidence-supported route.
