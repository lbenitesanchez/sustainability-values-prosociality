# Sustainability attitudes, human values and prosocial self-descriptions among business students: An exploratory secondary analysis

Materials accompanying the expanded M2 manuscript and supplement. The package contains **four main tables, 29 supplementary tables (some divided into panels), and one model figure**. The public analytical materials are available in this repository: https://github.com/lbenitesanchez/sustainability-values-prosociality. Original participant records and private fitted objects are not included.

## Start here

Run from this directory in a clean R session:

```sh
Rscript --vanilla R/00_master.R --mode=check_exports
```

This checks the table inventory, generator paths and preserved source files. It does **not** fit models. See [the execution guide](documentation/README_RERUN.md) for private model extraction and artificial-data modes, the execution order and required private inputs.

The complete mapping is [documentation/table_map.csv](documentation/table_map.csv). S28 has four panels and S29 two; therefore 33 logical tables are stored in 37 formatted CSVs. [documentation/figure_map.json](documentation/figure_map.json) describes Figure 1.

## Contents

- `manuscript/`: the expanded English Word documents.
- `tables/`: the exact formatted tables and panels used in Word.
- `results/complete/`: full-precision aggregate parameters, sample summaries and numerical checks.
- `data/`: the original-variable and analytical-control dictionaries, data-access instructions, historical summary matrices and explicitly artificial demonstration data.
- `figures/`: vector PDF/SVG, 600-dpi PNG, captions and public structural inputs for Figure 1.
- `R/`: documented extraction, figure, synthetic and preserved analytical scripts.
- `documentation/`: table mapping, source hashes, environment information and scope of verified reproduction.

To redraw the figure from the public manifest, without participant data or fitted RDS:

```sh
Rscript --vanilla R/02_draw_M2.R
```

Check that script's output-directory option before rerunning if the supplied figure files are to remain untouched.

## What was verified

The added outputs were extracted from the existing M2 CFA, factor-PSB SEM and individual-response fits. Sample descriptives were calculated from the original labelled SPSS file. The existing 24 tables were compared with their frozen references; the additional tables were checked against their own fitted objects. No model, exclusion rule, indicator allocation or original multiplicity family was changed. The latest update did not re-estimate SEM/CFA or retrain prediction models.

The full historical fitting workflow is preserved as source code but is **not independently executable using only these public aggregates**. It requires authorised individual records and private dependencies described in `data/README_ACCESS.md` and the execution guide. Original SAV, private RDS/RData, individual predictions and resampling assignments are excluded.

## Interpretation

The analysis uses the same 653 questionnaires as the source ASP-S study. One institution represents each country; institutional and country effects cannot be separated. Estimates concern cross-sectional associations and self-descriptions, not observed behaviour or causal effects. Regression Holm adjustments, pointwise intervals and the original standardised indirect-effect Holm family are identified separately. Artificial examples do not reproduce empirical estimates.
