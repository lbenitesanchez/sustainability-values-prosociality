# Regenerate Figure 1 (M2) without confidential inputs

From the root of this package:

```r
install.packages(c("lavaan", "digest")) # only if not already installed
source("R/02_draw_M2.R", encoding = "UTF-8")
```

Or use a clean process:

```sh
Rscript --vanilla R/02_draw_M2.R
```

The script finds the package relative to its own path and writes to `figures/regenerated/`. Set `M2_PUBLIC_INPUT_DIR` or `M2_FIGURE_DIR` to override input/output folders.

## Public inputs

- `inputs/M2_parameter_manifest.csv`: parameter table exported from the archived M2 SEM, including fixed/free status and labels. Contains aggregate model results, not respondent records.
- `inputs/M2_model_metadata.csv`: N, estimator/robust standard error settings, and typed logical identification options.
- `inputs/M2_puentes2_syntax.txt`: actual M2 syntax.

No SPSS file, respondent-level dataset, `.RData`, or fitted `.rds` is required for this mode. The public script verifies 21 focal regressions, 49 control regressions, 24 loadings, seven free residual covariances, and three marker loadings fixed to 1. The four values remain observed variables; PSB remains latent with PP1–PP3 indicators.

## Scope and outputs

This is a **regeneration of the figure from exported model structure**, not re-estimation or raw-data reproduction. `lavaanify()` parses syntax; it does not fit a model. The output is a vector PDF, a vector SVG, a 7800 × 5700 PNG with 600 dpi metadata, captions, checks, hashes, and `sessionInfo()`.

The SVG and both PNG files were byte-identical when regenerated using the public inputs and when extracted from the private archived model on the recorded environment. `public_regeneration_comparison.csv` records the hashes. PDF creation timestamps may differ although the drawing is the same. Font rendering can differ across operating systems; diagram coordinates do not.

The seven controls and residual covariance structure are documented in the caption, not drawn as extra nodes. Complete numerical loadings and structural results are in the supplement and exported CSV tables. Arrow colors/styles distinguish path families, not significance, causality or time order.
