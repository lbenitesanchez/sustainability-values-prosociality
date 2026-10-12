# R package validation record

The public R entry point was tested on 11 October 2026 using R 4.6.1. These checks concern execution and artifact consistency, not scientific model adequacy or full raw-data reproduction.

- All 40 R files then present parsed successfully.
- The 39 preserved original/extension source files matched their indexed hashes; no duplicate source contents were included.
- Synthetic mode produced 800 artificial rows and 49 columns and passed all 17 built-in checks. It read no participant data and fitted no model.
- Private-extraction mode ran in a fresh `Rscript --vanilla` process with explicit private project/reference/Word/library paths. The historical 24-table stage passed 48 comparisons and 2,710 cell comparisons, with zero unexplained differences and two documented M0 label changes. It fitted zero new models.
- The complete-M2 stage extracted existing fitted objects and calculated sample aggregates; the assembled output contained 33 canonical table IDs and 37 CSV table/panel files. This validates the private extraction route, not a rerun of original SEM/CFA/CV fitting.
- No SAV, RDS, RData, workbook, ZIP, or participant-level file was added under public `R/` or `documentation/` during preparation.

The final `check_exports` default passed against the enriched 33-ID/37-CSV inventory, including all mapped generator paths and all 39 source-hash checks. The latest private-extraction run also passed. Comparing its 37 resulting table/panel CSVs with the public tables found no numerical changes; only three S22 cells intentionally shorten the repeated CSR1 row label from "CSR1: attended ethics/CSR course at university" to "CSR1 course attendance". Layout whitespace was normalized for this comparison. Any later artifact change requires rerunning the relevant check.

No new empirical model fitting was performed by these R-package checks. Original model producers are preserved for inspection and authorized reconstruction of their private dependency tree, and are not launched by the master script.

Metadata screening of the assembled public package found no SAV, RDS or RData files, symlinks or common credential patterns. The `id` field in the figure parameter manifest identifies model parameters, not participants. The `fold` field in the archived prediction-metric CSV identifies an aggregate fold summary, not an individual assignment. This screening is not a formal privacy guarantee; archived prediction outputs are distinct from the current grouped-table evidence.
