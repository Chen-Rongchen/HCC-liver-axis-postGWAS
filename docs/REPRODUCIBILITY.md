# Reproducibility notes

The archive contains corrected mapped-input scores, unmasked model members and design matrices, existing model estimates, Rule A membership/results, source tables and nomenclature validation outputs. CURRENT_MANIFEST.tsv identifies the files by SHA-256. Source files cited in supplementary tables S4/S7 are retained in this snapshot; metadata/retained_source_files.tsv records their original hashes. The historical pixi manifest is separate from the current environment manifest.

`pixi run check` validates 16 design matrices. `pixi run replay-models` also reproduces existing unmasked coefficients and cluster-t intervals from the derived inputs. These commands do not rerun raw-GWAS processing, bootstrap, leave-one-block-out or Rule A fits. Rule A assignments and saved results are in S25/S26; the corresponding project-relative source scripts are under revision_source/. The recorded numerical environment is metadata/current_environment.json. Source snapshots have project-relative dependencies and are not a certified portable end-to-end raw-GWAS workflow.

Rule A was defined after the original findings and changes the analyzed gene population; it is not a causal decomposition. Source-mapped support does not establish an effector gene, and distinct minimum-P inputs do not establish independent signals. Historical sample-size/effect fields with unverified semantics are not certified RSS sample sizes or log-odds effects. Historical coloc/QTL results do not support reliable current signal-sharing inference.

The archive supplies code and shareable derived research data. Full raw GWAS datasets and individual-level reference resources are not redistributed; obtain them under the providers' access and reuse terms.

The four supplementary display workbooks are mapped to their original TSV schemas by metadata/submission_display_schema.json. Figure source data are in figure_source_data/; report-display identities and validation outputs are in nomenclature/.
