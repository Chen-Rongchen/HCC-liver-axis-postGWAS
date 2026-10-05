# HCC liver-trait evidence

Code and derived data for studying recurrent liver-trait gene evidence and regional variant inputs in hepatocellular carcinoma.

## Run

```sh
pixi run check
pixi run replay-models
```

These commands validate the archived model inputs and reproduce existing unmasked model estimates.

## Files

- `code/`, `config/`: analysis code and settings.
- `derived/`, `supplementary_source_tsv/`, `figure_source_data/`: derived inputs and results.
- `metadata/`, `source_snapshot/`, `revision_source/`, `nomenclature/`: source records and validation files.

See [reproducibility notes](docs/REPRODUCIBILITY.md) for scope and source details. Raw GWAS data remain subject to their original access terms.
