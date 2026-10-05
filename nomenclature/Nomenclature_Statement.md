# Nomenclature statement

Frozen analytical keys are retained separately from reporting descriptions. The archived HGNC map contains 18,067 exact approved-symbol matches, 237 unique alias/previous-symbol mappings and 17 unresolved or ambiguous entries across 18,321 genes. The 283 selected candidates have no unresolved symbols in that snapshot; C2orf16, FAM166C and ZBED9 map respectively to SPATA31H1, CIMIP2C and SCAND3. These display mappings do not change analytical membership. The 17 universe-wide symbol cases are separate from the six candidates lacking TSS annotation.

Reference-oriented descriptions were verified for all 11 Figure 3 inputs against GRCh37 (NC_000006.11 and NC_000019.9). The ten SNPs passed Mutalyzer normalization on 2 October 2026; the indel was resolved against its source panel and normalized on 3 October 2026. Reported_Variant_Nomenclature.tsv links the frozen keys to reporting descriptions; Nomenclature_Validation_Results.json contains the validation outputs and source evidence. Analytical keys, candidate membership and pair counts are unchanged.

| Frozen coordinate and unordered allele key (GRCh37) | Reference-oriented reporting description |
|---|---|
| 6:25685878:G:GTATAAT | NC_000006.11:g.25685895_25685900dup |
| 6:25710571:A:G | NC_000006.11:g.25710571G>A |
| 6:25715657:A:G | NC_000006.11:g.25715657G>A |
| 6:25878848:A:G | NC_000006.11:g.25878848A>G |
| 6:25918225:C:T | NC_000006.11:g.25918225T>C |
| 6:25918855:A:G | NC_000006.11:g.25918855G>A |
| 6:26093141:A:G | NC_000006.11:g.26093141G>A |
| 6:26123502:C:T | NC_000006.11:g.26123502T>C |
| 6:26485717:A:C | NC_000006.11:g.26485717A>C |
| 6:26599509:A:G | NC_000006.11:g.26599509A>G |
| 19:19379549:C:T | NC_000019.9:g.19379549C>T |

The PSC input chr6:25685878 G/GTATAAT is rs71544650. Although the summary-statistic header defines risk/other alleles, the [1000 Genomes Phase 1 integrated v3 record](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz) explicitly gives chr6:25685878, rs71544650, REF=G and ALT=GTATAAT, matching the source coordinate, identifier and alleles. The [PSC Genotype Imputation methods](https://discovery.ucl.ac.uk/id/eprint/10082444/1/668122.pdf) identify this panel as part of the combined reference. The source risk allele GTATAAT therefore corresponds to an insertion of TATAAT after GRCh37 position 25685878. Mutalyzer normalizes this insertion to NC_000006.11:g.25685895_25685900dup; the coordinate shift reflects HGVS representation within the repeat, not a change to the analytical key. The dbSNP 1000 Genomes observation corroborates the insertion. This source-panel match resolves the previously open deletion/duplication ambiguity for this reported input.

Reference and validation services: [Ensembl GRCh37 sequence API](https://grch37.rest.ensembl.org/documentation/info/sequence_region), [Mutalyzer](https://mutalyzer.nl/api/) and [HGVS reference and description checklist](https://hgvs-nomenclature.org/stable/recommendations/checklist/).

Additional input named in the text: rs2642438 corresponds to GRCh37 chr1:220970028, with the A/G pair used by the recorded DUSP10 support events. The GRCh37 Ensembl record contains A/C/G/T at this multiallelic site; the query and raw response are retained in Additional_Reported_Input_Identity.json. The rs identifier and frozen allele pair identify the reported input without inferring effect direction or assigning a new HGVS expression. Figure 3 still displays the original eleven verified inputs.
