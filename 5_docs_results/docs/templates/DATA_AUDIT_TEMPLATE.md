# DATA AUDIT — <source name>

| Field | Value |
|---|---|
| source_id | |
| Source (dataset/API name, exact files/endpoints) | |
| Publishing organization | |
| Retrieval date (UTC) | |
| Source version / release | |
| Source update date / cadence | |
| License / access conditions | |
| Unit of observation | |
| Sample size (actual, as ingested) | |
| Geography (resolution, vintage) | |
| Person-level? | yes/no |
| Geographic? | yes/no |
| Omics? | yes/no |
| Wearable? | yes/no |
| True participant linkage across modalities? | yes/no + evidence |

## Files / endpoints retrieved
Exact URLs, file names, byte sizes, sha256 (see MANIFEST.json), query parameters.

## Key variables
Variable name → meaning → units → how used downstream.

## Missingness
Measured, per key variable (counts and %), not guessed.

## Linkage strategy
Which identifiers join which tables; what is NOT joinable and why.

## Limitations and caveats
Selection, coverage, measurement, modeling (e.g. small-area estimates are modeled, not observed), vintage mismatches.

## Processed outputs
Tables written to data/processed and their row counts.

## Reproduce
`uv run python -m measure_it.<module>`
