# Figure captions

One entry per published figure in `results/figures/`: the caption, then what data and code it comes from.
Drafts in `results/figures/drafts/` are not captioned here.

**Shared style.** Every figure uses `measure_it.figures.style`. It sets the palette: Okabe-Ito colourblind-safe
hues, one fixed colour per `data_layer`, validated for colour-vision deficiency (numbers in the module docstring).
It sets the fonts: DejaVu Sans and DejaVu Sans Mono, which ship with matplotlib, so no network fonts are used.
Sizes follow journal column widths (183 mm double column). Output is PNG at 300 dpi plus SVG with text converted to
paths. Rebuilding a figure from unchanged inputs gives byte-identical files.

Layer colours used throughout: person-level blue, condition-level molecular bluish green, population/geographic
orange, facility/provider/research reddish purple, ontology sky blue, measurement/technology vermillion,
derived dark grey. Each layer is also named in text, so colour is never the only cue.

Regenerate Figures 1, 2 and 10: `uv run python -m measure_it.figures.schematics` (about 3 s).
Tests: `uv run pytest tests/test_figures_schematics.py`. The tests measure every text element, so any text
outside the figure, overlapping other text, or outside its box fails. Inside a card, text must also stay at least
0.045 in from the side borders. Lines stacked inside a card are set at 1.2 × the font size, so an underscore never
touches the line below. They also check the product language and
confirm that each figure matches the registry or code it is drawn from.

---

## Figure 1. System architecture (`fig01_system_architecture.png`, `.svg`)

**Caption.** How the public-data alpha of the Measure It to Cure It engine turns a person-level phenotype into a
candidate deployment opportunity. The top row is the question path: person → phenotype → measurable biology →
technology → geography → clinic → deployment. The numbered badges show where each step is computed. Boxes are
grouped by the kind of join that connects them:

- **Person-level (solid arrows).** Layer 1 data (NHANES 2011-2014, the Stanford smartwatch cohorts, the NIH
  post-infectious ME/CFS omics deposits) are joined only on official participant IDs within one dataset.
  Datasets are never pooled. The MeasurementAdapter is the only modality-specific step. It emits per-participant
  phenotype scores: physiology against a reference, not diagnoses.
- **Concept-level (dashed).** Phenotype axes reach conditions through HPO terms and Mondo. Conditions reach two
  things through ontology IDs, never through people: condition-level molecular evidence (Layer 2), and candidate
  measurement classes and device categories drawn from ClinicalTrials.gov, NIH RePORTER and openFDA.
- **Ecological (dotted).** Burden and population context (Layer 3) and facilities, providers and research sites
  (Layer 4) are joined on place (FIPS, ZCTA, distance). They describe places and facilities, never individuals.

Every arrow inside a row uses that row's line style. The arrow from facilities to the deployment opportunity is
dotted because clinic capacity and research readiness are aggregated per geography.

A candidate deployment opportunity is a weighted sum. Its components (burden, vulnerability, diagnostic desert,
clinic capacity, research readiness) are kept as columns, and the ranks are re-run across weight sets. The
MCP / agent panel lists the 15 SPEC tools beside the layer each one reads. A tool with no data returns
UNKNOWN / NOT AVAILABLE. The bottom strip lists the provenance columns that every processed row carries. Object
ids let `trace_evidence` follow a recommendation back to its source rows.

**Derived from.** `measure_it.figures.schematics.figure1`. Some elements are generated at build time:
- per-layer source counts: `data_layer` in SOURCE_REGISTRY.yaml;
- adapter class names: `measure_it.measurements.adapters.list_adapters()`;
- MCP tool names: the AGENTIC MCP LAYER block of `docs/SPEC.md`;
- provenance column names: `measure_it.provenance.PROVENANCE_COLUMNS`.

Other elements are curated in the builder, from `docs/SPEC.md`, `docs/CONVENTIONS.md`, `docs/ADAPTERS.md` and
`docs/SCORING.md`: the box descriptions, the table names in the box footers, the join labels, and the mapping of
tools to layers (`TOOL_BANDS`).

Status at build time: `src/measure_it/mcp/` held only an empty `__init__.py`, so the agent panel is labelled
"tool contract from SPEC; server module in build". The label disappears on the next rebuild after a server
module lands (`mcp_server_present`).

---

## Figure 2. Public-data source architecture (`fig02_public_data_sources.png`, `.svg`)

**Caption.** Every public source in SOURCE_REGISTRY.yaml, grouped by data layer. Each row shows:
- the source's publisher and unit of observation;
- its status, as a glyph plus a word;
- the geographic levels it provides (point, ZCTA/tract, county, state, national);
- any access gating;
- arrows to the canonical processed tables its `processed_outputs` names.

Partitioned outputs (`<table>__<source>`) are folded into their table. Boxes in bold are integration tables
written by two or more sources. Each of these boxes is tall enough that its arrowheads stay separate, so the
number of sources writing a table can be counted (for example, 5 for `geo_context`).

Snapshot used (registry `generated_at` 2026-09-24T12:50:29+00:00, written by `registry_merge` in run 20260924T124442Z):
- **Sources:** 27, of which 26 ingested, 1 partial and 0 blocked.
- **Tables:** 59 canonical tables in 26 boxes, reached by 45 source-to-box arrows. 10 boxes integrate two or more
  sources.
- **NIH post-infectious ME/CFS deposits (partial):** the mapMECFS data files require an account. All 29 anonymous
  file probes returned HTTP 403 "Login Required" (`docs/MAP_MECFS_LINKAGE_AUDIT.md`). The per-participant
  actigraphy, HRV, orthostatic and CPET files were therefore not used. The open omics deposits enter the molecular
  layer as condition-level evidence only; this is the long arrow into `condition_molecular_evidence`.
- **Census ACS and Data.gov:** the keyed API endpoints were not used. The keyless ACS Summary File and the
  keyless Data.gov catalogue origin were used instead.
- **HGNC and Reactome:** reference data used by the analyses. They write no processed table.

Pending registry merge: the `data/raw/*/registry_entry.yaml` fragments are ahead of SOURCE_REGISTRY.yaml, and the
ClinicalTrials.gov fragment now writes `ctgov_facility_summary` instead of
`research_site_registry__clinicaltrials_gov`. Merging changes that table name and `generated_at`. The counts above
do not change (checked on an in-memory merge).

The geography dots summarise the registry's free-text `geographic_resolution`. For example, ZCTA/tract includes
ZIP codes placed at their ZCTA centroid.

**Derived from.** `measure_it.figures.schematics.figure2(registry)`, which reads SOURCE_REGISTRY.yaml and nothing
else. Field by field:
- `data_layer` sets the band. `status` sets the glyph and the counts.
- `publisher` is shown with parentheticals removed.
- `unit_of_observation` is shown as its first clause.
- `geographic_resolution` is mapped to levels by the keyword rules in `resolution_levels`.
- `access_conditions` gives the gating pills (`access_flags`).
- `processed_outputs` is converted by `canonical_table`. An entry that is not a table name produces "no processed
  table written".

Arrows are grouped by the exact set of sources that write each table (`table_chips`). Source names are the only
hand-written text: abbreviations in `DISPLAY_NAMES`, with the registry `name` used for any source not listed
there. After any registry merge, rebuild the figure to update the snapshot. The tests check that every source,
every canonical table and every status count appears in the figure.

---

<!-- data-figures:start (generated by measure_it.figures.data_figures; do not edit by hand) -->

**Data figures (3-9).** Built by `measure_it.figures.data_figures` (modules `fig_wearable`, `fig_molecular`, `fig_geography`, `fig_drilldown`; maps `data_maps`, interactive pages `fig_maps_html`). Every number on a figure and in these captions is read from the processed or result tables, or returned by a measure_it query function, when the figure is built; none is typed in. Regenerate: `uv run python -m measure_it.figures.data_figures` (about 20 s; also written by the pipeline step `figures_data`). Tests: `uv run pytest tests/test_figures_data.py`. The build refuses to save a figure whose text leaves the canvas, overlaps other text or uses a banned product-language phrase. Maps use equal-area projections (contiguous U.S. EPSG:5070; Alaska and Hawaii insets in their own projections, not to scale). Interactive versions of Figures 6-9 are in `results/maps/` (Plotly; one shared `plotly.min.js`, so they open offline).

---

## Figure 3. Wearable phenotype differences (`fig03_wearable_phenotype_differences.png`, `.svg`)

**Caption.** Objective activity and heart-rate phenotypes in public person-level wearable data. **(a)** NHANES 2011-2014 wrist accelerometry (26 features): difference between cases and controls in SD units, adjusted for age, age² and sex, with 95% CI; filled markers pass BH-FDR 5% across the 26 features. The ME/CFS-like group is a constructed PROXY (unexplained fatigue with functional limitation; 132 cases / 7,664 controls; 14 of 26 features pass FDR), not ME/CFS. Fibromyalgia is defined by a prescription reason-for-use code (Rx-defined; 2013-14 cycle only; 34 cases, below the 50-case power threshold; 13 of 26 pass FDR), so it may reflect the medication as much as the condition. Functional limitation (2,340 cases; 24 of 26) and all-cause mortality, the positive control (770 deaths by 2019-12-31; 21 of 26) have the same sign for 23 of 26 features. These are associations in one cross-sectional sample, not diagnostic performance; NHANES 2011-2014 has no heart rate. **(b)** Stanford acute COVID-19 cohorts (Mishra 2020 + Alavi 2022), labelled acute COVID-19, NOT Long COVID: mean resting-HR and step z-scores against each person's pre-infection baseline around the true onset and around the same 84 placebo-eligible participants' placebo onsets (83 with step data; placebo windows end at day +14 because they must precede the infection), with participant-bootstrap 95% CIs. Resting HR peaks at 1.03 SD on day 0, but the pre-specified days 0-14 contrast with placebo is null (+0.11 SD, 95% CI -0.20 to +0.45; p = 0.53); steps fall (-1.10 SD, -1.33 to -0.88; p < 0.001), which is expected behaviour when ill rather than a latent physiological signal. **(c)** Uwakwe 2025, the only release with a Long COVID label (31 self-reported Long COVID vs 95 recovered): age/sex-adjusted differences of seven whole-record resting-HR features; 0 of 7 have p < 0.05 (unadjusted: 1 of 7, 0 after BH-FDR). The wearable-only model has AUROC 0.50 (95% CI 0.40 to 0.61; label-permutation p = 0.30): a null result.

**Derived from.** `results/tables/nhanes_phenotype_signatures.csv` and `nhanes_signature_summary.csv` (`measure_it.wearables.signatures`); `stanford_acute_trajectory_by_day.csv` (group `pooled_placebo_eligible`), `stanford_acute_endpoints.csv`, `stanford_phenotype_signatures.csv` (rows `*__adj_age_sex` and the unadjusted Hedges' g rows) and `stanford_uwakwe_model_performance.csv` (`measure_it.wearables.stanford_models`). The feature grouping and short labels are display choices (`fig_wearable.FEATURE_GROUPS`); the test checks that the 26 grouped features are exactly the table's.

---

## Figure 4. Clinical-only versus clinical + wearable model performance (`fig04_clinical_vs_wearable_models.png`, `.svg`)

**Caption.** Does wearable physiology add to demographic and clinical variables? **(a)** Cross-validated AUROC with participant-bootstrap 95% CIs for five feature sets per NHANES target (penalised logistic regression, 5 x 5 CV) and for the three sets available in the Uwakwe Long COVID release (L2 logistic, 10 x 5 CV; its baseline is age, sex, BMI, race and vaccination, and there is no clinical set). Wearable features alone have a lower point estimate than demographics for 6 of 6 NHANES targets (for example functional limitation 0.731 vs 0.755). **(b)** Paired increments: Test 1 (demographics + wearable minus demographics) and Test 2 (clinical + wearable minus clinical), with the primary model (filled circle), the other two model families (grey) and the 2011-12 -> 2013-14 temporal split (open diamond). The reading column applies the rules of `results/WEARABLE_NHANES_RESULTS.md`: supported when the primary CI excludes 0 and the matching negative control has p <= 0.05; for Test 2, 'robust' also needs CI > 0 in all three families, in the temporal split and a Holm-adjusted bootstrap p <= 0.05 across the six targets (post hoc). Test 1 is supported for 5 of 6 targets (functional limitation, fatigue, fair/poor health, depressive symptoms, mortality); Test 2 is supported and robust only for functional limitation, meets the rule but is not robust for fatigue, depressive symptoms, and is null for fair/poor health, ME/CFS-like proxy, mortality. In Uwakwe Long COVID, adding the wearable block changes AUROC by -0.024 (95% CI -0.076 to +0.027; wearable-block permutation p = 0.47): null. **(c)** Label-permutation negative controls. NHANES stores only the summary of the 200 permutations of each wearable-only model (null mean, 95th percentile, maximum); 6 of 6 observed AUROCs lie above their null maximum (smallest margin: ME/CFS-like PROXY 0.587 vs 0.581), p = 1/201 each. The Uwakwe panels show all 1,000 stored draws: the wearable-only model sits inside its null (p = 0.30), the demographic/clinical baseline above it (p = 0.001).

**Derived from.** `results/tables/nhanes_model_performance.csv`, `nhanes_increments.csv` (Holm p from `measure_it.wearables.nhanes_report.holm_boot`), `nhanes_permutation_null.csv`, `nhanes_shuffle_control.csv`, `nhanes_targets.csv`; `stanford_uwakwe_model_performance.csv`, `stanford_uwakwe_model_increments.csv`, `stanford_uwakwe_permutation_nulls.csv`, `stanford_uwakwe_permutation_null_draws.csv`. The readings are recomputed by `fig_wearable.nhanes_verdicts` with the write-up's rules; the test checks them against `results/WEARABLE_NHANES_RESULTS.md`.

---

## Figure 5. Condition -> molecular evidence -> physiological system -> measurement class (`fig05_molecular_measurable_physiology.png`, `.svg`)

**Caption.** Condition-level molecular enrichment from public databases (Open Targets, GWAS Catalog, mapMECFS published analytes, Reactome), not participant-linked and not patient multi-omics. **(a)** The demo cluster (Long COVID, ME/CFS, POTS): each condition's top enriched Reactome pathways per evidence source (15 pathways, 18 data-derived condition -> pathway edges, colour = source), each pathway's physiological system (Reactome hierarchy plus the curated anchor table), and the curated system -> measurement-class links (dashed; 15 classes; a curated link is an assumption, not evidence that a measurement detects the condition). 9 of the 18 condition -> pathway edges are post-hoc flagged (dotted): most of the pathway's genes are drug-target-only genes of trialled drugs, where one drug mechanism can list a whole gene family. Bold measurement classes belong to the wearable autonomic/activity bundle. The three cells beside each system give its support tier for Long COVID, ME/CFS and POTS. **(b)** Support tier of every condition -> system link: >= 2 source families, one non-literature family, literature co-mention only, tested but not supported, or no / too few genes; † marks a post-hoc flag and the arrow gives the tier without the flagged sources (16 flagged cells of 84). Agreement between sources is not independent replication: they share curated knowledge and literature.

**Derived from.** `results/tables/fig5_nodes.csv` (238 nodes) and `fig5_edges.csv` (475 edges), written by `measure_it.omics.graph`; `data/processed/measurable_biology.parquet` (support tiers, post-hoc flags, tier without flagged sources). Only the top 3 pathways per source and condition are drawn; the tables hold every edge. The wearable bundle membership comes from `configs/relevance.yaml`.

---

## Figure 6. U.S. burden and unmet need for Long COVID (`fig06_long_covid_burden_unmet_need.png`, `.svg`)

**Caption.** **(a)** Adults currently experiencing long COVID, BRFSS 2023 small-area model aggregated to states (modelled; states absent from the 2023 file, Kentucky and Pennsylvania, are model predictions), by state: burden evidence level A (a direct measure) at STATE resolution, all 51 jurisdictions (50 states + DC). **(b)** The same estimates ranked with their model 95% intervals (from 4.7% in DC to 9.1% in WV); most intervals overlap, so the state ranking is uncertain. **(c)** County diagnostic desert for long COVID: the mean percentile of burden, SVI, low relevant-provider density and few trials within 50 km (3,144 counties). Printed on the figure: no observed county long-COVID measure exists; the burden component is a BRFSS 2023 multilevel-regression and post-stratification estimate whose within-state differences come only from modelled county composition, urbanicity and covariates (results/BRFSS_SAE_RESULTS.md); it is not observed county prevalence. The evidence-level key names the levels A-D and marks the one used.

**Derived from.** `measure_it.geography.query.get_condition_burden('long_covid', 'state')` and `('long_covid', 'county')` (inherited flag); `data/processed/geo_condition_features.parquet` (`diagnostic_desert` and its four percentile components, `desert_burden_resolution` = 'county (modelled small-area estimate)'). Source: Measure-It MRP small-area estimate from CDC BRFSS 2023 + ACS 2020-2024 (BRFSS 2023 LLCP public file + ACS 2020-2024 5-year; model M1_plus_acs_covariates; run 2026-10-07T18:54:56+00:00). Interactive version: `results/maps/fig06_long_covid_burden_unmet_need.html`.

---

## Figure 7. Provider and research infrastructure overlay (`fig07_infrastructure_overlay.png`, `.svg`)

**Caption.** **(a)** Relevant-provider density by county: individual NPIs per 100,000 residents in the long-COVID specialty groups (primary care, pulmonology, cardiology, neurology, PM&R, infectious disease; curated in configs/relevance.yaml), quintile classes. **(b)** The specialist-only variant: the same groups without primary care (zero class plus quartiles of the non-zero counties). Specialists are 13.7% of the relevant NPIs, and 1,354 of 3,144 counties have none, so panel a mostly maps primary-care supply. A self-reported taxonomy does not mean a provider evaluates or treats the condition. **(c)** 17,818 active HRSA health-center service-delivery sites (17,744 geocoded and drawn); 263 Long COVID or ME/CFS trials with a U.S. site (condition named literally; 724 site rows), drawn at 194 CITY CENTROIDS because ClinicalTrials.gov publishes a city geoPoint, not an address; and 146 NIH-funded organisations (131 geocoded: 115 at the awardee ZIP centroid, 16 at a source point) with 974 projects whose title or abstract names either condition. A trial or grant is a research-activity signal, not evidence that a measurement works.

**Derived from.** `geo_condition_features` (`relevant_providers_per_100k`, `relevant_specialists_per_100k`, condition long_covid); `facilities__hrsa` (`is_service_delivery_site`, `is_active`); `measure_it.ingestion.clinicaltrials.trials_for_condition(..., literal_only=True, us_sites_only=True)` for long_covid and me_cfs, then `trial_sites`; `measure_it.facilities.query.find_relevant_research_centers` for both conditions, NIH project ids merged per facility. Interactive version: `results/maps/fig07_infrastructure_overlay.html`.

---

## Figure 8. Deployment opportunity map (`fig08_deployment_opportunities.png`, `.svg`)

**Caption.** Candidate deployment opportunities for the primary query, Long COVID or ME/CFS x wearable autonomic/activity monitoring, by county. **(a)** Equal-weight composite (mean percentile of burden, vulnerability, diagnostic desert, clinic capacity and research readiness) for the 2,395 ranked counties (counties of the 50 states + DC with population ≥ 10,000). Grey: 725 counties scored but not ranked (population < 10,000). Hatched: 24 counties with an incomplete set burden (no CMS value for ME/CFS; 9 of them are the Connecticut planning regions, which CMS publishes only on the legacy counties), not ranked. Numbered: the top 10. The long-COVID burden is a modelled BRFSS 2023 small-area estimate (not observed county prevalence) and the ME/CFS burden a claims proxy (level C). **(b)** The top 10 with their 5th-95th percentile ranks over 1,000 Monte Carlo draws (burden error, inherited-burden deviation, Dirichlet weights): interval widths 38-239 ranks, and 4 of the 10 are in the top 10 in at least half of the draws. **(c)** Rank of the same counties under the 9 named weight sets (incl. `evidence_weighted`, which adds the measurement's own performance evidence and expected yield) and their share of top-25 placements over 1,000 random weight vectors. Kendall's W across the named sets is 0.829 (0.855 without burden only; 0.708 across random weights): the overall ordering is concordant, but which counties make the short list depends on the weights. A ranked county is a candidate for pilot evaluation, not a validated diagnostic pathway.

**Derived from.** `measure_it.scoring.recommend.rank_deployment_opportunities('Long COVID or ME/CFS', 'wearable autonomic monitoring', 'county', top_n=10)` (top 10, rank intervals, P(top 10), ranks under each weight set); `data/processed/deployment_opportunities.parquet` (composite, `rank_eligible`, `burden_incomplete`); `results/tables/test7_region_stability_primary.csv` and `test7_summary.csv`. Interactive version: `results/maps/fig08_deployment_opportunities.html`.

---

## Figure 9. Example geography drill-down: San Diego County, California (`fig09_san_diego_drilldown.png`, `.svg`)

**Caption.** **(a)** San Diego County (population 3,288,774) with the 50 km candidate-pool radius around its Census internal point, HRSA health-center sites, Long COVID / ME/CFS trial sites (city centroids) and NIH-funded organisations. **(b)** The 18 candidate implementation or study partners returned by `find_candidate_clinics` for Long COVID and ME/CFS x wearable autonomic/activity monitoring, merged by facility, each with the characteristic ranking that selected it, the query that returned it (Long COVID, ME/CFS or both) and its counts (matched specialty/implementer groups, trials, NIH core projects, technology-experience trials, community-access flags; maxima over the two queries). The six characteristics are ranked separately and never combined; each facility has characteristics suggesting it may be a viable implementation or study partner, which is not evidence of capability for this measurement. For ME/CFS no facility in the pool has NIH research activity (absence in these public registries, not proof that none exists). **(c)** Burden by condition with evidence level and source resolution: Long COVID is a modelled BRFSS 2023 small-area estimate (not observed county prevalence); 8 conditions have no usable burden estimate (level D, UNKNOWN / NOT AVAILABLE). **(d)** SVI (published percentile) and the county's five opportunity components for the primary query: burden 0.26, vulnerability 0.71, diagnostic desert 0.09, clinic capacity 0.83, research readiness 0.98; composite 0.57, rank 971 of 2,395 (Monte Carlo 5-95%: 243-1,620): low burden and desert percentiles keep a county with high research readiness out of the short list. **(e)** Individual NPIs per 100k in the relevant non-primary-care specialty groups, with the U.S. rate (sum of county counts over the 50 states + DC divided by their population) as a tick. **(f)** Trials with a site in the county by condition (condition named literally), NIH RePORTER projects whose title or abstract names any registry condition (60, 27 core, at 12 organisations; Long COVID 20, ME/CFS 0) and HRSA sites (258).

**Derived from.** `measure_it.geography.query.get_geographic_context('San Diego County, California')` and `get_condition_burden(condition, 'county', geo_id='06073')`; `measure_it.scoring.recommend.candidate_sites('06073', ['long_covid', 'me_cfs'], 'wearable_autonomic_activity_monitoring')` (the demo recommendation's call to `facilities.matching.find_candidate_clinics`, default radius from configs/scoring.yaml) with facility coordinates from `facilities`; the county's row of `deployment_opportunities`; `provider_density_county` for the U.S. reference rates. Map layers as in Figure 7. Interactive version: `results/maps/fig09_san_diego_drilldown.html`.

<!-- data-figures:end -->

---

## Figure 10. Modularity: wearables today, MeasurementAdapter, CAPRIO capillaroscopy tomorrow (`fig10_adapter_modularity.png`, `.svg`)

**Caption.** One interface separates the measurement modality from everything downstream.

**Top: the two implemented adapters,** as their own `describe_measurement()` reports them. Both use public data.
- **WearableAdapter:** wrist accelerometry; trained on NHANES 2011-2014; 38 embedding features; four phenotype
  score columns scored against an age/sex reference.
- **WearableHeartRateAdapter:** consumer-wearable resting heart rate, plus steps where released; trained on the
  Stanford cohorts; 22 features; two scores against the person's own baseline.

Each card also shows the adapter's reproduction check. `embed()` on raw NHANES minutes reproduces all 38 stored
features for 523 participants (maximum absolute difference 2.9e-11). `embed()` on the released raw minute files
reproduces the 15 stored longitudinal features for 126 participants (maximum absolute difference 0).

**Middle: the interface.** The abstract `MeasurementAdapter` has four methods: `preprocess`, `embed`,
`phenotype_score` and `describe_measurement`. The panel gives their signatures and docstrings, plus the 11
`MeasurementDescription` fields that are all the downstream layers read about a modality.

**Bottom: the `CapillaroscopyAdapter` stub.** A single dashed "implements (stub)" arrow joins it to the
interface. The card beside it probes the same class and is not a second adapter. No CAPRIO or other private data
are used or fabricated, and `requires_private_data` is True. Its `preprocess`, `embed` and `phenotype_score` raise NotImplementedError, and
each error message names what CAPRIO must supply. `describe_measurement()` already works, so implementer
groups, deployment complexity, FDA context and linked conditions resolve before any data exist.

**Right: the downstream path.** The adapter-agnostic contract in `measure_it.measurements.adapters` receives two
inputs:
- phenotype scores go to `scores_to_long` and `summarize_scores`;
- the description goes to `validate_description` through `measurement_context`.

The contract feeds the downstream layers, which do not change. Swapping in CAPRIO changes the measurement ids
(capillaroscopy, microvascular_function), not the code. With today's public data the swap barely changes the
candidate deployment opportunities: in sensitivity analysis S5 the `nailfold_capillaroscopy` bundle, scored through
the identical code path, gives a county ranking with composite Spearman 0.994 and top-10 overlap 9 of 10 against the
wearable bundle (`results/tables/test7_sensitivity_comparisons.csv`), because burden, vulnerability and desert do not
depend on the measurement. The panel's closing sentence reads these two numbers from that table.

**Derived from.** `measure_it.figures.schematics.figure10()`, drawn from the code at build time:
- **Adapters:** each registered adapter is instantiated through `measure_it.measurements.adapters`
  (`list_adapters`, `get_adapter`), and its card shows its `describe_measurement()` output.
- **Interface:** the methods come from `MeasurementAdapter.__abstractmethods__`, with `inspect` signatures and the
  first docstring line of each. The fields come from `dataclasses.fields(MeasurementDescription)`.
- **Stub behaviour** is probed, not assumed: `_probe_stub_methods` calls each data method with None and records
  the NotImplementedError. The "CAPRIO supplies ..." lines come from `capillaroscopy_adapter.CAPRIO_REQUIREMENTS`.
- **Validation lines** come from `results/tables/wearable_adapter_validation_summary.json` and
  `wearable_hr_adapter_validation_summary.json`. A line is omitted if its file is missing or does not report
  `all_pass`.
- **Curated text:** the list of downstream layers and the explanatory sentences, from `docs/ADAPTERS.md`
  sections 3 and 6.
