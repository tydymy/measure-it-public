# SPEC — "Measure It to Cure It" public-data alpha (verbatim project brief, 2026-09-23)

# OBJECTIVE

Build a fully reproducible, public-data-only alpha of the **"Measure It to Cure It" AI translation and deployment engine** for the TOPx HHS Invisible Illness challenge.

The system should demonstrate that we can:

**person-level multimodal phenotype → objective wearable phenotype → disease/phenotype ontology → molecular/multi-omics context → candidate measurable physiology/technology → geographic burden/unmet need → clinic/research-site readiness → recommended deployment opportunity**

The current private CAPRIO/MAESTRO capillaroscopy datasets MUST NOT be required for this demonstration.

Instead, use publicly available wearable and physiologic datasets as the measurement modality.

Capillaroscopy should appear only as an example of a future interchangeable measurement adapter. The architecture must make it straightforward to replace the wearable module later with CAPRIO nailfold-capillaroscopy embeddings without rewriting the geographic, clinic, ontology, molecular, or agentic layers.

The scientific thesis is:

> Invisible illnesses are heterogeneous and poorly represented by diagnosis codes alone. Objective physiological measurements can reveal latent phenotypes. Public biomedical data can identify which measurements may be informative, while U.S. federal open data can identify the regions, facilities, and clinicians where those measurements could have the greatest deployment value.

The product is therefore NOT merely: a disease prevalence map, a wearable classifier, a provider directory, or an LLM medical chatbot.

It is a **modular measurement-validation and deployment engine**.

---

# NON-NEGOTIABLE SCIENTIFIC GUARDRAILS

Do not optimize for an impressive-looking result at the expense of scientific validity.

Do not fabricate joins between unrelated datasets.

There are four fundamentally different data layers:

1. **Person-level biomedical data** — wearable measurements; symptoms/questionnaires; clinical measurements; labs; diagnoses where available; omics only where truly linked to the same participant
2. **Condition-level molecular evidence** — genes; variants; pathways; transcriptomics; proteomics; metabolomics; published biomarkers
3. **Population/geographic data** — disease or symptom burden; social vulnerability; demographic context; health-care utilization
4. **Facility/provider/research infrastructure** — physicians; specialties; health centers; clinical-trial sites; NIH-funded institutions; deployable measurement technologies

NEVER imply that a participant from layer 1 is the same individual as anyone represented in layers 3 or 4.

Geographic joins are ecological/population-level joins.

Do not infer individual patient location from public clinical datasets.

Do not infer that an omics signature belongs to a wearable participant unless the source explicitly contains shared participant identifiers linking those modalities.

If omics and wearables are not person-linked, call the result: **condition-level molecular enrichment** and NOT: **patient multi-omics**

Every derived table must contain provenance fields.

Every score must expose its component variables.

Do not use an LLM-generated number as data.

Do not silently substitute proxy measures for disease prevalence.

If a condition lacks direct prevalence data, label the variable as a proxy.

---

# PRIMARY DEMONSTRATION QUESTION

Build the system so it can answer queries such as:

> "Where in the United States would deploying wearable-based autonomic or activity monitoring be most valuable for detecting objective abnormalities associated with Long COVID, ME/CFS, dysautonomia/POTS-like phenotypes, or related invisible illnesses, and which clinics or research sites are best positioned to implement the technology?"

And ultimately:

> "Given phenotype X and candidate measurement technology Y, where should we deploy Y to generate the most useful clinical evidence?"

The system should be general enough that `wearable physiology` can later be replaced with `nailfold capillaroscopy` without changing the downstream architecture.

---

# PUBLIC DATA SOURCES

## A. Person-level wearable + phenotype layer

### 1. NHANES 2011–2014

Use as the large population-scale public clinical/wearable dataset.

Relevant components should include, where available: demographic variables; medical-condition questionnaires; cardiovascular questionnaires; diabetes; sleep; physical functioning; prescription medications; blood pressure; BMI/body composition; laboratory measurements; seven-day wrist accelerometry; minute/hour/day accelerometry summaries.

Prefer the released minute/hour/day accelerometry tables for the primary pipeline. Raw 80-Hz acceleration can be supported as an optional advanced mode.

Join NHANES tables ONLY using official participant identifiers.

Construct one row per participant with derived wearable features. Examples: mean daily activity; total activity; sedentary fraction; activity fragmentation; day-to-day variability; circadian regularity; relative amplitude; interdaily stability; intradaily variability; active-period duration; low-activity burden; sleep-window proxies if scientifically supportable from available variables; weekday/weekend differences; valid wear time.

Do NOT invent heart-rate or HRV data for NHANES if it is not present.

### 2. Stanford Snyder Lab Long COVID wearable dataset

Use this as an invisible-illness-specific wearable dataset. Read the actual data dictionary before defining features. Do not assume particular columns from the publication alone.

Where available, derive: heart-rate baselines; resting heart rate; deviations from individual baseline; activity/step features; recovery trajectory after infection; longitudinal variability; wearable-derived Long COVID discrimination features.

Use the associated published analysis as a reproduction benchmark where practical.

### 3. NIH post-infectious ME/CFS deep-phenotyping dataset / MapME/CFS

Attempt to use this because it provides unusually rich invisible-illness phenotyping. Potential modalities reported by the study include: actigraphy; heart-rate variability; orthostatic challenge; CPET; strength; metabolic testing; immune measurements; clinical labs; metabolomics; proteomics; lipidomics; mitochondrial sequencing; transcriptomics; stool/metagenomic measurements. Several molecular datasets are publicly deposited in GEO/SRA.

IMPORTANT: First determine whether participant IDs actually allow modalities to be joined. Create `docs/MAP_MECFS_LINKAGE_AUDIT.md`. Document: downloaded datasets; participant identifiers; number of subjects per modality; whether IDs overlap; exact overlap counts; whether wearable/actigraphy and omics can truthfully be linked.

If yes, create a true linked multimodal subset. If no, keep the molecular data as condition-level evidence. The pipeline MUST work even if MapME/CFS access cannot be automated.

### 4. Optional public wearable benchmarks

Use only if useful: PhysioNet wearable datasets; CovIdentify; public Apple Watch/accelerometry sleep datasets. These may validate signal-processing components. Do not treat healthy sleep cohorts as invisible-illness cohorts.

# B. DISEASE AND PHENOTYPE ONTOLOGY

Create a reusable disease ontology layer. Start with: Long COVID / post-COVID condition; ME/CFS; POTS; dysautonomia; fibromyalgia; Ehlers-Danlos syndrome / hypermobility spectrum disorders; mast-cell-related syndromes where appropriate; Lyme disease / post-infectious syndromes; migraine; IBS; gastroparesis; endometriosis.

Build `condition_registry.parquet`. Fields should include: canonical_condition_id; preferred_name; aliases; ICD-10-CM codes where appropriate; MONDO IDs; EFO IDs where available; MeSH IDs where available; HPO phenotype terms; parent_condition; infectious/post-infectious flag; autoimmune flag; autonomic flag; vascular flag; pain flag; fatigue/PEM flag; GI flag; neurologic flag; ontology provenance.

Use public ontology resources preferentially: Mondo; EFO; HPO; Monarch mappings; EMBL-EBI OLS. Mondo explicitly includes cross-references to ICD-10-CM and numerous biomedical disease ontologies. Do not require proprietary terminology services for the MVP.

# C. MOLECULAR / MULTI-OMICS EVIDENCE GRAPH

Construct `condition_molecular_evidence.parquet` and/or a small graph database.

Sources:
- **Open Targets** — official Open Targets MCP and/or GraphQL API. Retrieve where available: target-disease associations; genetics; variants; GWAS evidence; pathways; known drugs; biological evidence.
- **GWAS Catalog** — current supported REST API. Retrieve: traits; studies; variants; mapped genes; EFO terms; ancestry metadata.
- **GEO** — search for relevant disease transcriptomic/proteomic studies. Do not attempt to download every study. Create a curated query mechanism for target conditions.
- **SRA** — use when needed for public sequencing datasets associated with target diseases.
- **MapME/CFS molecular datasets** — use as a high-value disease-specific source when publicly accessible.

For every molecular edge store: condition; ontology ID; evidence type; gene/variant/pathway/analyte; source database; source accession; evidence direction if defined; study population; sample size if easily available; date retrieved.

DO NOT create an arbitrary "omics score" from unrelated evidence. The immediate purpose of this layer is to connect **clinical phenotype → possible measurable biology** rather than to claim mechanistic proof.

# D. MEASUREMENT / TECHNOLOGY DISCOVERY LAYER

This is essential. The engine should not assume wearables are the only solution. Construct `measurement_registry.parquet`. The system should discover objective measurement approaches being used for invisible illnesses. Potential measurement classes include: wearable HR; ECG; PPG; HRV; accelerometry; posture detection; continuous oxygen saturation; respiratory rate; sleep measurements; continuous temperature; CPET; autonomic testing; endothelial-function testing; vascular imaging; capillaroscopy; retinal imaging; digital gait analysis; blood biomarkers; metabolomics; proteomics; immune assays; microvascular measurements.

## ClinicalTrials.gov
Search target conditions. Extract: intervention names; device names; outcome measures; physiological endpoints; eligibility criteria; trial phase/type; recruiting status; sponsors; facility names; facility cities/states; measurement technologies. Use this to infer what objective measurements researchers are actively attempting to deploy. Do NOT interpret use in a clinical trial as proof of effectiveness.

## NIH RePORTER
Search grants using target conditions plus terms such as: wearable; biomarker; diagnostic; sensor; imaging; autonomic; microvascular; capillary; digital health; point of care. Extract: project title; abstract; PI; institution; city/state/ZIP; funding; project dates; disease; measurement technology. Use this as a research-capability signal.

## openFDA medical device APIs
Search: device classification; 510(k); PMA where relevant; device registration/listing. For candidate technologies, determine whether a corresponding regulated device category exists and whether specific devices have FDA clearance/approval records. Store: device name; product code; device class; specialty; clearance identifier if applicable; applicant/manufacturer; decision date; intended device category. A regulatory record is a **deployment-readiness signal**, not evidence that a technology diagnoses the target illness.

# E. PUBLIC GEOGRAPHIC BURDEN LAYER

Canonical geography should preferably use: state FIPS; county FIPS; ZCTA where available; latitude/longitude for facilities. Do not force all data to a finer resolution than the original source supports. Create `geo_condition_features.parquet`.

- **CDC PLACES** — county, place, census tract, and ZCTA health measures. Include relevant measures such as: self-rated health; disability-related measures; depression where applicable; physical inactivity; chronic disease; obesity; sleep; access to care. These are population context variables and proxies where appropriate.
- **CDC Long COVID / Household Pulse Survey** — public Long COVID prevalence estimates. Respect the geographic resolution provided by the dataset. Do not downscale state estimates to counties and present them as observed county prevalence.
- **CMS Mapping Medicare Disparities** — county/state chronic-condition prevalence, utilization, hospitalization, and related measures where applicable. For each condition define `burden_evidence_level`: A = direct condition measure; B = closely matching coded condition; C = symptom/comorbidity proxy; D = no usable burden estimate. Always expose this label.
- **CDC disease-specific surveillance** — use where applicable, for example Lyme/tickborne disease.
- **CDC/ATSDR Social Vulnerability Index** — join at the appropriate FIPS geography. Store: overall SVI; socioeconomic theme; household characteristics; racial/ethnic minority status where present in the published SVI construct; housing/transportation theme. Do not treat vulnerability as disease prevalence.
- **Census ACS** — population denominator; age structure; disability; insurance; poverty; income; transportation; broadband/internet access where relevant; rural/urban context if appropriately derived.

# F. PROVIDER AND FACILITY LAYER

Create `clinic_registry.parquet`.

**NPPES** — use current downloadable NPPES V2 data. Include: NPI; organization/person type; provider name; taxonomy; specialty; practice address; city; state; ZIP; latitude/longitude after geocoding where legally/permissibly possible. Create specialty groupings relevant to invisible illness: cardiology; electrophysiology; neurology; rheumatology; immunology/allergy; infectious disease; pulmonology; gastroenterology; pain medicine; physical medicine/rehabilitation; primary care; gynecology; vascular medicine/surgery. Do not assume an NPI specialty means that provider treats Long COVID, ME/CFS, or POTS.

**HRSA** — add health-center and medically underserved infrastructure. Use the public service-delivery-site dataset. Store: facility; location; facility type; available metadata; geographic identifiers. This can identify community sites where deployment could improve access.

# G. CLINICAL RESEARCH READINESS

Create `research_site_registry.parquet`. Combine:
- **ClinicalTrials.gov site history** — for each facility/location: number of relevant trials; number recruiting; conditions studied; technologies studied; study recency; device/digital-health experience.
- **NIH RePORTER** — for each organization/geography: number of relevant grants; active projects; cumulative funding if appropriate; investigators; research topics.

Research readiness must remain separate from patient burden.

# H. DATA.GOV / U.S. OPEN DATA DISCOVERY AGENT

Build a reusable Data.gov discovery component. The current Data.gov catalog API is a metadata catalog. Use it to: search for relevant U.S. government datasets; inspect metadata; identify publisher; identify download/API resources; route retrieval to the original agency.

Implement `search_us_open_data(query)`. Return: title; agency; description; access level; resource URLs; format; modified date. Do NOT assume Data.gov itself hosts the underlying records.

---

# CORE DATA MODEL

Use DuckDB + Parquet as the canonical analytical store. Recommended tables:

```text
participants
participant_conditions
participant_labs
participant_wearable_features
participant_phenotype_embeddings

conditions
condition_ontology_mappings
condition_molecular_evidence

measurements
condition_measurement_evidence
measurement_regulatory_status

geographies
geo_condition_burden
geo_context
geo_vulnerability

providers
facilities
clinical_trials
trial_sites
nih_projects

deployment_candidates
```

Every table should include where applicable: source_name; source_record_id; source_version; retrieved_at; source_geographic_resolution; evidence_type; evidence_level; provenance_notes.

---

# PHASE 1 — BUILD THE WEARABLE PHENOTYPE

Start simple. Do not immediately train a deep neural network. Construct interpretable wearable features from NHANES and the Long COVID dataset.

Create baselines: 1. logistic regression; 2. regularized logistic regression; 3. random forest or gradient boosting. Only add a temporal/deep model if it improves appropriately held-out performance.

Targets can include: disease labels where available; symptom burden; functional limitation; Long COVID status; physiologic cluster membership.

For unsupervised analyses, compare: PCA; UMAP only for visualization; Gaussian mixture / HDBSCAN or another defensible clustering method. Do not describe a cluster as a biological subtype without independent evidence.

Use participant-level train/test splitting. Avoid temporal leakage.

Report: AUROC; AUPRC; sensitivity; specificity; calibration; confidence intervals where practical. Perform permutation-label negative controls.

# PHASE 2 — REPRESENT THE "DIGITAL PERSON"

Create a common patient representation: Demographics + Symptoms/conditions + Labs + Medications if available + Wearable physiology + Functional status → Digital Phenotype Vector.

If a public dataset contains truly linked omics: Digital Phenotype Vector + participant-linked omics → Multimodal Person Representation. Otherwise: Digital Phenotype Vector → Condition ontology → Condition-level molecular evidence. Label the second approach clearly as **molecular enrichment**, not individual multi-omics.

# PHASE 3 — IDENTIFY WHAT COULD MAKE THE PHENOTYPE VISIBLE

For each disease/phenotype cluster, construct a measurement-evidence table. Example:

```text
Phenotype: orthostatic/autonomic dysfunction
Observable signals: HR response; HRV; posture-associated HR change; activity intolerance; sleep disruption
Potential technologies: ECG wearable; PPG wearable; accelerometer; posture sensor; continuous SpO2; autonomic testing
Evidence: public patient datasets; ClinicalTrials.gov; NIH RePORTER; Open Targets / molecular context; public literature-linked omics
```

The engine should calculate separate dimensions: phenotype_signal_strength; measurement_evidence_strength; technology_maturity; regulatory_visibility; deployment_complexity. Do not collapse them immediately into one black-box number.

# PHASE 4 — BUILD GEOGRAPHIC OPPORTUNITY FEATURES

For geography `g`, condition `c`, and measurement `m`, calculate transparent components such as: burden(g,c); vulnerability(g); diagnostic_desert(g,c); clinic_capacity(g,m); research_readiness(g,c,m); technology_saturation(g,m).

diagnostic_desert = high burden + high vulnerability + low relevant-provider density + low relevant-trial availability. Keep all raw components.

Opportunity(g,c,m) = w1·burden + w2·vulnerability + w3·diagnostic_desert + w4·clinic_capacity + w5·research_readiness. Use configurable weights. Default to equal weighting unless a defensible rationale exists. Perform sensitivity analysis across weight choices. The ranking must not depend on one arbitrary weight configuration. If burden is only a proxy, propagate that uncertainty into the result.

# PHASE 5 — CLINIC MATCHING

For every high-opportunity geography: identify candidate facilities within configurable radii. Rank site characteristics separately: clinical_specialty_match; relevant_trial_history; NIH_research_activity; community_access; distance_to_target_population; technology_experience.

The output should NOT say "this is the best doctor." It should say: "This facility has characteristics suggesting it may be a viable implementation or study partner." Provide the underlying reasons.

# PHASE 6 — TECHNOLOGY-TO-CLINIC MATCHING

CONDITION → PHENOTYPE → OBJECTIVE MEASUREMENT → TECHNOLOGY → GEOGRAPHY → CLINIC / RESEARCH SITE.

Example: Long COVID / dysautonomia-like phenotype → orthostatic HR + activity abnormalities → wearable ECG/PPG + accelerometer → candidate FDA-visible wearable class → county/state with high Long COVID burden and low specialty/trial coverage → nearby cardiology/autonomic/rehab/community health infrastructure.

The final statement should be framed as **candidate deployment opportunity**, not **validated diagnostic pathway**.

---

# AGENTIC MCP LAYER

Build an MCP server using FastMCP or an equivalent implementation. Expose at least:

```python
search_condition(condition: str)
normalize_condition(condition_or_code: str)
get_patient_phenotype_signature(condition: str)
get_molecular_context(condition: str)
discover_candidate_measurements(condition: str, phenotype: str | None = None)
get_measurement_evidence(measurement: str, condition: str)
get_regulatory_context(technology: str)
get_condition_burden(condition: str, geography_level: str)
get_geographic_context(geography_id: str)
find_candidate_clinics(geography_id: str, condition: str, measurement: str)
find_relevant_trials(condition: str, measurement: str | None = None)
find_relevant_research_centers(condition: str)
rank_deployment_opportunities(condition: str, measurement: str, geography_level: str = "county")
search_us_open_data(query: str)
trace_evidence(object_id: str)
```

The LLM/agent must never generate facts when a tool has no data. Return `UNKNOWN / NOT AVAILABLE` instead.

# DEPLOYMENT-OPPORTUNITY OUTPUT SCHEMA

```json
{
  "condition": {},
  "phenotype": {},
  "measurement": {},
  "technology": {"name": "", "regulatory_context": "", "technology_evidence": []},
  "geography": {"name": "", "fips": "", "burden": null, "burden_evidence_level": "", "vulnerability": null, "diagnostic_desert": null},
  "candidate_sites": [],
  "research_evidence": [],
  "molecular_context": [],
  "uncertainties": [],
  "provenance": [],
  "recommended_next_step": ""
}
```

# DASHBOARD

Build a simple interactive dashboard (Streamlit, Dash, or a lightweight React/FastAPI interface). Required pages:

1. **National Opportunity Map** — disease burden; diagnostic deserts; deployment opportunities; facilities; trial sites. Allow toggling layers independently.
2. **Condition Explorer** — input `Long COVID`; return phenotype features; wearable abnormalities from public data; molecular evidence; candidate measurements; relevant technologies; prevalence/burden sources.
3. **Measurement Explorer** — input `wearable autonomic monitoring`; return conditions where this is studied; physiological signals; public supporting datasets; relevant trials; FDA device classes; regions with unmet need.
4. **Geography Explorer** — input `San Diego County, California`; return target-condition burden; vulnerability; provider density; relevant specialists; active/completed trials; NIH-funded research; candidate technologies; nearby health centers.
5. **Deployment Recommendation** — example query: "Find 10 U.S. regions where wearable monitoring of autonomic/activity abnormalities in Long COVID or ME/CFS would be useful to evaluate, and identify clinics/research sites that could plausibly participate." Return transparent ranked table; map; component scores; evidence; limitations.

# CAPRIO / CAPILLAROSCOPY ADAPTER

Create an explicit interface:

```python
class MeasurementAdapter:
    def preprocess(self, raw_data): ...
    def embed(self, raw_data): ...
    def phenotype_score(self, embedding): ...
    def describe_measurement(self): ...
```

Implement `WearableAdapter` for the public-data alpha. Create a stub `CapillaroscopyAdapter` with no private data. Document that CAPRIO can later supply: image → microvascular embedding → capillary phenotype → same disease/ontology/geography/clinic pipeline. This is critical to showing that CAPRIO is a validation program built on top of a generalizable infrastructure rather than the infrastructure itself.

# DATA AUDIT REQUIREMENTS

For EVERY data source create `data/raw/<source>/DATA_AUDIT.md`. Include: source; publishing organization; retrieval date; license/access conditions; update date; unit of observation; sample size; geography; key variables; missingness; limitations; linkage strategy; whether person-level; whether geographic; whether omics; whether wearable; whether true participant linkage exists.

Create `SOURCE_REGISTRY.yaml`. The application should read from this registry.

# REPRODUCIBILITY

Use: Python; uv; DuckDB; Parquet; Polars or pandas; GeoPandas where necessary; scikit-learn; FastAPI; FastMCP; Plotly/Leaflet/MapLibre or equivalent; pytest. Cache API requests. Do not repeatedly hit federal services unnecessarily. Record timestamps and source versions.

# VALIDATION AND ANTI-SYCOPHANCY TESTS

The project is unsuccessful if it only produces attractive maps. Explicitly test:

1. **Wearable signal** — Can a wearable-derived phenotype distinguish or meaningfully stratify a relevant condition in at least one public dataset? Compare against demographic-only baselines.
2. **Multimodal increment** — Does adding wearable physiology improve over clinical variables alone? Report both positive and negative results.
3. **Molecular coherence** — Do disease-linked molecular resources return plausible biological systems? Do not judge this by intuition alone. Report source counts and agreement/disagreement across databases.
4. **Geography** — Are hotspots driven by actual burden variables rather than population size alone? Repeat analyses using prevalence and absolute counts where possible.
5. **Deployment** — Does adding facility/research information change where the system would deploy a technology compared with disease burden alone? This is a core product demonstration.
6. **Negative controls** — Shuffle: phenotype labels; geography labels; clinic locations where appropriate. The pipeline should lose meaningful structure.
7. **Ranking robustness** — Re-run opportunity rankings under multiple reasonable weight sets. Report regions whose ranking is unstable.

# REQUIRED FINAL DATASETS

participant_wearable_features.parquet; participant_clinical_features.parquet; condition_registry.parquet; condition_molecular_evidence.parquet; condition_measurement_evidence.parquet; measurement_registry.parquet; geo_condition_features.parquet; clinic_registry.parquet; research_site_registry.parquet; deployment_opportunities.parquet; measure_it_public.duckdb

# REQUIRED FIGURES

Approximately 8–10 polished figures. At minimum:
1. System architecture: Person → phenotype → measurable biology → technology → geography → clinic → deployment
2. Public-data source architecture
3. Example wearable phenotype differences
4. Clinical-only versus clinical + wearable model performance
5. Condition → molecular evidence → measurable-physiology graph
6. U.S. burden / unmet-need map for one target condition
7. Provider/research infrastructure overlay
8. Deployment opportunity map
9. Example geography drill-down
10. Modularity: Wearables today ↕ Measurement Adapter ↕ CAPRIO capillaroscopy tomorrow

# PRIMARY DEMO STORY

Prefer Long COVID / ME/CFS / dysautonomia because public wearable data make that demonstration feasible.

1. Public patient-level wearable data show an objectively measurable activity/autonomic phenotype.
2. Clinical and ontology data connect that phenotype to one or more invisible illnesses.
3. Public molecular resources show associated biological pathways and candidate measurable biology.
4. ClinicalTrials.gov, NIH RePORTER, and openFDA identify technologies and objective measurement strategies already under investigation or deployable.
5. Federal population data identify where relevant burden and vulnerability are concentrated.
6. NPPES, HRSA, ClinicalTrials.gov, and RePORTER show where clinical/research infrastructure does or does not exist.
7. The engine identifies **measurement deserts** and **candidate deployment sites**.
8. An agent can answer: "What should we measure, where should we measure it, and who could realistically deploy it?"
9. CAPRIO is then presented as one future device-specific implementation of the same engine.

# KEY PRODUCT LANGUAGE

Prefer: candidate; evidence-supported; measurable phenotype; objective signal; deployment opportunity; measurement desert; translational site; research readiness; implementation candidate; public-data evidence; physiological phenotype.

Avoid (unless directly supported): diagnosed by AI; proves; cures; confirms mechanism; patient has; definitive biomarker; best clinic; optimal treatment.

# SUCCESS CRITERIA

1. At least one real public patient-level wearable dataset has been processed.
2. Wearable features are linked to real clinical or disease information.
3. At least one invisible-illness example is demonstrated.
4. A disease/phenotype ontology connects clinical concepts to molecular evidence.
5. Candidate objective measurement technologies can be discovered automatically.
6. U.S. public geographic data can be queried and normalized.
7. Providers/facilities/trial sites are mapped.
8. The system generates reproducible geographic deployment opportunities.
9. Every recommendation can be traced back to source data.
10. CAPRIO could later be inserted as a new MeasurementAdapter without restructuring the platform.

# FIRST IMPLEMENTATION PRIORITY

Build one complete vertical slice first: Long COVID / ME/CFS-like phenotype → public wearable physiology → objective phenotype model → Mondo/EFO/HPO → Open Targets/GWAS/GEO molecular context → wearable/autonomic measurement candidates → ClinicalTrials.gov + RePORTER + openFDA → CDC/CMS/ACS/SVI geographic context → NPPES + HRSA + trial sites → deployment-opportunity map. Then generalize.

# FINAL REPORT

Produce `results/FINAL_REPORT.md` including: executive summary; exact datasets used; actual sample sizes; linkage audit; wearable methods; phenotype modeling; molecular enrichment; technology discovery; geographic methods; facility matching; scoring formula; validation; negative results; limitations; screenshots/figures; example MCP queries; example deployment recommendation; explanation of how CAPRIO plugs in later.

Do not describe unfinished components as finished. Do not inflate performance. Do not infer patient-level geography. Do not call condition-level molecular enrichment "personal multi-omics."

The end product should demonstrate:

> **A reusable AI system that translates fragmented invisible-illness phenotypes into objective measurement opportunities and uses U.S. public data to determine where those technologies can be deployed and which clinics or research sites are positioned to act.**
