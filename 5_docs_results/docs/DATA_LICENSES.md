# Data licences and redistribution

This page records, for each of the 37 sources in `SOURCE_REGISTRY.yaml`, the provider's licence, whether **derived
processed tables** (the Parquet/DuckDB tables the pipeline writes to `data/processed`) could be redistributed publicly,
for example as a data-bundle release asset, and the attribution the provider asks for.

* **Current policy: code only.** No data bundle is published. Everyone rebuilds the data locally with
  `uv run measure-it pipeline` (README "Building the data"). The MIT licence in `LICENSE` covers the project's code and
  documentation, **not** any third-party data.
* **What the repository does contain** is described in section 3: provenance sidecars and the derived result tables,
  figures and write-ups of the reference build.
* Checked 2026-09-28 against the providers' own pages (the "evidence" column lists the pages read). The DATA_AUDIT.md
  of each source (`data/raw/<source_id>/DATA_AUDIT.md`) records what the pipeline saw at retrieval time. This is a
  working record, **not legal advice**; where the terms are silent or ambiguous the verdict says **uncertain**.

Verdicts: **yes** (no condition beyond good practice), **yes, with attribution**, **yes, with conditions** (listed),
**no, rebuild locally**, **uncertain**.

## 1. Per-source table

### Person-level sources (layer 1)

| source_id | licence / terms | redistribute derived tables? | required attribution / conditions | evidence |
|---|---|---|---|---|
| `nhanes_2011_2014`, `nhanes_2003_2006` | U.S. public-use files under the NCHS Data User Agreement (PHS Act 308(d), CIPSEA) | **yes, with conditions**: the agreement does not address redistribution; per-SEQN rows containing only public-use variables may be shared if the terms travel with them | Pass on the agreement: statistical reporting and analysis only; no attempt to learn anyone's identity; no linkage with individually identifiable data from other datasets; no re-identification research. Cite "National Center for Health Statistics, NHANES <cycle>" (no official citation text found) and note that use does not imply CDC endorsement | https://www.cdc.gov/nchs/policy/data-user-agreement.html ; https://www.cdc.gov/other/agencymaterials.html |
| `stanford_longcovid_wearables`: Uwakwe 2025 | ODC-By 1.0 (Stanford Digital Repository) | **yes, with attribution** | "Uwakwe, C., Rangan, E., Kumar, S., Brooks, A., Maguire, P., and Snyder, M. (2025). Snyder Lab Long COVID Study Dataset. Version 1. Stanford Digital Repository. Available at https://purl.stanford.edu/cb174pb4851/version/1. https://doi.org/10.25740/cb174pb4851." Use terms: content "will not be used to identify or to otherwise infringe the privacy or confidentiality rights of individuals" | https://purl.stanford.edu/cb174pb4851 ; https://snyderlabs.stanford.edu/long-covid |
| `stanford_longcovid_wearables`: Mishra 2020 | none: public download link in the paper, no licence in or next to the files; the article (Nat Biomed Eng 4:1208) is not open access | **no, rebuild locally** (per-person rows); aggregate results uncertain | Cite Mishra T et al. 2020, doi:10.1038/s41551-020-00640-6. Redistributing per-participant rows would need the Snyder lab's permission | https://www.nature.com/articles/s41551-020-00640-6 |
| `stanford_longcovid_wearables`: Alavi 2022 | article CC BY 4.0 (Nat Med 28:175); raw-data ZIP has no licence | **no, rebuild locally** for rows from the ZIP (uncertain); article Source Data: yes, with attribution | Cite Alavi A et al. 2022, doi:10.1038/s41591-021-01593-2 | https://www.nature.com/articles/s41591-021-01593-2 |
| `mapmecfs_nih_pi_mecfs` (Walitt 2024) | article and supplements CC BY 4.0; Pennsieve dataset 356 "Creative Commons Attribution"; GEO/SRA NCBI policy; **mapMECFS files: data-use agreement** | article supplements, Pennsieve and GEO-derived tables: **yes, with attribution**. mapMECFS-hosted files: **no** ("redistribution of data in any manner is prohibited"); none was downloaded | Cite Walitt B et al. 2024, Nat Commun, doi:10.1038/s41467-024-45107-3, and the Pennsieve DOI doi:10.26275/ile7-wrsk; CC BY: link the licence and indicate changes | https://www.nature.com/articles/s41467-024-45107-3 ; https://api.pennsieve.io/discover/datasets/356 ; https://www.mapmecfs.org/about ; mapMECFS DUA (https://www.mapmecfs.org/files/mapMECFS_DUA_26jan2021.pdf) |
| `charlton_lc_mecfs_cpet_source` (MUSCLE-ME) | CC BY 4.0 (Nat Commun article; data also on Zenodo doi:10.5281/zenodo.20642303) | **yes, with attribution** (see note A) | Cite Charlton et al. 2026, doi:10.1038/s41467-026-75725-y; link CC BY 4.0; indicate changes | https://www.nature.com/articles/s41467-026-75725-y |
| `appelman_lc_pem_source` | CC BY 4.0 (Nat Commun) | **yes, with attribution** (note A) | Cite Appelman B et al. 2024, Nat Commun 15:17, doi:10.1038/s41467-023-44432-3 | https://www.nature.com/articles/s41467-023-44432-3 |
| `klein2023_mylc_ml_table` (MY-LC) | CC BY 4.0 (Nature) | **yes, with attribution** (note A) | Cite Klein J et al. 2023, Nature 623:139-148, doi:10.1038/s41586-023-06651-y | https://www.nature.com/articles/s41586-023-06651-y |
| `fm_thermography` | CC BY 4.0 (PLOS ONE; S1 Data) | **yes, with attribution** | Cite Sempere-Rubio N et al. 2021, PLOS ONE 16(6):e0253281, doi:10.1371/journal.pone.0253281 ("provided the original author and source are credited") | https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0253281 ; https://journals.plos.org/plosone/s/licenses-and-copyright |
| `endo_arg1_repod` | **CC BY 4.0** (set per file on RepOD; the Dataverse API's dataset-level `license` field is null, which is why the registry says "not stated") | **yes, with attribution** | "Czystowska-Kuzmicz, Malgorzata, 2024, "Endometriosis patients cohort - ARG1 ELISA and activity data and clinical data", https://doi.org/10.18150/PY1P9X, RepOD, V1"; also cite Pliszkiewicz M et al., J Clin Med 2024;13:1489 | https://repod.icm.edu.pl/dataset.xhtml?persistentId=doi:10.18150/PY1P9X |
| `heds_hsd_olink_serum_cinquina2026` | **CC BY-NC-ND 4.0** (Clinical Proteomics 2026, PMC13081554), **not** CC BY 4.0 as `SOURCE_REGISTRY.yaml` and its DATA_AUDIT.md say | **no, rebuild locally** ("You do not have permission under this licence to share adapted material derived from this article or parts of it"; non-commercial only). Whether group-level statistics computed from the file count as adapted material is **uncertain** | Cite Cinquina V et al. 2026, Clinical Proteomics (PMC13081554) | Europe PMC full-text XML of PMC13081554 (license element) |
| `geo_cohorts` | NCBI policy: NCBI places no restrictions; submitters may claim rights on individual series | **yes, with conditions**: check each GEO series (GSE) for submitter terms | Cite each GSE accession and its publication | https://www.ncbi.nlm.nih.gov/home/about/policies/ |

### Condition-level molecular, ontology and literature sources (layer 2 and reference)

| source_id | licence / terms | redistribute derived tables? | required attribution / conditions | evidence |
|---|---|---|---|---|
| `open_targets` | data CC0 1.0 (code Apache-2.0) | **yes** | Requested: cite the latest Open Targets Platform publication; component sources keep their own licences | https://platform-docs.opentargets.org/licence |
| `gwas_catalog` | EMBL-EBI terms of use (only summary statistics are CC0; graphics CC BY 4.0) | **yes, with attribution** | "attribution ... in accordance with good scientific practice": cite the GWAS Catalog (Sollis E et al., NAR 2023) and the release date; EBI adds no restriction beyond those of the original data owners | https://www.ebi.ac.uk/gwas/docs/about ; https://www.ebi.ac.uk/about/terms-of-use |
| `ncbi_geo_sra` (metadata only) | NCBI policy (as `geo_cohorts`) | **yes** | Cite accessions | https://www.ncbi.nlm.nih.gov/home/about/policies/ |
| `reactome` | data **CC0** (illustrations and icons CC BY 4.0; software Apache-2.0); the registry's "CC BY 4.0" is out of date | **yes** | "Attribution is encouraged but not required": cite Reactome release 97 | https://reactome.org/license |
| `hgnc` | CC0 | **yes** | Recommended: Seal RL et al., "Genenames.org: the HGNC resources in 2023", NAR; RRID:SCR_002827 | https://www.genenames.org/about/license/ |
| `mondo_ontology`: Mondo | CC BY 4.0 | **yes, with attribution** | Cite Mondo with its release version; link CC BY 4.0 | https://mondo.monarchinitiative.org/ ; Mondo GitHub LICENSE |
| `mondo_ontology`: Monarch KG disease-phenotype associations | no single licence; each component source's licence applies (HPO annotations: HPO licence) | **uncertain** (per component source) | Follow the HPO row for HPO annotations | Monarch app terms page and docs (GitHub monarch-initiative); KG metadata.yaml has no licence field |
| HPO terms and annotations (via `mondo_ontology`, `ebi_ols4`) | HPO licence (custom) | **yes, with conditions**: the HPO content and its logical relationships may not be altered; whether filtered or joined subsets count as altered is **uncertain**, so label them as derived and state the version | "This service/product uses the Human Phenotype Ontology (version information). Find out more at http://www.human-phenotype-ontology.org"; cite Gargano MA et al., NAR 2024, doi:10.1093/nar/gkad1005; show the HPO version/date | https://hpo.jax.org/license ; `hp.obo` header |
| `ebi_ols4` (EFO, MeSH, DOID, ORDO lookups) | OLS: EMBL-EBI terms of use; each ontology its own licence: EFO Apache-2.0; DOID CC0; ORDO CC BY 4.0; MeSH NLM terms | **yes, with attribution** | EFO: keep the Apache notice; ORDO: "Orphadata Science: Free access data from Orphanet. © INSERM 1999. Available on http://sciences.orphadata.com/. Data version [date]." and indicate changes; MeSH: "Courtesy of the U.S. National Library of Medicine", no implied endorsement, state that the data may not be current | https://www.ebi.ac.uk/about/terms-of-use ; https://www.orphadata.com/legal-notice/ ; https://www.nlm.nih.gov/databases/download/mesh.html ; `efo.obo`, `doid.obo` headers |
| `cms_icd10cm` | U.S. government file; ICD-10-CM is NCHS's clinical modification of WHO's ICD-10 (WHO copyright) | **uncertain**: widely redistributed, but neither the CMS nor the NCHS page states a reuse licence | Cite CMS/NCHS ICD-10-CM <fiscal year> | https://www.cms.gov/medicare/coding-billing/icd-10-codes ; https://www.cdc.gov/nchs/icd/icd-10-cm/index.html |
| `published_device_evidence`, `published_lab_evidence` | abstracts and full texts stay under publisher/author copyright; Europe PMC metadata has no stated licence | **yes, with conditions**: short verbatim quotes with citations rely on fair use / fair dealing, not a licence; do not redistribute full abstracts or texts | Cite every article (the tables carry DOI/PMID) | https://dev.europepmc.org/Copyright |

### Population / geographic sources (layer 3)

| source_id | licence / terms | redistribute derived tables? | required attribution / conditions | evidence |
|---|---|---|---|---|
| `cdc_places` | Public domain (data.cdc.gov metadata) | **yes, with attribution** | "Centers for Disease Control and Prevention, National Center for Chronic Disease Prevention and Health Promotion, Division of Population Health"; CDC general terms: attribute the agency, do not imply endorsement by CDC, ATSDR, HHS or the U.S. Government, do not change substantive content | https://data.cdc.gov/api/views/swc5-untb.json ; https://www.cdc.gov/other/agencymaterials.html |
| `cdc_svi` | Public domain (CDC/ATSDR) | **yes, with attribution** | "Centers for Disease Control and Prevention/ Agency for Toxic Substances and Disease Registry/ Geospatial Research, Analysis, and Services Program. CDC/ATSDR Social Vulnerability Index 2022 Database US." plus access date; CDC general terms | https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html |
| `cdc_long_covid` (Household Pulse) | "Public Domain U.S. Government" (data.cdc.gov; some tables show licence "None") | **yes, with attribution** | Cite NCHS, Household Pulse Survey long COVID estimates; CDC general terms | data.cdc.gov catalog metadata ; https://www.cdc.gov/other/agencymaterials.html |
| `cdc_lyme` | "Public Domain U.S. Government" | **yes, with attribution** | Cite CDC Lyme disease surveillance data; CDC general terms | data.cdc.gov catalog metadata |
| `census_acs`, `census_geography` | U.S. government work | **yes, with attribution** | Cite "the Census Bureau as the source of the original data only", e.g. U.S. Census Bureau, American Community Survey 5-year estimates 2020-2024, Summary File; U.S. Census Bureau, cartographic boundary files 2024; no implied endorsement; do not present modified content as the Census Bureau's | https://www.census.gov/about/policies/citation.html ; https://www.census.gov/data/developers/about/terms-of-service.html |
| `usda_rucc` | U.S. government work | **yes, with attribution** | "U.S. Department of Agriculture, Economic Research Service. (January 2024). Rural-Urban Continuum Codes." | https://www.ers.usda.gov/data-products/rural-urban-continuum-codes/documentation |
| `cms_mmd` | U.S. government work (aggregate, cell-suppressed) per the registry; the data.cms.gov terms page could not be machine-read | **uncertain** (very likely yes: federal aggregate data) | Cite CMS Office of Minority Health, Mapping Medicare Disparities Tool | https://data.cms.gov/tools/mapping-medicare-disparities-by-population ; https://data.cms.gov/terms-of-use (JavaScript-only) |
| `data_gov_catalog` (metadata only) | federal content not under copyright (17 U.S.C. 105); non-federal records may differ | **yes** | Requested: "Data retrieved from Data.gov (https://www.data.gov/)" | https://data.gov/privacy-policy/ |

### Facility, provider and research sources (layer 4)

| source_id | licence / terms | redistribute derived tables? | required attribution / conditions | evidence |
|---|---|---|---|---|
| `clinicaltrials_gov` | U.S. government database; NLM terms (some third-party copyright; copyright outside the U.S.) | **yes, with conditions** | "Courtesy of the U.S. National Library of Medicine"; attribute the source as ClinicalTrials.gov; display the date the data were processed; state any modifications; state conspicuously that the data may not be current; no implied NLM endorsement; no marketing use of contact e-mails | https://clinicaltrials.gov/about-site/terms-conditions ; https://www.nlm.nih.gov/databases/download/terms_and_conditions.html |
| `nih_reporter` | RePORTER FAQ calls the data "available in the public domain"; no licence page; project abstracts are written by grantees | **yes for project metadata; abstracts uncertain** | Cite NIH RePORTER with the retrieval date | https://report.nih.gov/faqs ; https://api.reporter.nih.gov/ |
| `openfda_device` | CC0 1.0 | **yes** | None required; do not imply FDA endorsement; GMDN terms inside device records need their own licence for commercial use | https://open.fda.gov/license/ |
| `nppes` | FOIA-disclosable CMS data | **yes** (no restriction found) | Cite CMS NPPES NPI Registry, <month> file | https://download.cms.gov/nppes/NPI_Files.html ; https://www.cms.gov/medicare/regulations-guidance/administrative-simplification/data-dissemination |
| `cdc_brfss` | U.S. government work (CDC public-use survey file; no restriction on the public LLCP files) | **yes** | Cite CDC Behavioral Risk Factor Surveillance System Survey Data, <year>; the small-area estimates are this project's model, not CDC estimates | https://www.cdc.gov/brfss/annual_data/annual_data.htm |
| `cms_physician_service` | U.S. government work (aggregate, CMS suppresses cells < 11 beneficiaries) | **uncertain** (very likely yes: federal aggregate data; data.cms.gov terms page is JavaScript-only) | Cite CMS, Medicare Physician & Other Practitioners - by Provider and Service, <data year> | https://data.cms.gov/provider-summary-by-type-of-service/medicare-physician-other-practitioners/medicare-physician-other-practitioners-by-provider-and-service ; https://data.cms.gov/terms-of-use |
| `cms_partd_prescriber` | U.S. government work (aggregate, CMS suppresses counts < 11) | **uncertain** (very likely yes: federal aggregate data) | Cite CMS, Medicare Part D Prescribers - by Provider and Drug, <data year> | https://data.cms.gov/provider-summary-by-type-of-service/medicare-part-d-prescribers/medicare-part-d-prescribers-by-provider-and-drug |
| `nucc_taxonomy` | © American Medical Association; NUCC permission to use and distribute (royalty-free, U.S.) | **yes, with conditions (restrictive)**: no changes except format/style and **no derivative works**; keep the AMA copyright notice and all notices and disclaimers; commercial use needs a licence. Safest: ship taxonomy **codes** only, not an edited copy of the code table | Keep the AMA copyright notice and NUCC disclaimers | https://www.nucc.org/index.php/code-sets-mainmenu-41/provider-taxonomy-mainmenu-40/csv-mainmenu-57 and its "Permission to use and distribute" page |
| `hrsa_health_centers`, `hrsa_hpsa` | U.S. government data; "Usage limitations: None" | **yes** | Cite HRSA Data Warehouse, retrieval date | https://data.hrsa.gov/data/download |

**Note A (journal supplements).** Nature and Nature Communications apply CC BY 4.0 to "this article" and material
"included in the article"; no explicit statement was found that Supplementary Information and Source Data files are
covered by the same licence. Treat them as CC BY 4.0 with attribution (the usual reading), but this is **uncertain**.

**Discovery probes** (`data/raw/device_candidates/*`, `data/raw/lab_candidates/*`, `data/raw/snyder_candidates/*`)
are not registry sources: files were downloaded only to verify candidates, and only their MANIFEST.json is tracked.
Each candidate's licence is in `results/tables/device_dataset_candidates.csv` and `lab_dataset_candidates.csv`
(`licence` column). The three `excluded_privacy` records were never downloaded into `data/raw`.

## 2. What a data bundle could contain (for a future decision)

If a bundle (`uv run measure-it bundle`) is ever published, the verdicts above imply:

| can be included | leave out (rebuild locally) or get permission first |
|---|---|
| tables from U.S. federal sources (CDC, Census, USDA, HRSA, CMS NPPES, ClinicalTrials.gov with its credit and currency notice, RePORTER metadata, openFDA, Data.gov) | `participant_*__stanford_covid` rows from Mishra 2020 and Alavi 2022 (no licence) |
| NHANES per-person tables, with the NCHS Data User Agreement terms in the bundle | any table built from the Cinquina 2026 file (CC BY-NC-ND) |
| Uwakwe 2025 rows (ODC-By, with the preferred citation) | an edited copy of the NUCC taxonomy table (no derivative works; codes alone are fine) |
| CC BY supplements (Charlton, Appelman, Klein, Walitt/Pennsieve, PLOS fibromyalgia thermography, RepOD ARG1), with citations and a change notice | HPO-derived tables unless labelled as derived with the HPO version and citation (uncertain) |
| Open Targets, HGNC, Reactome, DOID (CC0), Mondo/ORDO (CC BY), EFO (Apache-2.0), GWAS Catalog (EBI terms) | ICD-10-CM text and CMS MMD tables until their terms are confirmed (uncertain) |

A published bundle would also need an `ATTRIBUTION.md` inside it collecting the attribution texts of section 1, and a
statement of the retrieval dates (ClinicalTrials.gov and MeSH require a "may not be current" notice).

## 3. What the repository contains today

| path | content | licence position |
|---|---|---|
| `data/raw/**/{DATA_AUDIT.md,MANIFEST.json,registry_entry.yaml,QUERY_LOG.json,queries.json}`, `data/processed/*.meta.json`, `SOURCE_REGISTRY.yaml` | provenance written by this project (URLs, sizes, hashes, counts, audit prose) | project documentation (MIT); no third-party data rows |
| `results/*.md`, `results/tables/*.csv|json`, `results/figures`, `results/maps` | results of the reference build: aggregate statistics, rankings of counties and facilities from federal data, molecular enrichment, short quoted claims with citations | derived from the sources above; attribution per section 1 applies. Four tables hold **per-person derived values** (participant id + model output or detection window, no demographics or locations): `stanford_uwakwe_oof_predictions.csv` (Uwakwe, ODC-By: fine with attribution), `fm_thermography_oof_predictions.csv` (PLOS CC BY: fine), and `stanford_acute_participants.csv`, `stanford_acute_detection_windows.csv` (Mishra 2020 / Alavi 2022: **no licence; remove before the repository goes public, or ask the Snyder lab**). The six `heds_hsd_olink_serum_cinquina2026_*` tables are group-level statistics from a CC BY-NC-ND file (**uncertain**; see docs/GITHUB_RELEASE_PLAN.md) |
| `results/maps/plotly.min.js` | Plotly.js bundle | MIT (its own header) |

## 4. Corrections this audit found in the source registry

To be fixed in the modules that write the registry fragments (then `uv run measure-it merge-registry`):

* `heds_hsd_olink_serum_cinquina2026`: licence is **CC BY-NC-ND 4.0**, not "CC BY 4.0 (article and additional files,
  BMC open access)" (`src/measure_it/labs/heds_hsd_olink_serum_cinquina2026.py`, its registry entry and DATA_AUDIT.md).
* `reactome`: data are **CC0** (illustrations CC BY 4.0), not "CC BY 4.0 (Reactome data)"
  (`src/measure_it/omics/reference_data.py`).
* `endo_arg1_repod`: licence is **CC BY 4.0** (per file on RepOD), not "not stated".
* `nucc_taxonomy`: the terms forbid derivative works and changes other than format; the registry's "free public
  download; NUCC terms of use" understates this.
* `gwas_catalog`: only summary statistics are CC0; the catalogue itself is under the EMBL-EBI terms of use (as the
  registry says).
