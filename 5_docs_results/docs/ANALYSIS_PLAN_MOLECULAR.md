# Analysis plan — molecular coherence (Test 3) and the phenotype -> measurable-biology graph

Written 2026-09-23 **before** any agreement, enrichment or control statistic was computed.
Modules: `src/measure_it/omics/coherence.py` (Test 3), `src/measure_it/omics/graph.py`
(`measurable_biology`, Figure 5 data). Outputs: `results/tables/test3_*.csv`,
`results/tables/fig5_*.csv`, `results/figures/drafts/fig5_*`, `results/MOLECULAR_COHERENCE.md`,
processed table `measurable_biology`. Deviations made after this plan was written are listed in
the last section; nothing above that section is edited after results were seen.

Everything here is **condition-level molecular enrichment**: gene lists that public databases attach
to a disease concept. Nothing is measured on, or linked to, any participant; no result is a
statement about a person, and no combined "omics score" is built.

## 1. Question

SPEC validation test 3: *Do disease-linked molecular resources return plausible biological systems?
Report source counts and agreement/disagreement across databases; do not judge by intuition.*

Operationalised as three questions, each with a negative control:

* **Q1 gene level** — for the same condition, do independent public sources name overlapping genes
  more than the same sources do for *different* conditions?
* **Q2 pathway level** — do the Reactome pathways enriched in each source's gene list agree across
  sources for the same condition more than across conditions?
* **Q3 system level** — which physiological systems (immune, neuronal, vascular/endothelial,
  autonomic/cardiac, mitochondrial/metabolic, connective tissue, other) are enriched per condition,
  are they condition-specific, and do they match the condition's curated clinical flags more than
  shuffled flags do?

## 2. Inputs (all read with `store.read_table`)

| table | use |
|---|---|
| `condition_registry` | conditions, primary Mondo id, `grouping_only`, clinical flags (Q3 control) |
| `condition_molecular_evidence` (+ partitions) | gene sets (below) |
| `condition_molecular_coverage` | why a source has no rows (id absent vs zero hits) |
| `geo_study_catalog`, `sra_study_catalog` | source counts only (study metadata, no genes) |

Reference data fetched through `measure_it.http` (disk-cached; version recorded in every output):

* Reactome (current release; version from `ContentService/data/database/version`):
  `ReactomePathways.gmt.zip` (gene-symbol membership of every human pathway, hierarchy-inclusive),
  `ReactomePathways.txt`, `ReactomePathwaysRelation.txt` (parent -> child), and the
  `AnalysisService/identifiers/projection` API for a cross-check only (section 5.4).
* HGNC complete set (`hgnc_complete_set.txt`) for symbol harmonisation and the gene background.

Audit traps taken into account (from the DATA_AUDIT files of open_targets, gwas_catalog,
ncbi_geo_sra, mapmecfs_nih_pi_mecfs):

* Open Targets association calls are **indirect by default** (descendant diseases included). Broad or
  component ids (autonomic nervous system disorder, Ehlers-Danlos syndrome, post-infectious
  disorder) are dominated (80-100 %) by descendant diseases such as neoplasms and monogenic subtypes.
* **POTS** maps in Open Targets only to MONDO:0011479, which Mondo labels *POTS due to NET
  deficiency* (monogenic, predicate narrow); Open Targets labels it plain POTS. Results for `pots`
  carry this caveat and are repeated without `pots` in the sensitivity analysis.
* **MCAS** has no identifier present in Open Targets or the GWAS Catalog and zero GEO/SRA hits:
  it is reported as "no molecular evidence available", never imputed.
* GWAS Catalog includes sub-threshold (p < 1e-5) associations and multi-trait (pleiotropy/MTAG)
  analyses that attach the same locus to two conditions (e.g. "Endometriosis or migraine").
* Open Targets `genetic_association` includes `gwas_credible_sets`, which are derived from GWAS
  studies that are also in the GWAS Catalog: OT-genetic vs GWAS-Catalog agreement is **partly shared
  underlying data**, not independent replication.
* Open Targets literature evidence is Europe PMC co-mention text mining (a study-attention signal).
* mapMECFS (Walitt 2024) differential genes/proteins: **no gene- or protein-level row passes FDR**
  (fdr_significant is False or UNKNOWN for every such row); the lists are nominal p < 0.05 results
  from ~14-18 per group, and the PBMC/muscle tables list only significant rows (the tested universe
  is not published). Some gene symbols in the supplementary tables were converted to dates by the
  spreadsheet (e.g. "2022-03-03"); these are recovered only through an Ensembl id when present,
  otherwise dropped and counted.
* GEO/SRA rows are study metadata from text queries (dysautonomia GEO includes familial dysautonomia,
  MSA and Parkinson's series): counted, never used as gene evidence.

## 3. Gene sets (unit = condition x source; genes = HGNC approved protein-coding symbols)

Harmonisation: Ensembl gene id -> HGNC symbol when an Ensembl id is given (Open Targets, mapMECFS
RNA-seq); otherwise symbol -> approved symbol, else a *unique* previous symbol, else a *unique*
alias; ambiguous or unmapped symbols are dropped and counted per source. Only HGNC
`locus_group == "protein-coding gene"` is kept (non-coding genes are not comparably represented
across the sources or in Reactome); the number dropped is reported.

| source id | definition (primary) |
|---|---|
| `ot_genetic` | Open Targets association rows (`overall_association`) with `ot_datatype_score__genetic_association > 0` |
| `gwas_catalog` | GWAS Catalog variant rows: every gene in `mapped_genes` (all curated associations) |
| `ot_literature` | Open Targets association rows with `ot_datatype_score__literature > 0` (Europe PMC co-mention) |
| `ot_other` | Open Targets association rows with any of clinical (drug targets), animal_model, rna_expression, genetic_literature, somatic_mutation, affected_pathway > 0 |
| `mapmecfs` | me_cfs only: gene-level analytes (PBMC + muscle RNA-seq, serum + CSF SomaScan) from the all-participant tables SD16A, SD19A, SD17A, SD17B with published nominal p < 0.05 (no FDR-significant rows exist) |

Open Targets identifier scope (primary): rows whose registry predicate (`ontology_match`) is exact
or narrow (broad ids excluded), and, for ids with descendant diseases in Open Targets
(`n_descendants_in_source > 0`), only targets with `has_direct_association == True`.
Note: the per-datatype scores are from the indirect call; for ids with descendants a target kept
by the direct flag could still owe its *genetic* score to a descendant — the sensitivity analysis
with all rows brackets this.

Conditions: all 14 registry conditions are reported. A (condition, source) set enters overlap
statistics when it is non-empty and Reactome enrichment when it has >= 5 genes in the Reactome
background; smaller sets are reported as "too few genes" (UNKNOWN), not as negative results.
`post_infectious_syndrome` (grouping node, broad id only) has no primary-scope rows and appears only
in the sensitivity analysis.

## 4. Q1 gene-level agreement

For each condition and each unordered source pair (A, B), both non-empty:
overlap k, |A|, |B|, Jaccard, expected E = |A||B|/N, log2 fold = log2((k + 0.5)/(E + 0.5)),
one-sided hypergeometric p = P(X >= k). BH-FDR across all same-condition tests.

* **Background (stated)**: N = HGNC approved protein-coding genes (primary). Sensitivity: N = union
  of protein-coding genes named by any source for any condition in this project (smaller, more
  conservative).
* **Control (cross-condition baseline)**: the same statistics for every pair (A of condition c,
  B of condition c'), c != c'. Per source pair: mean same-condition log2 fold vs mean cross-condition
  log2 fold; **condition-label permutation test**: statistic T = mean over conditions of the
  same-condition log2 fold; null = T after permuting the condition labels of source B among the
  conditions where both sets are non-empty (exact enumeration when n! <= 50,000, otherwise 20,000
  random permutations with `config.SEED`); p = (1 + #{T_perm >= T_obs}) / (1 + n_perm) (exact:
  proportion including the identity). 95 % CI of the mean difference (same - mean cross, per
  condition) by bootstrap over conditions (2,000 resamples, SEED). BH across the source pairs.
* **Per condition**: empirical rank of the same-condition log2 fold among the 2(n-1) cross pairs
  that share one side with it.
* **mapMECFS** (one condition only): rank of me_cfs among all conditions when the fixed mapMECFS
  set is compared with each condition's `ot_*` / `gwas_catalog` set.

Primary statistic is the continuity-corrected log2 fold (size-adjusted); Jaccard is reported
alongside because it is easy to read but depends on set sizes.

## 5. Q2 pathway-level agreement (Reactome)

### 5.1 Over-representation analysis (ORA)
Library: Reactome GMT, human, symbols harmonised to HGNC protein-coding; pathways with 10-500
genes after harmonisation are tested. Background = union of harmonised genes over all Reactome
pathways (the Reactome-annotated protein-coding genome); input sets are intersected with it.
One-sided hypergeometric test; BH-FDR within each (condition, source) set. A pathway is
**enriched** when FDR < 0.05 and overlap k >= 2.

### 5.2 Agreement
For each condition and source pair where both ORAs ran: Spearman rho of -log10 p across pathways
tested in both (primary); overlap and Jaccard of the enriched sets (secondary; hypergeometric p
against the tested-pathway universe reported but flagged anti-conservative because Reactome
pathways are nested). Control: identical to 4 (cross-condition pairs, label permutation of source B,
bootstrap CI, BH across source pairs).

### 5.3 Physiological-system classification (Reactome hierarchy)
Each pathway's ancestors (ReactomePathwaysRelation) are walked to the top-level pathways. Systems
are assigned from the **anchor table below**: if any ancestor-or-self is a *sub-pathway anchor*,
the pathway takes the system(s) of those anchors; otherwise the system(s) of its top-level
ancestor(s); anything else is `other`. A pathway with several parents can belong to several
systems (counted in each). The top-level rows follow Reactome's own top level; the sub-pathway rows
are curated additions (Reactome has no autonomic top level, so autonomic/cardiac depends entirely
on them). Sensitivity: classification from the top-level rows only.
Cross-check: the hierarchy-derived top level vs the `reactome_top_level_term` Open Targets reports
for the same pathway ids (percent agreement).

<!-- REACTOME_SYSTEM_ANCHORS:BEGIN -->
| reactome_id | reactome_name | level | physiological_system | rationale |
|---|---|---|---|---|
| R-HSA-168256 | Immune System | top | immune | Reactome top level for innate/adaptive immunity and cytokine signalling |
| R-HSA-112316 | Neuronal System | top | neuronal | Reactome top level for synaptic transmission, neuronal ion channels |
| R-HSA-109582 | Hemostasis | top | vascular_endothelial | platelet activation, coagulation, cell-surface interactions at the vascular wall, NO/cGMP platelet homeostasis |
| R-HSA-1474244 | Extracellular matrix organization | top | connective_tissue | collagen/elastic-fibre formation and matrix degradation |
| R-HSA-1430728 | Metabolism | top | mitochondrial_metabolic | Reactome top level for intermediary and energy metabolism |
| R-HSA-5576891 | Cardiac conduction | sub | autonomic_cardiac | cardiac action potential, pacemaking and ion homeostasis (under Muscle contraction) |
| R-HSA-390696 | Adrenoceptors | sub | autonomic_cardiac | sympathetic effector receptors (under Signal Transduction) |
| R-HSA-390648 | Muscarinic acetylcholine receptors | sub | autonomic_cardiac | parasympathetic effector receptors (under Signal Transduction) |
| R-HSA-181430 | Norepinephrine Neurotransmitter Release Cycle | sub | autonomic_cardiac | sympathetic noradrenergic transmission (under Neuronal System) |
| R-HSA-209905 | Catecholamine biosynthesis | sub | autonomic_cardiac | synthesis of the sympathetic transmitters (under Metabolism) |
| R-HSA-202131 | Metabolism of nitric oxide: NOS3 activation and regulation | sub | vascular_endothelial | endothelial NO synthase (under Metabolism) |
| R-HSA-194138 | Signaling by VEGF | sub | vascular_endothelial | endothelial growth/permeability signalling (under Signal Transduction) |
| R-HSA-1592230 | Mitochondrial biogenesis | sub | mitochondrial_metabolic | under Organelle biogenesis and maintenance |
| R-HSA-5205647 | Mitophagy | sub | mitochondrial_metabolic | under Autophagy |
| R-HSA-1268020 | Mitochondrial protein import | sub | mitochondrial_metabolic | under Protein localization |
| R-HSA-5368287 | Mitochondrial translation | sub | mitochondrial_metabolic | under Metabolism of proteins |
| R-HSA-8949215 | Mitochondrial calcium ion transport | sub | mitochondrial_metabolic | under Transport of small molecules |
| R-HSA-9675108 | Nervous system development | sub | neuronal | axon guidance and neural development (under Developmental Biology) |
<!-- REACTOME_SYSTEM_ANCHORS:END -->

All other top levels (Signal Transduction, Gene expression, Metabolism of proteins/RNA, Cell Cycle,
DNA repair/replication, Chromatin organization, Developmental Biology, Disease, Drug ADME,
Transport of small molecules, Vesicle-mediated transport, Cellular responses to stimuli, Programmed
Cell Death, Autophagy, Organelle biogenesis, Protein localization, Cell-Cell communication, Circadian
clock, Digestion and absorption, Reproduction, Sensory Perception, Muscle contraction outside
Cardiac conduction) map to `other`. Sensory Perception is deliberately `other`: it is dominated by
olfactory-receptor clusters that GWAS loci hit positionally, and no candidate measurement class
targets it. Disease pathways are `other` (they re-annotate normal processes in disease states).

### 5.4 Reactome AnalysisService cross-check
Three pre-specified sets are also submitted to the Reactome AnalysisService (projection to human,
no interactors): me_cfs/ot_genetic, migraine/gwas_catalog, long_covid/ot_literature. Reported:
Spearman rho between local and service p-values over shared pathways, and overlap of FDR < 0.05
sets. Informational (validates the local implementation); it does not gate any result.

## 6. Q3 system-level enrichment

Per (condition, source, system): gene-level hypergeometric test of the source set against the
system gene universe (union of genes of all pathways assigned to the system), within the Reactome
background; BH across all (condition, source, system) tests; plus the number of enriched pathways
(5.1) in the system. A system is **supported by a source** when system-level FDR < 0.05 **and**
the source has >= 1 enriched pathway in that system.

* **Condition specificity**: for each source and system, the percentile of the condition's log2 fold
  among all conditions with a tested set for that source; `condition_specific` = above the median of
  the other conditions.
* **Plausibility control (curated flags vs data)**: expected systems per condition from
  `condition_registry` flags, fixed here before looking at results:
  immune <- infectious_or_post_infectious OR autoimmune; autonomic_cardiac <- autonomic;
  vascular_endothelial <- vascular; neuronal <- neurologic; mitochondrial_metabolic <- fatigue_pem
  (connective_tissue has no flag and is excluded). Statistic: total over conditions of |expected ∩
  observed| (observed = systems supported by >= 1 source), and the same with observed restricted to
  condition-specific support. Null: permute which condition's observed profile goes with which flag
  profile (exact or 20,000 permutations, SEED). The flags are themselves curated clinical
  assumptions (`configs/conditions.yaml`), so this is a consistency check, not ground truth.

Support tiers used downstream: `multi_source` (>= 2 source families: {ot_genetic, gwas_catalog} count
as one family "genetic" because they share GWAS data; ot_literature; ot_other; mapmecfs),
`single_source`, `literature_only`, `not_supported`, `too_few_genes`/`no_molecular_data` (UNKNOWN).

## 7. Measurable biology (`measurable_biology`, processed table)

One row per condition x physiological system x candidate measurement class. The
condition -> system link is **data-derived** (section 6); the system -> measurement-class link is
**curated** (table below, `evidence_type = curated_config` for that link). Rows are written for every
condition and every mapped (system, class) pair, including `not_supported` and
`no_molecular_data`, so absence is explicit. A second, data-derived row family links a condition to
an omics measurement class when public study catalogues show that the class has already been
applied to the condition (GEO in-scope series by type -> transcriptomics / proteomics; mapMECFS
published analytes -> transcriptomics / proteomics / metabolomics / immune_assays); these carry
`physiological_system = cross_system_assay_platform` and are research-readiness context, not
molecular support.

<!-- SYSTEM_MEASUREMENT_MAP:BEGIN -->
| physiological_system | measurement_class | rationale |
|---|---|---|
| immune | immune_assays | immune-cell phenotyping, autoantibody and activation assays measure the immune processes these pathways describe |
| immune | proteomics | cytokine, complement and immunoglobulin proteins are quantifiable in plasma/CSF by affinity proteomics |
| immune | blood_biomarkers | CRP, cytokines, complement and cell counts are standard clinical-laboratory readouts of immune activity |
| vascular_endothelial | endothelial_function | FMD / peripheral arterial tonometry measure NO-dependent endothelial vasodilation |
| vascular_endothelial | microvascular_function | laser Doppler / NIRS reactivity measures microvascular regulation |
| vascular_endothelial | capillaroscopy | nailfold capillaroscopy images capillary density/morphology (the future CAPRIO adapter) |
| vascular_endothelial | blood_biomarkers | coagulation and platelet-activation markers (D-dimer, fibrinogen, micro-clot assays) read out hemostasis pathways |
| autonomic_cardiac | hrv | heart-rate variability indexes cardiac vagal/sympathetic modulation |
| autonomic_cardiac | ecg_ambulatory | ambulatory ECG measures rate, rhythm and conduction |
| autonomic_cardiac | autonomic_testing | tilt / active stand / Valsalva / QSART test autonomic reflexes directly |
| autonomic_cardiac | wearable_heart_rate | resting and posture-related heart-rate changes are the wearable-accessible autonomic signal |
| mitochondrial_metabolic | cpet | VO2 peak and ventilatory threshold measure whole-body oxidative capacity |
| mitochondrial_metabolic | metabolomics | direct measurement of the products of metabolic pathways |
| neuronal | small_fiber_testing | skin-biopsy / corneal confocal quantification of small peripheral nerve fibres |
| neuronal | digital_cognitive | reaction time / processing speed as a functional readout of central neuronal systems |
| connective_tissue | vascular_imaging | echocardiographic aortic-root / vascular screening is the standard objective measurement in heritable connective-tissue disorders (weak link; no connective-tissue class exists in configs/measurements.yaml) |
<!-- SYSTEM_MEASUREMENT_MAP:END -->

`other` maps to no measurement class. Measurement class ids must exist in
`configs/measurements.yaml` (tested).

`object_id = molbio:<condition_id>|<system>|<measurement_class>` (namespace to be added to
CONVENTIONS section 7 by its owner; reported, not edited).

## 8. Figure 5 data

`fig5_nodes.csv` / `fig5_edges.csv`: condition -> source evidence (top 5 genes per source by the
source's own statistic: Open Targets datatype score, GWAS minimum p, mapMECFS minimum p; top 3
enriched pathways per source by FDR) -> physiological system (Reactome hierarchy) -> measurement
class (curated). Edge attributes carry the evidence type (data-derived vs curated), source and
counts. Draft figures in `results/figures/drafts/fig5_*`: a layered graph for the demo cluster
(long_covid, me_cfs, pots) and a condition x system heatmap of supporting sources.

## 9. Sensitivity analyses (pre-specified)

S1 Open Targets indirect scope (all ids incl. broad, no direct filter; adds post_infectious_syndrome).
S2 GWAS genome-wide significant only (p < 5e-8). S3 GWAS single-trait only (drop reported traits
matching pleiotropy / MTAG / "or" / "and/or"). S4 gene background = project gene universe. S5 without
`pots`. S6 mapMECFS including the sex-stratified tables (SD16B/C, SD19B/C). S7 system classification
from top-level anchors only. For each: the Q1/Q2 permutation p per source pair and the number of
supported condition-system links.

## 10. What would count as a negative result (stated in advance)

* Same-condition agreement not above the cross-condition baseline for a source pair (permutation
  p >= 0.05) = no condition-specific coherence between those sources, reported as such.
* A system "supported" for (nearly) every condition by a source but not condition-specific =
  source-level bias (e.g. literature attention to immune genes), not disease biology.
* Flag concordance not above the shuffled-flag null = the molecular systems do not track the
  curated clinical profile.
* Presence/absence disagreements (one source has genes, another none) are listed per condition.

## 11. Language

Product language only: candidate, evidence-supported, measurable phenotype, condition-level
molecular enrichment. Never "proves", "confirms mechanism", "definitive biomarker", "patient has".

## Deviations from this plan

Recorded 2026-09-23 after the analysis ran. The pre-specified definitions above were not changed;
everything below is additional and labelled as such in the outputs.

1. **Post-hoc diagnostics (not pre-specified).** After the primary run, two features of the evidence
   were seen: (a) the only GWAS-supported system (gastroparesis, immune) rests on three HLA class II
   genes from one MHC locus of one study; (b) Open Targets `ot_other` support is often driven by the
   clinical (drug-target) datatype, where one trialled mechanism lists a whole gene family (e.g. 51
   complex I subunits for a complex I inhibitor in ME/CFS; 26 calcium-channel genes for a
   calcium-channel modulator). Added: `n_independent_gwas_loci` (GWAS positions merged within 1 Mb) and
   `share_overlap_drug_target_only`, both computed on the genes of a system's enriched pathways
   (`test3_system_enrichment`); `post_hoc_flags` and an exploratory
   `support_status_without_flagged_sources` in `measurable_biology` (flag rule: single GWAS locus;
   >= 50 % drug-target-only genes). The pre-specified `support_status` is unchanged.
2. **Post-hoc sensitivity S8**: `ot_other` without the clinical (drug-target) datatype.
3. **Exploratory flag-concordance variant** `non_literature_sources` (observed systems supported by a
   non-literature source) added next to the two pre-specified definitions.
4. Reporting rules added when writing the report: a permutation p < 0.05 whose bootstrap CI of the
   mean difference includes 0 is reported as "above by the permutation test, but the CI includes 0";
   the mapMECFS comparison is reported as the rank of me_cfs (no threshold).
5. Study-catalogue rows of `measurable_biology` also count SRA RNA studies naming the condition
   (transcriptomics) and map mapMECFS transposable-element rows to transcriptomics and lipid-class /
   stool-metabolite rows to metabolomics.
6. `results/tables/test3_reactome_ora.csv` stores only pathways with overlap > 0 (all tested pathways
   were used for FDR and Spearman). Full source-version strings are in
   `results/tables/test3_run_metadata.json`; per-row `source_version` is abbreviated.
7. The mapMECFS partition was given the common evidence columns before the analysis (task
   instruction); the union was rebuilt with `store.union_partitions`.

### Changes made in independent review (2026-09-23, after results; none alters a pre-specified statistic)

8. **Correction: per-condition disagreement criterion.** `test3_disagreements` labelled a pair
   `overlap_not_condition_specific` when its per-condition empirical p (section 4) was >= 0.05. With
   n_cross = 16-18 cross pairs the smallest attainable p is 1/(1 + n_cross) = 0.053-0.059. That
   threshold could never be met, so 10 pairs whose same-condition fold is above every cross pair were
   listed as disagreements. The criterion is now p <= max(0.05, 1/(1 + n_cross)), and the table carries
   `n_cross_pairs` and `min_attainable_empirical_p`.
9. **Flag concordance, reporting.** The two pre-specified definitions (section 6) are Holm-adjusted
   against each other (`perm_p_holm_prespecified`), and each concordant pair lists its supporting
   sources (`concordant_pairs`). Two exploratory definitions were added, labelled "added in review":
   `condition_specific_non_literature` and `condition_specific_non_literature_unflagged`. Rationale:
   Europe PMC co-mention and the targets of trialled drugs both reflect what is already believed
   clinically about a condition, and the curated flags encode the same beliefs. The flag-concordance p
   of every definition is also written per sensitivity analysis (`test3_sensitivity`, level
   `system_flag_concordance`).
10. **Drug records named.** `test3_system_enrichment.top_drug_mechanisms_behind_overlap` (Open Targets
    'other' rows) names the drug records (ChEMBL mechanism) whose targets account for most of the
    genes behind a system's enriched pathways. `measurable_biology.caveats` repeats it for flagged rows.
11. **Figure 5 drafts and data.** `fig5_edges.post_hoc_flag` marks condition->gene, condition->pathway
    and condition->system edges that rest on flagged evidence (same rule as deviation 1). The draft
    figures draw these edges and cells differently.
12. **Query tool.** `get_measurable_biology` gains a `system_level` block (UNKNOWN when no
    condition->system link can be tested, e.g. ptlds, which otherwise returned only study-catalogue
    rows as "available"). It returns status `not_supported` rather than UNKNOWN when molecular data
    exist, were tested, and support no system mapped to the requested measurement.
13. `test3_gene_sets.csv` rows are sorted, so the file is byte-stable across runs (set iteration order
    is randomised per process).
