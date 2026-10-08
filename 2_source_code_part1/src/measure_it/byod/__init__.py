"""Bring your own data (BYOD): add an outside group's person-level data to the engine without editing core code.

    uv run measure-it byod validate <dir>     # schema + privacy checks; refuses identifying / geographic content
    uv run measure-it byod ingest <dir>       # namespaced person-layer partitions, provenance, registry entry, audit
    uv run measure-it byod evaluate <dir>     # locks the manifest's primary analysis, then CV / bootstrap / permutation
    uv run measure-it byod subgroup <dir>     # locks the subgroup plan; device-defined subgroup -> computable phenotype
    uv run measure-it byod deploy <dir>       # metric link + evidence_weighted scoring with the new record
    uv run measure-it byod remove <id>        # removes the dataset and every record it produced
    uv run measure-it byod list               # user datasets on this machine
    uv run measure-it byod show <id>          # locked plan, models, performance record, deploy report (JSON)

Contract: docs/BRING_YOUR_OWN_DATA.md; template folder templates/byod/; worked example docs/WALKTHROUGH.md.

What crosses from the person layer into the ranking is exactly one aggregate performance record per dataset (AUROC and
CI, sensitivity at the pre-specified specificity, n, label basis, comparator), in the table `byod_performance_records`
that `measure_it.scoring.metric_link` reads. Person rows never carry a geographic column (refused at validation,
re-checked after ingestion) and are never joined to a place, facility or another dataset's participants.

Modules
    common      paths, constants, manifest loading, file discovery
    checks      validation (schema, privacy / geography refusals, group sizes) -> report
    adapters    TabularFeatureAdapter (pre-computed features) and loading a partner's own MeasurementAdapter
    ingest      person-layer partitions, registry entry, DATA_AUDIT.md, canonical unions, digital person
    person      the dataset's Digital Phenotype Vector spec (read by measure_it.wearables.digital_person)
    evaluate    plan lock + the pre-specified analyses -> phenotype_signatures__byod_<id>, performance record
    subgroup    device-defined subgroup inside a label's cases -> computable phenotype (measure_it.harmonize)
    records     performance records for measure_it.scoring.metric_link (demo records excluded by default)
    deploy      metric link, re-scoring of the affected combinations, ranking-change report, heatmap / embedding
    remove      clean removal and restoration of the rankings
    tools       list_user_datasets / get_user_dataset_metric (tool-facade envelopes)
"""
