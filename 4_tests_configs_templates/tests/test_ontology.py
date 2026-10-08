"""Tests for the ontology layer: measure_it.ontology.build_registry and measure_it.ontology.normalize.

Unit tests use synthetic inputs only. Tests marked `data` read the processed outputs
(run `uv run python -m measure_it.ontology.build_registry` first).
"""
from __future__ import annotations

import io
import json
import zipfile

import pandas as pd
import pytest
import yaml

from measure_it.config import RAW, UNKNOWN, load_config
from measure_it.ontology import build_registry as br
from measure_it.ontology import normalize as nz
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.registry import REQUIRED_FIELDS

# --------------------------------------------------------------------------- text / id helpers


def test_text_keys():
    assert nz.compact_key("Post-exertional malaise") == nz.compact_key("Postexertional malaise")
    assert nz.norm_text("ME/CFS") == "me cfs"
    assert nz.singular_key("Ehlers-Danlos syndromes") == nz.singular_key("Ehlers Danlos syndrome")
    assert nz.singular_key("illness") == "illness"  # 'ss' endings are kept


@pytest.mark.parametrize("text,expected", [
    ("PASC", True), ("CFS", True), ("hEDS", True), ("ME/CFS", True), ("POTS", True),
    ("migraine", False), ("Lyme disease", False), ("Fibromyalgia", False),
])
def test_is_acronym(text, expected):
    assert nz.is_acronym(text) is expected


@pytest.mark.parametrize("raw,expected", [
    ("mesh:D015673", "MESH:D015673"), ("Orphanet:1983", "ORPHANET:1983"), ("SCTID:51771007", "SNOMEDCT:51771007"),
    ("MONDO_0005404", "MONDO:0005404"), ("icd10cm:Q79.6", "ICD10CM:Q79.6"), ("ICD10CM:G9332", "ICD10CM:G93.32"),
    ("EFO_0000555", "EFO:0000555"),
])
def test_canonical_curie(raw, expected):
    assert nz.canonical_curie(raw) == expected


@pytest.mark.parametrize("raw,expected", [
    ("G9332", "G93.32"), ("g93.32", "G93.32"), ("U099", "U09.9"), ("G90A", "G90.A"), ("N80", "N80"),
    ("B90-B94", "B90-B94"), ("not a code", None),
])
def test_icd_dotted(raw, expected):
    assert nz.icd_dotted(raw) == expected


def test_compose_relations():
    assert br.compose("exact", "broad") == "broad"
    assert br.compose("broad", "exact") == "broad"
    assert br.compose("narrow", "narrow") == "narrow"
    assert br.compose("close", "narrow") == "narrow"
    assert br.compose("broad", "narrow") == "related"


def test_mondo_xref_relation():
    assert br.mondo_xref_relation(["MONDO:equivalentTo"], None) == ("exact", "current")
    assert br.mondo_xref_relation(["MONDO:equivalentObsolete"], None)[1] == "obsolete"
    assert br.mondo_xref_relation(["MONDO:relatedTo"], None) == ("related", "current")
    assert br.mondo_xref_relation(["MONDO:mondoIsNarrowerThanSource"], None) == ("broad", "current")
    # skos predicate wins over the qualifier when present
    assert br.mondo_xref_relation(["MONDO:equivalentTo"], "skos:closeMatch") == ("close", "current")


def test_description_match():
    assert "contains name" in br.description_match("Postural orthostatic tachycardia syndrome [POTS]",
                                                   ["Postural orthostatic tachycardia syndrome"])
    assert "contains CMS description" in br.description_match("Hypermobility syndrome", ["joint hypermobility syndrome"])
    assert br.description_match("Mast cell activation, unspecified", ["mast cell activation syndrome"])
    assert br.description_match("Secondary mast cell activation", ["mast cell activation syndrome"]) == ""
    # most specific (longest) name is reported
    assert "hypermobile" in br.description_match("Hypermobile Ehlers-Danlos syndrome",
                                                 ["Ehlers-Danlos syndromes", "hypermobile Ehlers-Danlos syndrome"])


def test_repair_flow_split():
    # YAML `[endometriosis of uterus, adenomyosis (related, distinct)]` parses into three items
    parsed = yaml.safe_load("a: [endometriosis of uterus, adenomyosis (related, distinct)]")["a"]
    assert parsed == ["endometriosis of uterus", "adenomyosis (related", "distinct)"]
    assert br.repair_flow_split(parsed) == ["endometriosis of uterus", "adenomyosis (related, distinct)"]
    assert br.repair_flow_split(["PTLDS", "chronic Lyme disease (contested term)"]) == \
        ["PTLDS", "chronic Lyme disease (contested term)"]
    assert br.repair_flow_split([]) == []
    assert br.repair_flow_split(["x (never closed"]) == ["x (never closed"]


def test_owl_version_iri(tmp_path):
    p = tmp_path / "h.xml"
    p.write_text('<owl:Ontology rdf:about="http://purl.obolibrary.org/obo/mondo.owl">\n'
                 '  <owl:versionIRI rdf:resource="http://purl.obolibrary.org/obo/mondo/releases/2026-09-01/mondo.owl"/>')
    assert br.owl_version_iri(p) == "http://purl.obolibrary.org/obo/mondo/releases/2026-09-01/mondo.owl"


def test_mondo_icd10cm_xref_index():
    terms = {"MONDO:1": {"name": "a", "obsolete": False, "xref": [("ICD10CM:D89.41", [])]},
             "MONDO:2": {"name": "b", "obsolete": True, "xref": [("ICD10CM:D89.43", [])]}}
    rev = br.mondo_icd10cm_xref_index(terms)
    assert rev["ICD10CM:D89.41"] == [("MONDO:1", "a")]
    assert "ICD10CM:D89.43" not in rev  # obsolete classes do not count


def test_name_variants_split_preferred_name():
    c = {"id": "x", "preferred_name": "Dysautonomia (autonomic nervous system disorder)",
         "aliases": ["POTS", "adenomyosis (related, distinct)"], "search_terms": ["dysautonomia"]}
    v = br.name_variants(c)
    texts = [x["text"] for x in v]
    assert texts[0] == "Dysautonomia"
    assert "autonomic nervous system disorder" in texts
    assert "adenomyosis" in texts and "adenomyosis (related, distinct)" not in texts
    assert [x for x in v if x["text"] == "POTS"][0]["acronym"] is True
    assert texts.count("Dysautonomia") == 1  # the search term is a duplicate


OBO = """format-version: 1.2
data-version: releases/2026-09-01

[Term]
id: MONDO:1
name: partial androgen insensitivity syndrome
synonym: "PAIS" EXACT []

[Term]
id: MONDO:2
name: myalgic encephalomeyelitis/chronic fatigue syndrome
synonym: "chronic fatigue syndrome" EXACT [DOID:8544]
synonym: "CFS" EXACT ABBREVIATION []
xref: EFO:0004540 {source="MONDO:equivalentTo", source="MONDO:EFO"}
xref: Orphanet:1983 {source="MONDO:equivalentObsolete"}
property_value: skos:exactMatch EFO:0004540

[Term]
id: MONDO:3
name: obsolete thing
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


def test_parse_obo_and_match(tmp_path):
    p = tmp_path / "t.obo"
    p.write_text(OBO)
    header, terms = br.parse_obo(p)
    assert header["data-version"] == "releases/2026-09-01"
    assert set(terms) == {"MONDO:1", "MONDO:2", "MONDO:3"}
    assert terms["MONDO:3"]["obsolete"] is True
    xr = {x["target_id"]: x for x in br.mondo_xrefs(terms["MONDO:2"])}
    assert xr["EFO:0004540"]["class_relation"] == "exact"
    assert xr["ORPHANET:1983"]["target_status"].startswith("obsolete")

    idx = br.build_mondo_index(terms)
    c = {"id": "x", "preferred_name": "Post-infectious syndromes (grouping)", "aliases": ["PAIS", "CFS"],
         "search_terms": ["chronic fatigue syndrome"]}
    cands = br.match_mondo(br.name_variants(c), idx, terms)
    by = {(x["text"], x["mondo_id"]): x for x in cands}
    assert by[("PAIS", "MONDO:1")]["acronym"] is True          # homonym hit is flagged as acronym-only
    assert by[("chronic fatigue syndrome", "MONDO:2")]["match_type"] == "exact synonym"
    assert not any(x["mondo_id"] == "MONDO:3" for x in cands)  # obsolete classes never match


def test_parse_cms_order(tmp_path):
    lines = [
        "00001 G90     0 Disorders of autonomic nervous system                        Disorders of autonomic nervous system",
        "00002 G90A    1 Postural orthostatic tachycardia syndrome [POTS]             Postural orthostatic tachycardia syndrome [POTS]",
    ]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("Code Descriptions/icd10cm_order_2026.txt", "\r\n".join(lines) + "\r\n")
    z = tmp_path / "c.zip"
    z.write_bytes(buf.getvalue())
    df = br.parse_cms_order(z)
    assert df["code"].tolist() == ["G90", "G90.A"]
    assert df["billable"].tolist() == [False, True]
    assert df.loc[1, "description"] == "Postural orthostatic tachycardia syndrome [POTS]"


# --------------------------------------------------------------------------- normalize (synthetic index)


def _toy_index() -> nz.OntologyIndex:
    reg = pd.DataFrame([
        dict(canonical_condition_id="long_covid", preferred_name="Long COVID (post-COVID-19 condition)",
             aliases=["PASC", "post-COVID-19 condition"], search_terms=["long COVID"], primary_mondo_id="MONDO:0100233",
             mondo_ids=["MONDO:0100233"], icd10cm_codes=["U09.9"], efo_ids=[], mesh_ids=[], primary_mondo_relation="exact",
             mondo_labels=["long COVID-19"], mondo_exact_synonyms=["PASC"]),
        dict(canonical_condition_id="me_cfs", preferred_name="Myalgic encephalomyelitis/chronic fatigue syndrome",
             aliases=["ME/CFS", "CFS", "chronic fatigue syndrome"], search_terms=[], primary_mondo_id="MONDO:0005404",
             mondo_ids=["MONDO:0005404"], icd10cm_codes=["G93.32"], efo_ids=[], mesh_ids=["MESH:D015673"],
             primary_mondo_relation="exact", mondo_labels=[], mondo_exact_synonyms=[]),
        dict(canonical_condition_id="pots", preferred_name="Postural orthostatic tachycardia syndrome",
             aliases=["POTS"], search_terms=[], primary_mondo_id="MONDO:0001315", mondo_ids=["MONDO:0001315"],
             icd10cm_codes=["G90.A"], efo_ids=[], mesh_ids=[], primary_mondo_relation="broad",
             mondo_labels=["orthostatic intolerance"], mondo_exact_synonyms=[]),
        dict(canonical_condition_id="eds_hsd", preferred_name="Ehlers-Danlos syndromes / hypermobility spectrum disorders",
             aliases=["hEDS", "hypermobile Ehlers-Danlos syndrome"], search_terms=[], primary_mondo_id="MONDO:0020066",
             mondo_ids=["MONDO:0020066"], icd10cm_codes=["Q79.62"], efo_ids=[], mesh_ids=[], primary_mondo_relation="narrow",
             mondo_labels=[], mondo_exact_synonyms=[]),
        dict(canonical_condition_id="endometriosis", preferred_name="Endometriosis", aliases=[], search_terms=[],
             primary_mondo_id="MONDO:0005133", mondo_ids=["MONDO:0005133"], icd10cm_codes=["N80"], efo_ids=[],
             mesh_ids=[], primary_mondo_relation="exact", mondo_labels=[], mondo_exact_synonyms=[]),
    ])

    def m(cid, ont, tid, pred, role=None, status="current", kind=None, src="test"):
        return dict(canonical_condition_id=cid, target_ontology=ont, target_id=tid, predicate_condition=pred,
                    mondo_role=role, target_status=status, code_kind=kind, in_registry_list=status == "current",
                    mapping_source=src, target_label=None)

    maps = pd.DataFrame([
        m("long_covid", "ICD10CM", "ICD10CM:U09.9", "exact", kind="code"),
        m("me_cfs", "ICD10CM", "ICD10CM:G93.32", "exact", kind="code"),
        m("pots", "ICD10CM", "ICD10CM:G90.A", "exact", kind="code"),
        m("eds_hsd", "ICD10CM", "ICD10CM:Q79.62", "narrow", kind="code"),
        m("endometriosis", "ICD10CM", "ICD10CM:N80", "exact", kind="code"),
        m("me_cfs", "MONDO", "MONDO:0005404", "exact", role="primary"),
        m("pots", "MONDO", "MONDO:0001315", "broad", role="primary"),
        m("me_cfs", "EFO", "EFO:0004540", "exact", role="primary", status="obsolete in OLS4 (replaced by MONDO:0005404)"),
        m("me_cfs", "MESH", "MESH:D015673", "exact", role="primary"),
        {**m("pots", "MESH", "MESH:D054972", "narrow", role="narrower"), "target_label": "Neurocirculatory Toy Label"},
    ])
    phen = pd.DataFrame([dict(canonical_condition_id="eds_hsd", hpo_id="HP:0001382", hpo_label="Joint hypermobility",
                              association_scope="direct", negated=False)])
    axes = pd.DataFrame([dict(axis_id="fatigue", axis_label="Fatigue", hpo_id="HP:0012378", hpo_label="Fatigue")])
    codes = pd.DataFrame([dict(code_nodot="Z0000", code="Z00.00", description="Encounter for general adult medical "
                               "examination without abnormal findings", billable=True, in_fy2026=True)])
    return nz.OntologyIndex.from_tables(reg, maps, phen, axes, codes)


@pytest.fixture(scope="module")
def toy():
    return _toy_index()


@pytest.mark.parametrize("query,cid", [
    ("PASC", "long_covid"), ("U09.9", "long_covid"), ("u099", "long_covid"), ("CFS", "me_cfs"),
    ("G90.A", "pots"), ("hEDS", "eds_hsd"), ("G93.32", "me_cfs"), ("N80.101", "endometriosis"),
    ("MONDO:0005404", "me_cfs"), ("MONDO_0005404", "me_cfs"), ("MESH:D015673", "me_cfs"), ("D015673", "me_cfs"),
    ("post-COVID-19 condition", "long_covid"), ("Long Covid", "long_covid"), ("POTS", "pots"),
    ("HP:0001382", "eds_hsd"),
])
def test_normalize_toy(toy, query, cid):
    res = nz.normalize_condition(query, index=toy)
    assert res["status"] != UNKNOWN, res
    assert res["matches"][0]["canonical_condition_id"] == cid, res
    assert res["matches"][0]["match_reason"]


def test_normalize_unknowns(toy):
    for q in ["banana", "", "Z00.00", "HP:0000001", "MONDO:9999999"]:
        res = nz.normalize_condition(q, index=toy)
        assert res["status"] == UNKNOWN, (q, res)
        assert res["reason"]
    assert "valid FY2026 code" in nz.normalize_condition("Z00.00", index=toy)["reason"]
    fat = nz.normalize_condition("fatigue", index=toy)
    assert fat["status"] == UNKNOWN and fat["phenotype_axis"]["hpo_id"] == "HP:0012378"


def test_normalize_obsolete_efo_is_flagged(toy):
    res = nz.normalize_condition("EFO:0004540", index=toy)
    assert res["matches"][0]["canonical_condition_id"] == "me_cfs"
    assert "obsolete" in res["matches"][0]["match_reason"]


def test_broad_primary_mondo_class_is_not_a_match(toy):
    # pots' primary class MONDO:0001315 is curated as BROADER than POTS: its role must not override the predicate
    res = nz.normalize_condition("MONDO:0001315", index=toy)
    assert res["status"] == "partial" and res["matches"][0]["canonical_condition_id"] == "pots"
    assert res["matches"][0]["confidence"] == nz.PRED_CONF["broad"]
    assert nz.normalize_condition("MONDO:0005404", index=toy)["status"] == "matched"


def test_narrow_id_label_is_partial(toy):
    # the label of a narrower MeSH id is a subtype name: an exact text hit is partial, not matched
    res = nz.normalize_condition("neurocirculatory toy label", index=toy)
    assert res["status"] == "partial" and res["matches"][0]["canonical_condition_id"] == "pots"


def test_icd_exact_beats_category(toy):
    res = nz.normalize_condition("N80", index=toy)
    assert res["status"] == "matched" and res["matches"][0]["confidence"] == 1.0
    res = nz.normalize_condition("N80.101", index=toy)
    assert "falls under mapped ICD-10-CM category N80" in res["matches"][0]["match_reason"]


def test_search_condition_acronyms_are_exact_only(toy):
    assert nz.search_condition("POT", index=toy) == [] or nz.search_condition("POT", index=toy)[0]["canonical_condition_id"] != "pots"
    top = nz.search_condition("long covid brain fog", index=toy)[0]
    assert top["canonical_condition_id"] == "long_covid" and "query contains" in top["match_reason"]


# --------------------------------------------------------------------------- processed outputs


@pytest.fixture(scope="module")
def tables():
    from measure_it.store import read_table, table_exists
    names = ["condition_registry", "condition_ontology_mappings", "condition_phenotypes", "phenotype_axes",
             "ontology_icd10cm_codes"]
    if not all(table_exists(n) for n in names):
        pytest.skip("ontology tables not built; run uv run python -m measure_it.ontology.build_registry")
    return {n: read_table(n) for n in names}


@pytest.mark.data
def test_registry_covers_config(tables):
    reg = tables["condition_registry"]
    cfg_ids = [c["id"] for c in load_config("conditions")["conditions"]]
    assert sorted(reg["canonical_condition_id"]) == sorted(cfg_ids)
    assert reg["canonical_condition_id"].is_unique
    for t in tables.values():
        assert set(PROVENANCE_COLUMNS) <= set(t.columns)
        assert (t["data_layer"] == "ontology").all()
    assert reg["mondo_ids"].map(len).gt(0).all()
    for col in br.FLAG_COLUMNS:
        assert reg[col].dtype == bool
    prov = json.loads(reg.iloc[0]["ontology_provenance"])
    assert prov["versions"]["mondo"].startswith("Mondo v")


@pytest.mark.data
def test_key_resolutions(tables):
    reg = tables["condition_registry"].set_index("canonical_condition_id")
    assert reg.at["me_cfs", "primary_mondo_id"] == "MONDO:0005404"
    assert reg.at["me_cfs", "primary_mondo_match_type"] == "exact synonym"
    assert reg.at["long_covid", "primary_mondo_id"] == "MONDO:0100233"
    assert reg.at["pots", "primary_mondo_relation"] == "broad"          # POTS only exists as a synonym
    assert reg.at["post_infectious_syndrome", "primary_mondo_match_type"] == "manual choice"
    assert "U09.9" in reg.at["long_covid", "icd10cm_codes"]
    assert "G93.32" in reg.at["me_cfs", "icd10cm_codes"]
    assert "G90.A" in reg.at["pots", "icd10cm_codes"]
    assert "MESH:D054972" in reg.at["pots", "mesh_ids"]
    assert json.loads(reg.at["pots", "id_predicates"])["MESH:D054972"] == "exact"
    assert "MONDO_0005404" in reg.at["me_cfs", "efo_short_forms"]


@pytest.mark.data
def test_icd_codes_validated_against_cms(tables):
    reg, maps, codes = tables["condition_registry"], tables["condition_ontology_mappings"], tables["ontology_icd10cm_codes"]
    fy26 = set(codes.loc[codes["in_fy2026"], "code"])
    listed = {c for lst in reg["icd10cm_codes"] for c in lst}
    assert listed and listed <= fy26
    icd = maps[maps["target_ontology"] == "ICD10CM"]
    curated = icd[icd["mapping_source"] == br.CMS_CURATED_SOURCE]
    assert len(curated) > 0
    assert curated["description_match"].str.len().gt(0).all()
    assert (curated["in_cms_fy2026"] == True).all()  # noqa: E712
    assert (curated["target_label"] == curated["cms_description"]).all()


@pytest.mark.data
def test_mapping_integrity(tables):
    maps, reg = tables["condition_ontology_mappings"], tables["condition_registry"]
    assert set(maps["predicate_condition"]) <= set(nz.PREDICATE_RANK)
    assert maps["mapping_source"].notna().all()
    assert maps["mapping_id"].is_unique
    # obsolete / not-found targets never reach the registry id lists
    bad = set(maps.loc[maps["target_status"] != "current", "target_id"])
    for col in ["efo_ids", "mesh_ids", "doid_ids", "orphanet_ids"]:
        for lst in reg[col]:
            assert not (set(lst) & bad)
    # every EFO id in the registry was confirmed present (and current) in EFO via OLS
    efo_ok = set(maps.loc[(maps["target_ontology"] == "EFO") & (maps["target_status"] == "current"), "target_id"])
    assert {e for lst in reg["efo_ids"] for e in lst} <= efo_ok


@pytest.mark.data
def test_phenotypes_and_axes(tables):
    ph, ax = tables["condition_phenotypes"], tables["phenotype_axes"]
    assert ph["hpo_id"].str.startswith("HP:").all()
    assert set(ph["association_scope"]) <= {"direct", "narrower_class", "descendant"}
    has_freq = ph["frequency_qualifier"].notna()
    assert ph.loc[has_freq, "frequency_label"].notna().all()
    assert ax["axis_id"].is_unique
    matched = ax[ax["hpo_id"].notna()]
    assert matched["hpo_id"].str.startswith("HP:").all() and matched["hpo_label"].notna().all()
    oi = ax.set_index("axis_id").loc["orthostatic_intolerance"]
    assert pd.isna(oi["hpo_id"]) and len(oi["closest_hpo_terms_if_unmatched"]) > 0


@pytest.mark.data
@pytest.mark.parametrize("query,cid", [
    ("PASC", "long_covid"), ("U09.9", "long_covid"), ("CFS", "me_cfs"), ("G90.A", "pots"), ("hEDS", "eds_hsd"),
    ("POTS", "pots"), ("G93.32", "me_cfs"), ("EFO:0000555", "ibs"),
    ("D015673", "me_cfs"), ("fibromyalgia syndrome", "fibromyalgia"), ("irritable bowel syndrome", "ibs"),
    ("N80.101", "endometriosis"), ("M79.7", "fibromyalgia"), ("post-Lyme disease syndrome", "ptlds"),
    ("D89.41", "mcas"), ("D89.43", "mcas"), ("M35.7", "eds_hsd"), ("hypermobility syndrome", "eds_hsd"),
])
def test_normalize_real(tables, query, cid):
    nz.get_index.cache_clear()
    res = nz.normalize_condition(query)
    assert res["status"] in ("matched", "ambiguous"), res
    assert res["matches"][0]["canonical_condition_id"] == cid, res


@pytest.mark.data
def test_normalize_real_broad_mondo_class(tables):
    # MONDO:0001315 'orthostatic intolerance' is POTS' primary class but curated as broader than POTS, and it is a
    # narrower class of dysautonomia: the result must be partial with both candidates, never a confident POTS match.
    nz.get_index.cache_clear()
    res = nz.normalize_condition("MONDO:0001315")
    assert res["status"] == "partial", res
    cands = {m["canonical_condition_id"] for m in res["matches"] + res["other_candidates"]}
    assert {"pots", "dysautonomia"} <= cands


@pytest.mark.data
def test_normalize_real_unknown(tables):
    for q in ["banana", "Z00.00", "fatigue", "distinct", "distinct)"]:
        assert nz.normalize_condition(q)["status"] == UNKNOWN, q


@pytest.mark.data
def test_registry_aliases_not_split(tables):
    reg = tables["condition_registry"]
    for col in ("aliases", "search_terms"):
        for lst in reg[col]:
            for a in lst:
                assert a.count("(") == a.count(")"), a
    endo = reg.set_index("canonical_condition_id").at["endometriosis", "aliases"]
    assert "adenomyosis (related, distinct)" in list(endo)


@pytest.mark.data
def test_curated_cms_codes_have_no_mondo_xref(tables):
    """CMS-description additions are allowed only where Mondo has no xref for the code (checked against Mondo SSSOM)."""
    maps = tables["condition_ontology_mappings"]
    sssom = pd.read_csv(RAW / br.SOURCE_MONDO / "mondo.sssom.tsv", sep="\t", comment="#", dtype=str)
    mondo_icd = {nz.canonical_curie(o) for o in sssom["object_id"] if str(o).upper().startswith("ICD10CM:")}
    curated = set(maps.loc[maps["mapping_source"] == br.CMS_CURATED_SOURCE, "target_id"])
    assert curated and not (curated & mondo_icd), curated & mondo_icd
    # the MCAS subclass codes and M35.7 arrive through Mondo xrefs
    icd = maps[maps["target_ontology"] == "ICD10CM"].set_index(["canonical_condition_id", "target_id"])
    for key, mondo_id in {("mcas", "ICD10CM:D89.41"): "MONDO:0033954", ("mcas", "ICD10CM:D89.43"): "MONDO:0100006",
                          ("eds_hsd", "ICD10CM:M35.7"): "MONDO:0001798"}.items():
        assert icd.loc[key, "mondo_id"] == mondo_id and "xref" in icd.loc[key, "mapping_source"]


@pytest.mark.data
def test_related_distinct_ids_stay_out_of_lists(tables):
    """DOID:288 'endometriosis of uterus' (DOID synonyms: adenomyosis) is Mondo's exact xref of the related-but-distinct
    adenomyosis class, so it must not reach endometriosis' id lists through the alias name search."""
    reg, maps = tables["condition_registry"].set_index("canonical_condition_id"), tables["condition_ontology_mappings"]
    assert "DOID:288" not in list(reg.at["endometriosis", "doid_ids"])
    rel = maps[maps["mondo_role"] == "related_distinct"]
    for cid, tid in rel.loc[rel["class_relation"] == "exact", ["canonical_condition_id", "target_id"]].itertuples(index=False):
        lists = [x for col in ("efo_ids", "mesh_ids", "doid_ids", "umls_ids", "snomedct_ids", "orphanet_ids", "ncit_ids",
                               "icd10cm_codes") for x in reg.at[cid, col]]
        assert tid not in lists and tid.split(":", 1)[1] not in lists, (cid, tid)


@pytest.mark.data
def test_mondo_version_iri_recorded(tables):
    prov = json.loads(tables["condition_registry"].iloc[0]["ontology_provenance"])
    assert prov["versions"]["mondo_version_iri"] == "http://purl.obolibrary.org/obo/mondo/releases/2026-09-01/mondo.owl"
    assert (RAW / br.SOURCE_MONDO / br.MONDO_OWL_HEADER_FILE).read_text().count("owl:versionIRI") >= 1


@pytest.mark.data
def test_audit_counts_match_tables(tables):
    ax = tables["phenotype_axes"]
    n_missing = int(ax["hpo_id"].isna().sum())
    assert f"Axes without an HPO term: {n_missing} " in (RAW / br.SOURCE_OLS / "DATA_AUDIT.md").read_text()
    entry = yaml.safe_load((RAW / br.SOURCE_MONDO / "registry_entry.yaml").read_text())
    assert entry["sample_size"]["mapping_rows"] == len(tables["condition_ontology_mappings"])
    assert entry["sample_size"]["phenotype_rows"] == len(tables["condition_phenotypes"])
    assert entry["sample_size"]["phenotype_distinct_hpo"] == tables["condition_phenotypes"]["hpo_id"].nunique()
    cms = yaml.safe_load((RAW / br.SOURCE_CMS / "registry_entry.yaml").read_text())
    codes = tables["ontology_icd10cm_codes"]
    assert cms["sample_size"]["fy2026_rows"] == int(codes["in_fy2026"].sum())
    assert cms["sample_size"]["fy2027_rows"] == int(codes["in_fy2027"].sum())


@pytest.mark.data
def test_registry_entries_and_audits():
    for sid in (br.SOURCE_MONDO, br.SOURCE_OLS, br.SOURCE_CMS):
        d = RAW / sid
        entry = yaml.safe_load((d / "registry_entry.yaml").read_text())
        assert not [f for f in REQUIRED_FIELDS if f not in entry]
        audit = (d / "DATA_AUDIT.md").read_text()
        assert "| source_id |" in audit and "## Missingness" in audit
        assert (d / "MANIFEST.json").exists()


@pytest.mark.data
def test_every_mondo_candidate_reviewed(tables):
    """A full-name Mondo match must be used or explicitly curated; manual choices must name what they override."""
    cands = pd.read_csv(RAW / br.SOURCE_MONDO / "condition_mondo_candidates.tsv", sep="\t")
    unreviewed = cands[cands["status"].str.startswith("not used")]
    assert unreviewed.empty, unreviewed[["canonical_condition_id", "name", "mondo_id", "mondo_label"]].to_dict("records")
    assert set(cands.loc[cands["status"].str.startswith("rejected"), "acronym"]) <= {True}
    maps = tables["condition_ontology_mappings"]
    manual = maps[(maps["target_ontology"] == "MONDO") & (maps["mondo_role"] == "primary")
                  & (maps["match_type"] == "manual choice")]
    for r in manual.itertuples():
        auto = cands[(cands["canonical_condition_id"] == r.canonical_condition_id) & ~cands["acronym"]]
        if len(auto):
            assert "overrides automatic best match" in r.rationale
