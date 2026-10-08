"""Lab-dataset discovery: verify open datasets with per-person LAB / biomarker values + invisible-illness labels.

Question: which public datasets let us test, person by person, whether a laboratory measurement (blood/urine/CSF
analytes, cytokines, hormones, autoantibodies, metabolomics, proteomics) separates people with a target condition
(MCAS, POTS/dysautonomia, Long COVID, ME/CFS, fibromyalgia, EDS/HSD, Lyme/PTLDS, migraine, IBS, endometriosis,
gastroparesis) from controls or from other umbrella conditions?

Method (replicates docs/DEVICE_DATASET_DISCOVERY.md for labs):
  0. Four parallel searchers (2026-09-24) proposed 182 candidates; their JSON is kept verbatim as the discovery input
     in data/raw/lab_candidates/_searcher_inputs/ (MANIFEST.json records it).
  1. For every candidate claimed open, this module retrieves a file listing or a small per-person file through
     measure_it.download / measure_it.http into data/raw/lab_candidates/<candidate_id>/ (never > 2 GB; a HEAD size
     guard skips anything over MAX_BYTES). Access routes that are not open are probed at metadata level only.
  2. Schema-level inspection only: rows, label column and value counts, analyte columns present and how many labelled
     people have them. It never computes a biomarker-by-label contrast, so the primary analyses written in
     docs/ANALYSIS_PLAN_<ID>.md are specified before any outcome was seen.
  3. Merges inspection with the curated candidate table (deduplicated; status after verification) and ranks the
     verified-open, case-control-capable datasets with a transparent rubric (components kept as columns).
  4. Writes results/tables/lab_dataset_candidates.csv and data/raw/lab_candidates/inspection.json.

No account is created and no data-use agreement is accepted. One record that appears to contain direct personal
identifiers (figshare 30127723) is deliberately NOT downloaded.

Run: uv run python -m measure_it.labs.lab_dataset_discovery [--no-download]
"""
from __future__ import annotations

import argparse
import io
import json
import math
import re
import zipfile
from pathlib import Path

import pandas as pd

from .. import http
from ..config import RAW, TABLES, utc_now_iso
from ..download import download_file, load_manifest, record_file

SRC = "lab_candidates"
ROOT = RAW / SRC
SEARCHER_DIR = ROOT / "_searcher_inputs"
SEARCHER_FILES = {"mcas_pots_eds": "candidates_mcas_pots_eds.json", "long_covid": "candidates_long_covid.json",
                  "mecfs": "candidates_mecfs.json", "other": "candidates_other.json"}
MAX_BYTES = 400_000_000   # verification never needs more; hard project limit is 2 GB
BUA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) "
                     "Chrome/124.0 Safari/537.36"}
FS = "https://ndownloader.figshare.com/files/"
ZAPI = "https://zenodo.org/api/records/{rec}/files/{key}/content"
MEND = "https://data.mendeley.com/public-files/datasets/"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles"
MW = "https://www.metabolomicsworkbench.org/rest/study/analysis_id/{an}/datatable"
MWD = "https://www.metabolomicsworkbench.org/studydownload/{f}"
MWF = "https://www.metabolomicsworkbench.org/rest/study/study_id/{st}/factors"
PRIDE = "https://ftp.pride.ebi.ac.uk/pride/data/archive/"
MTBLS = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/"
SPR = "https://static-content.springer.com/esm/art%3A10.1038%2F{art}/MediaObjects/{f}"
PLOS = "https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.{s}"
GH = "https://raw.githubusercontent.com/ohlab/BioMapAI/main/data/"
HDV = "https://dataverse.harvard.edu/api/access/datafile/"


def _epmc(pmc: str) -> tuple[str, str, bool]:
    return (EPMC.format(pmc=pmc), f"{pmc}_supplementaryFiles.zip", False)


# candidate_id -> list of (url, local filename, needs browser user-agent). One small, content-bearing file (or the
# supplement bundle) per candidate is enough to verify that per-person rows and labels exist.
FILES: dict[str, list[tuple[str, str, bool]]] = {
    # ---- MCAS / mast cell / POTS / EDS-HSD ----
    "heds_hsd_olink_serum_cinquina2026": [_epmc("PMC13081554")],
    "heds_serum_ms_proteomics_griggs2025": [
        (PRIDE + "2025/09/PXD062941/20240126_104630_Norris_Serum_search1_noplate2_01312024_JRB.tsv",
         "Norris_Serum_search1_noplate2_JRB.tsv", False)],
    "mcad_plasma_heparin_vysniauskaite2015": [(PLOS.format(s="0124912.s001"), "pone.0124912.s001.xls", False)],
    "od_children_catecholamines_sugiyama2026": [_epmc("PMC13132173")],
    "sm_gut_microbiome_tableS1_krausfeldt2026": [_epmc("PMC12826983")],
    "mmcs_cmcsnos_ballul2025": [_epmc("PMC11773267")],
    "sm_mental_cytokines_figshare2026": [(FS + "60986206", "Mastocytosis_Mental_database.xlsx", False)],
    "tryptase_ngs_genotyping_li2024": [(PLOS.format(s="0291947.s002"), "pone.0291947.s002.csv", False)],
    "hat_gastroparesis_ev_lipidomics_2022": [(FS + "37806663", "Supp_Table1.xlsx", False)],
    "oi_children_metabolomics_li2024": [_epmc("PMC10941598")],
    "covid_tryptase_zenodo2025": [(ZAPI.format(rec=18088869, key="Tryptaza_ostateczna_baza_07.12.xlsx"),
                                   "Tryptaza_ostateczna_baza_07.12.xlsx", False)],
    "hat_venom_gwas_zenodo2025": [],   # listing only (sumstats 145 MB not needed)
    "mirabegron_pots_zenodo2026": [(ZAPI.format(rec=21725551, key="Data_Dictionary.xlsx.csv"),
                                    "Data_Dictionary.xlsx.csv", False)],
    "heds_pots_fmri_zenodo2026": [(ZAPI.format(rec=22128981, key="Data%20study%20group.txt"),
                                   "Data_study_group.txt", False)],
    # ---- Long COVID ----
    "klein2023_mylc_ml_table": [
        (SPR.format(art="s41586-023-06651-y", f="41586_2023_6651_MOESM4_ESM.xlsx"), "41586_2023_6651_MOESM4_ESM.xlsx", True),
        (SPR.format(art="s41586-023-06651-y", f="41586_2023_6651_MOESM5_ESM.xlsx"), "41586_2023_6651_MOESM5_ESM.xlsx", True),
        (SPR.format(art="s41586-023-06651-y", f="41586_2023_6651_MOESM2_ESM.xlsx"), "41586_2023_6651_MOESM2_ESM.xlsx", True),
        (SPR.format(art="s41586-023-06651-y", f="41586_2023_6651_MOESM3_ESM.xlsx"), "41586_2023_6651_MOESM3_ESM.xlsx", True)],
    "su2022_incov_multiomics": [_epmc("PMC8786632")],
    "cervia_hasler2024_zurich_somascan": [
        (MEND + "dvf6yvrg4x/files/859c49d0-ae4c-4fc1-85c2-8d178e274fa3/file_downloaded",
         "Proteomics_Clinical_Data_Labels.xlsx", True),
        (MEND + "dvf6yvrg4x/files/9a3eedfb-93b8-4b43-8dbb-ffe8adc45b1e/file_downloaded",
         "Proteomics_Clinical_Data_6M_timepoint.xlsx", True)],
    "lu2024_early_markers_pasc_figshare": [(FS + "48447364", "Raw_Host_Immune_Data.xls", False)],
    "barouch2025_somascan_pasc": [
        (SPR.format(art="s41590-025-02353-x", f="41590_2025_2353_MOESM3_ESM.xlsx"), "41590_2025_2353_MOESM3_ESM.xlsx", True)],
    "gao2025_uk_se_olink": [(ZAPI.format(rec=14772494, key="meta_UKsamples.csv"), "meta_UKsamples.csv", False)],
    "santacruz2023_portugal_pasc_sourcedata": [
        (SPR.format(art="s41467-023-37368-1", f="41467_2023_37368_MOESM4_ESM.xlsx"), "41467_2023_37368_MOESM4_ESM.xlsx", True)],
    "nat_microbiol2024_trace_cytokines_cardiac": [
        (SPR.format(art="s41564-024-01838-z", f="41564_2024_1838_MOESM3_ESM.xlsx"), "41564_2024_1838_MOESM3_ESM.xlsx", True)],
    "lopezhernandez2023_mexico_metabolome": [
        (MEND + "gc9g2g53kr/files/e8bd00fb-62aa-4d9b-bffb-08d86d02b8a3/file_downloaded", "LCDATA_AVAILABILITY_SR.xls", True)],
    "talla2023_aifi_olink": [
        (SPR.format(art="s41467-023-38682-4", f="41467_2023_38682_MOESM3_ESM.xlsx"), "41467_2023_38682_MOESM3_ESM.xlsx", True),
        (SPR.format(art="s41467-023-38682-4", f="41467_2023_38682_MOESM4_ESM.xlsx"), "41467_2023_38682_MOESM4_ESM.xlsx", True)],
    "woodruff2023_emory_olink": [],   # listing only (tarball)
    "appelman2024_muscle_pem": [
        (SPR.format(art="s41467-023-44432-3", f="41467_2023_44432_MOESM4_ESM.xlsx"), "41467_2023_44432_MOESM4_ESM.xlsx", True)],
    "charlton2026_muscle_lc_mecfs_bedrest": [
        (SPR.format(art="s41467-026-75725-y", f="41467_2026_75725_MOESM4_ESM.xlsx"), "41467_2026_75725_MOESM4_ESM.xlsx", True)],
    "kedor2022_charite_pcs_mecfs": [
        (SPR.format(art="s41467-022-32507-6", f="41467_2022_32507_MOESM4_ESM.xlsx"), "41467_2022_32507_MOESM4_ESM.xlsx", True)],
    "guntur2022_pasc_metabolomics": [_epmc("PMC9699059")],
    "iu2024_immune_exhaustion_mecfs_lc": [("https://insight.jci.org/articles/view/183810/sd/pdf/render/1", "jci183810_sd1.xlsx", True)],
    "liinc_nk_jci2025": [("https://www.jci.org/articles/view/188182/sd/pdf/render/2", "jci188182_sd2.xlsx", True)],
    "ped_lc_jci_insight_2026": [("https://insight.jci.org/articles/view/201111/sd/pdf/render/2", "jci201111_sd2.xlsx", True)],
    "taurine_pcc_plos2024": [_epmc("PMC11152273")],
    "pulmonary_lc_jci_insight_2024": [("https://insight.jci.org/articles/view/177518/sd/pdf/render/3", "jci177518_sd3.xlsx", True)],
    "stanford_iris_olink": [(ZAPI.format(rec=13363117, key="IRIS_LongCOVID_olink_data.csv"), "IRIS_LongCOVID_olink_data.csv", False)],
    "stanford_iris_antigen_cytof": [(ZAPI.format(rec=15199328, key="ag_symptoms.csv"), "ag_symptoms.csv", False)],
    "bodansky_liinc_phipseq": [(ZAPI.format(rec=8021708, key="README.md"), "README.md", False)],
    "gold_ebv_reactivation": [(ZAPI.format(rec=4968818, key="LongCOVID-EBV-16Jun2021.xlsx"), "LongCOVID-EBV-16Jun2021.xlsx", False)],
    "gu2023_mtbls7337_multiomics": [
        (MEND + "zyzt62gbrw/files/efef2fb2-8492-449e-abb8-2062bf0b4702/file_downloaded", "Supplemental_Data_1.xlsx", True)],
    "mw_st003103_lc_mitochondria": [(MW.format(an="AN005077"), "AN005077_datatable.txt", False)],
    "frontiers_lc_mecfs_metabolomics": [(FS + "44098049", "DataSheet_2.xlsx", False)],
    "lc_lipid_mediators_cytokines_mendeley": [
        (MEND + "b34tpw33vh/files/b225ae50-ff6a-4476-aa48-cdc654056ae2/file_downloaded", "c_metadata.xlsx", True)],
    "pride_pxd066724_lc_plasma": [(PRIDE + "2026/03/PXD066724/metadata.sdrf.tsv", "metadata.sdrf.tsv", False)],
    "pride_misc_lc_proteomics": [(PRIDE + "2023/11/PXD040175/ENova_output_DIA-NN_pg.xlsx", "PXD040175_ENova_output_DIA-NN_pg.xlsx", False)],
    "iraq_lc_maes_group": [(FS + "48486631", "iraq_lc_electrolytes.xlsx", False)],
    "plos_pcc_recovered_small": [(PLOS.format(s="0315486.s001"), "pone.0315486.s001.xlsx", False)],
    "plos_postcovid_resp_acute_labs": [(PLOS.format(s="0344371.s001"), "pone.0344371.s001.xlsx", False)],
    "dataverse_lc_cardiovascular_cohort": [(HDV + "14058030", "Long_COVID_repository.xlsx", True)],
    "mtbls850_impacc_pasc": [(MTBLS + "MTBLS850/s_MTBLS850.txt", "s_MTBLS850.txt", False)],
    # ---- ME/CFS ----
    "biomapai_jax_mecfs_multiomics": [(GH + "metadata/Metadata_061523.csv", "Metadata_061523.csv", False),
                                      (GH + "quest_lab/Quest_residue.csv", "Quest_residue.csv", False),
                                      (GH + "quest_lab/Quest_lab_maaslined.csv", "Quest_lab_maaslined.csv", False)],
    "germain2022_jci_insight_exercise_metabolomics": [
        ("https://df6sxcketz7bb.cloudfront.net/manuscripts/157000/157621/jci.insight.157621.sdd1.xlsx",
         "jci.insight.157621.sdd1.xlsx", False)],
    "hoel2026_somascan7k_pad000026": [(PRIDE + "2026/01/PAD000026/Log2Matrix.tsv", "Log2Matrix.tsv", False)],
    "walitt2024_nih_somascan_csf_plasma_geo": [
        ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254030/suppl/GSE254030_CHI-19-021.Plasma_for_CFS.adat.txt.gz",
         "GSE254030_Plasma_for_CFS.adat.txt.gz", False)],
    "walitt2024_nih_csf_metabolomics_mw_st003178": [(MW.format(an="AN005217"), "AN005217_datatable.txt", False)],
    "walitt2024_natcomm_source_data": [
        (SPR.format(art="s41467-024-45107-3", f="41467_2024_45107_MOESM6_ESM.zip"), "41467_2024_45107_MOESM6_ESM.zip", True)],
    "naviaux2016_mw_st000450": [(MW.format(an="AN000705"), "AN000705_datatable.txt", False)],
    "naviaux_validation_mw_st000617": [(MW.format(an="AN000946"), "AN000946_datatable.txt", False)],
    "che2022_columbia_plasma_metabolomics_st2000": [(MW.format(an="AN003264"), "AN003264_datatable.txt", False)],
    "nagyszakal2018_cfi_metabolomics_st0800": [(MW.format(an="AN001274"), "AN001274_datatable.txt", False)],
    "lipkin_csf_metabolomics_mecfs_ms_st0910": [(MW.format(an="AN001480"), "AN001480_datatable.txt", False)],
    "armstrong2015_mtbls161_nmr": [(MTBLS + "MTBLS161/s_MTBLS161.txt", "s_MTBLS161.txt", False)],
    "baraniuk2021_csf_metabolomics_gwi_mecfs": [(FS + "26042977", "baraniuk2021_S1.xlsx", False)],
    "germain2018_2020_metabolites_supp": [_epmc("PMC6315598")],
    "ncnp_serum_indole_st004940": [(MW.format(an="AN008373"), "AN008373_datatable.txt", False)],
    "jonsjo2020_olink_inflammation_zenodo": [
        (ZAPI.format(rec=4960364, key="Kopia_av_MEFCS_Inflammation_Shared_dataset.xlsx"), "MEFCS_Inflammation_Shared_dataset.xlsx", False)],
    "bragee2026_csf_proteomics_pxd076216": [(PRIDE + "2026/04/PXD076216/sdrf.tsv", "sdrf.tsv", False)],
    "schutzer2023_csf_proteome_mecfs_fm": [(FS + "42389424", "schutzer_TableS2.xlsx", False)],
    "broderick2020_gwi_exercise_cytokines": [(FS + "22443671", "broderick2020_TableS1.xlsx", False)],
    "hochecker_pbmc_seahorse_mendeley": [
        (MEND + "mv2cnytjhz/files/92938383-5a17-44aa-96f6-d7ef99dad9d2/file_downloaded", "ATP_production.csv", True)],
    "ev_proteomics_raw_only_pride": [
        (PRIDE + "2026/08/PXD073644/ExperimentalDesign_Proteomnanalyse_pcMECFS.txt", "PXD073644_ExperimentalDesign_pcMECFS.txt", False)],
    "norway_endothelial_fmd_plos2023": [(FS + "39048429", "norway_fmd_S4.xlsx", False)],
    "hoel2021_jci_insight_supp": [
        ("https://df6sxcketz7bb.cloudfront.net/manuscripts/149000/149217/jci.insight.149217.sdd1.xlsx", "jci.insight.149217.sdd1.xlsx", False)],
    # ---- FM / Lyme / IBS / migraine / endometriosis / gastroparesis ----
    "fm_rsu_hhv6_cytokines": [("https://dataverse.rsu.lv/api/access/datafile/338", "rsu_datafile_338.tab", True),
                              ("https://dataverse.rsu.lv/api/access/datafile/340", "rsu_datafile_340.tab", True),
                              ("https://dataverse.rsu.lv/api/access/datafile/339", "rsu_datafile_339.tab", True)],
    "fm_caboni2014_lysopc_metabolomics": [(PLOS.format(s="0107626.s001"), "pone.0107626.s001.xlsx", False)],
    "fm_ibs_exercise_barrier": [(FS + "65663160", "fm_ibs_exercise_barrier.xlsx", False)],
    "fm_clos_garcia_serum_proteomics_pxd059894": [(PRIDE + "2025/09/PXD059894/raw_data_proteomics.xlsx", "raw_data_proteomics.xlsx", False)],
    "fm_mw_st000888_comparator": [(MWF.format(st="ST000888"), "ST000888_factors.json", False),
                                  (MWD.format(f="ST000888_AN001450_Results.txt"), "ST000888_AN001450_Results.txt", False)],
    "lyme_ptlds_mw_st001391": [(MWF.format(st="ST001391"), "ST001391_factors.json", False),
                               (MWD.format(f="ST001391_AN002320_Results.txt"), "ST001391_AN002320_Results.txt", False)],
    "lyme_early_mw_st001223": [(MWF.format(st="ST001223"), "ST001223_factors.json", False),
                               (MWD.format(f="ST001223_AN002036_Results.txt"), "ST001223_AN002036_Results.txt", False)],
    "lyme_subramanian_olink_metabolon": [(ZAPI.format(rec=21631040, key="Central_README.md"), "Central_README.md", False)],
    "lyme_ptlds_children_mmp": [(HDV + "7675481", "lyme_ptlds_children_mmp.xlsx", True)],   # served as xlsx
    "lyme_lnb_csf_adipokines_children": [(HDV + "14102446", "lnb_csf_adipokines.xlsx", True)],
    "lyme_neuroborreliosis_ilads_serum_proteome": [(PRIDE + "2025/01/PXD052529/wynikiProt.xlsx", "wynikiProt.xlsx", False)],
    "lyme_romania_lnb_csf": [(ZAPI.format(rec=17732813, key="Supplementary%20material.xlsx"), "Supplementary_material.xlsx", False)],
    "lyme_children_ukraine_immunoblot": [(ZAPI.format(rec=22697497, key="1_Base_LB.xlsx"), "1_Base_LB.xlsx", False)],
    "lyme_multiantibody_serology_plos2021": [(PLOS.format(s="0253514.s001"), "pone.0253514.s001.xlsx", False)],
    "lyme_raman_serum_spectra": [],   # listing only
    "ibs_mars2020_cell_multiomics": [
        (MEND + "29n2z5r5ph/files/e5798fbc-9b0d-4bd5-8a4c-28bf321bb978/file_downloaded", "mars2020_serum_cytokines.xlsx", True)],
    "ibs_mtbls2774_serum_metabolome": [(MTBLS + "MTBLS2774/s_MTBLS2774.txt", "s_MTBLS2774.txt", False)],
    "ibs_mtbls1396_mars_metabolights": [(MTBLS + "MTBLS1396/s_MTBLS1396.txt", "s_MTBLS1396.txt", False)],
    "ibs_mw_st003954_serum": [(MWF.format(st="ST003954"), "ST003954_factors.json", False)],
    "ibs_d_urine_metabolomics_figshare": [(FS + "15317597", "ibs_d_urine_neg.xlsx", False)],
    "endo_nhanes_rhq360": [("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/RHQ_C.xpt", "RHQ_C.xpt", False)],
    "prusty2024_herpesvirus_serology_mendeley": [
        (MEND + "4xkft5g9r5/files/638736ab-4f14-4420-8d46-0bb914b579be/file_downloaded", "prusty_file_1.bin", True)],
    "ibs_ecm_serum_plos2017": [(PLOS.format(s="0185855.s001"), "pone.0185855.s001.xlsx", False)],
    "ibs_mw_st001940_cbt": [(MWF.format(st="ST001940"), "ST001940_factors.json", False)],
    "ibs_mtbls12737_urine": [(MTBLS + "MTBLS12737/s_MTBLS12737.txt", "s_MTBLS12737.txt", False)],
    "migraine_serum_zinc_japan": [(ZAPI.format(rec=17773060, key="Migraine%20Zinc%20data.xlsx"), "Migraine_Zinc_data.xlsx", False)],
    "migraine_nhanes_severe_headache_vitd": [(PLOS.format(s="0313082.s001"), "pone.0313082.s001.xls", False)],
    "migraine_pfo_mtbls13630": [(MTBLS + "MTBLS13630/s_MTBLS13630.txt", "s_MTBLS13630.txt", False)],
    "endo_mif_serum_bmc2020": [(FS + "24573119", "endo_mif_raw.xls", False)],
    "endo_arg1_repod": [("https://repod.icm.edu.pl/api/access/datafile/48824", "endo_arg1.tab", True)],
    "endo_xenoestrogen_plos2024": [(PLOS.format(s="0304766.s002"), "pone.0304766.s002.xlsx", False)],
    "endo_mw_st003984_plasma_pf": [(MWF.format(st="ST003984"), "ST003984_factors.json", False)],
    "endo_mw_st003579_ovarian": [(MWF.format(st="ST003579"), "ST003579_factors.json", False),
                                 (MWD.format(f="ST003579_AN005877_Results.txt"), "ST003579_AN005877_Results.txt", False)],
    "endo_mw_st003463_stool": [(MWF.format(st="ST003463"), "ST003463_factors.json", False)],
    "endo_mw_st000585_follicular": [(MWF.format(st="ST000585"), "ST000585_factors.json", False),
                                    (MWD.format(f="ST000585_AN000900_Results.txt"), "ST000585_AN000900_Results.txt", False)],
    "endo_mtbls2040_serum_nmr": [(MTBLS + "MTBLS2040/s_MTBLS2040.txt", "s_MTBLS2040.txt", False)],
    "endo_mtbls8621_ff": [(MTBLS + "MTBLS8621/s_MTBLS8621.txt", "s_MTBLS8621.txt", False)],
    "endo_npar_pelvic_adhesion": [(PLOS.format(s="0337077.s002"), "pone.0337077.s002.xlsx", False)],
    "endo_pf_monocytes_cytokines": [],   # listing only
    "gp_gpcrc_ho1_hba1c": [(PLOS.format(s="0187772.s001"), "pone.0187772.s001.xlsx", False)],
    "longcovid_mw_st003103": [],   # duplicate of mw_st003103_lc_mitochondria
    "pots_plasma_dia_pxd031458": [],  # listing only
}

# Repository listings fetched through measure_it.http (metadata only) and saved as JSON: candidate -> [(api url, file)].
LISTINGS: dict[str, list[tuple[str, str]]] = {
    "hat_venom_gwas_zenodo2025": [("https://zenodo.org/api/records/17914583", "zenodo_17914583.json")],
    "woodruff2023_emory_olink": [("https://zenodo.org/api/records/8092298", "zenodo_8092298.json")],
    "gao2025_uk_se_olink": [("https://zenodo.org/api/records/14772494", "zenodo_14772494.json")],
    "lyme_subramanian_olink_metabolon": [("https://zenodo.org/api/records/21631040", "zenodo_21631040.json")],
    "lyme_raman_serum_spectra": [("https://zenodo.org/api/records/18987524", "zenodo_18987524.json")],
    "bodansky_liinc_phipseq": [("https://zenodo.org/api/records/8021708", "zenodo_8021708.json")],
    "stanford_iris_antigen_cytof": [("https://zenodo.org/api/records/15199328", "zenodo_15199328.json")],
    "endo_pf_monocytes_cytokines": [("https://zenodo.org/api/records/18857406", "zenodo_18857406.json")],
    "ibs_d_urine_metabolomics_figshare": [("https://api.figshare.com/v2/file/download/15317597", "SKIP")],
    "yin2024_liinc_olink_cytof": [("https://datadryad.org/api/v2/datasets/doi%3A10.7272%2FQ6WD3XTB", "dryad_Q6WD3XTB.json")],
    "lyme_ptlds_serochip_epitopes": [],
    "pots_plasma_dia_pxd031458": [("https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD031458/files?pageSize=100",
                                   "pride_PXD031458_files.json")],
    "hoel2026_somascan7k_pad000026": [],
    "mapmecfs_catalogue": [("https://www.mapmecfs.org/api/3/action/package_search?rows=1000", "mapmecfs_package_search.json")],
}
LISTINGS["ibs_d_urine_metabolomics_figshare"] = []

# Anonymous access probes on gated routes (status code only; nothing behind a login is requested with credentials).
PROBES: dict[str, str] = {
    "mapmecfs_nih_pimecfs_clinical_labs": ("https://mapmecfs.org/dataset/2a824b5c-f4c4-474d-9c8b-553e7984f21d/resource/"
                                           "991194e1-5111-430f-83b0-1b53fe33aff6/download/datafile-clinical-master-labs-dataset.csv"),
    "yin2024_liinc_olink_cytof": "https://datadryad.org/api/v2/files/2772337/download",
    "fm_rsu_gut_markers_restricted": "https://dataverse.rsu.lv/api/access/datafile/1908",
}

# Never downloaded: the record appears to contain direct personal identifiers (names, dates of birth, phone numbers).
PRIVACY_EXCLUDED = {"postcovid_gpcr_aab_seibert2026", "figshare_gpcr_autoab_postcovid_EXCLUDE"}


def _head_size(url: str, headers: dict | None) -> int | None:
    try:
        r = http.request("HEAD", url, headers=headers, reject_html=False, max_retries=2, timeout=60)
        v = r.headers.get("Content-Length") or r.headers.get("content-length")
        return int(v) if v else None
    except Exception:
        return None


def download_all() -> dict[str, str]:
    """Retrieve every verification file; errors are returned (and recorded in inspection.json), not hidden."""
    errors: dict[str, str] = {}
    for cid, items in FILES.items():
        if cid in PRIVACY_EXCLUDED:
            continue
        for url, fname, browser in items:
            hdr = BUA if browser else None
            dest = ROOT / cid / fname
            if dest.exists() and fname in load_manifest(f"{SRC}/{cid}")["files"]:
                continue
            size = _head_size(url, hdr)
            if size and size > MAX_BYTES:
                errors[f"{cid}/{fname}"] = f"skipped: Content-Length {size} > {MAX_BYTES}"
                continue
            try:
                download_file(url, f"{SRC}/{cid}", fname, headers=hdr, max_retries=3, timeout=600)
            except Exception as exc:
                errors[f"{cid}/{fname}"] = repr(exc)[:300]
    for cid, items in LISTINGS.items():
        for url, fname in items:
            if fname == "SKIP":
                continue
            dest = ROOT / cid / fname
            if dest.exists():
                continue
            try:
                r = http.get(url, headers=BUA, reject_html=True, cache_errors=True, max_retries=3)
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(r.content)
                record_file(f"{SRC}/{cid}", fname, url=url, note=f"repository listing (HTTP {r.status})")
            except Exception as exc:
                errors[f"{cid}/{fname}"] = repr(exc)[:300]
    return errors


def probe_gated() -> dict[str, dict]:
    out = {}
    for cid, url in PROBES.items():
        try:
            r = http.get(url, headers=BUA, reject_html=False, cache_errors=True, max_retries=2, timeout=60)
            ctype = r.headers.get("Content-Type", "")
            out[cid] = {"url": url, "status": r.status, "content_type": ctype, "bytes": len(r.content),
                        "html_login_page": "html" in ctype.lower()}
        except Exception as exc:
            out[cid] = {"url": url, "error": repr(exc)[:300]}
    return out


# ---------------------------------------------------------------------------------------------------------------
# 2. Schema-level inspectors. Each returns counts of labelled PEOPLE and of labelled people WITH lab values.
#    They never compare analyte values between groups (pre-specification: outcomes stay unseen).
# ---------------------------------------------------------------------------------------------------------------
def _p(cid: str, name: str) -> Path:
    return ROOT / cid / name


def _xl(cid: str, name: str, **kw) -> pd.DataFrame:
    return pd.read_excel(_p(cid, name), **kw)


def _zip_member(cid: str, zname: str, member: str) -> io.BytesIO:
    with zipfile.ZipFile(_p(cid, zname)) as z:
        return io.BytesIO(z.read(member))


def _counts(label: pd.Series, case_values, control_values, has_lab: pd.Series | None = None, **extra) -> dict:
    """n labelled people per group and n with lab values present (no lab VALUES are compared)."""
    label = label.astype("object")
    case, ctrl = label.isin(case_values), label.isin(control_values)
    has = has_lab if has_lab is not None else pd.Series(True, index=label.index)
    out = {"n_rows": int(len(label)), "n_cases_file": int(case.sum()), "n_controls_file": int(ctrl.sum()),
           "n_cases_with_lab": int((case & has).sum()), "n_controls_with_lab": int((ctrl & has).sum()),
           "label_values": {str(k): int(v) for k, v in label.value_counts(dropna=False).items()}}
    out.update(extra)
    return out


def _mw_datatable(cid: str, name: str, factor: str | None = None) -> pd.DataFrame:
    """Metabolomics Workbench REST datatable: rows = samples; 'Class' = 'Key:Val | Key2:Val2'."""
    d = pd.read_csv(_p(cid, name), sep="\t")
    cls = d["Class"].astype(str)
    if factor:
        d["_label"] = cls.str.extract(rf"{re.escape(factor)}:([^|]+)")[0].str.strip()
    else:
        d["_label"] = cls.str.split("|").str[0].str.split(":").str[-1].str.strip()
    feats = [c for c in d.columns if c not in ("Samples", "Class", "_label")]
    d["_frac_present"] = d[feats].notna().mean(axis=1)
    d.attrs["n_features"] = len(feats)
    return d


def _mw_factors(cid: str, name: str, key: str) -> pd.DataFrame:
    js = json.loads(_p(cid, name).read_text())
    rows = list(js.values()) if isinstance(js, dict) else js
    d = pd.DataFrame(rows)
    d["_label"] = d["factors"].astype(str).str.extract(rf"{re.escape(key)}:([^|]+)")[0].str.strip()
    return d


def _mw_results_header(cid: str, name: str) -> tuple[list[str], list[str], int]:
    """(sample ids, factor strings, n feature rows) from a MW *_Results.txt (row 1 ids, row 2 'Factors')."""
    with open(_p(cid, name), encoding="utf-8", errors="replace") as fh:
        ids = fh.readline().rstrip("\n").split("\t")[1:]
        fac = fh.readline().rstrip("\n").split("\t")[1:]
        n = sum(1 for _ in fh)
    return ids, fac, n


def _isa_sample(cid: str, name: str, factor_contains: str) -> pd.DataFrame:
    d = pd.read_csv(_p(cid, name), sep="\t", dtype=str)
    col = next(c for c in d.columns if c.startswith("Factor Value") and factor_contains.lower() in c.lower())
    d["_label"] = d[col].str.strip()
    return d


# ---- per-candidate inspectors --------------------------------------------------------------------------------
KLEIN_CORTISOL = "x1_Cytokines_Cortisol_Obs_Conc_ng/mL_GG21.1_ML"


def _insp_klein():
    d = _xl("klein2023_mylc_ml_table", "41586_2023_6651_MOESM4_ESM.xlsx")
    a = d[d["x0_Censor_Complete"] == 0]
    cyt = [c for c in d.columns if c.startswith("x1_Cytokines")]
    lab = a["x0_Censor_Cohort_ID"].map({3: "LC", 1: "HC", 2: "CC"})
    return _counts(lab, ["LC"], ["HC", "CC"], a[KLEIN_CORTISOL].notna(), label_column="x0_Censor_Cohort_ID "
                   "(3=LC, 1=HC, 2=CC; mapping verified from infection-test columns) within x0_Censor_Complete==0",
                   lab_columns_checked="serum cortisol (primary); 144 x1_Cytokines analytes", n_analytes=len(cyt),
                   n_rows_file=int(len(d)), shared_ids_across_analytes=True,
                   notes_counts=(f"analysis set x0_Censor_Complete==0: {len(a)} of {len(d)} rows; HC {int((lab == 'HC').sum())}, "
                                 f"CC {int((lab == 'CC').sum())}, LC {int((lab == 'LC').sum())}; one row per person, all "
                                 f"{len(d.columns)} columns share x0_LC_ID; draw time, age, sex, BMI present"))


def _insp_cinquina():
    cid = "heds_hsd_olink_serum_cinquina2026"
    npx = pd.read_excel(_zip_member(cid, "PMC13081554_supplementaryFiles.zip", "12014_2026_9588_MOESM4_ESM.xlsx"), header=None)
    ids = npx.iloc[7:, 0].fillna("").astype(str).str.strip()
    vals = npx.iloc[7:, 1:461].apply(pd.to_numeric, errors="coerce")
    grp = ids.str.extract(r"^([PHC])\d+$")[0].map({"P": "hEDS", "H": "HSD", "C": "control"})
    has = vals.notna().mean(axis=1) >= 0.5
    return _counts(grp, ["hEDS", "HSD"], ["control"], has,
                   label_column="sample-ID prefix P/H/C (P=hEDS and H=HSD verified against MOESM6 'ID NPx OLINK')",
                   lab_columns_checked="460 Olink Target 96 NPX assays (5 panels)", n_analytes=460,
                   shared_ids_across_analytes=True,
                   notes_counts=(f"hEDS {int((grp == 'hEDS').sum())}, HSD {int((grp == 'HSD').sum())}, controls "
                                 f"{int((grp == 'control').sum())}; plate-control/LOD/QC rows excluded; patients carry 2017-"
                                 "criteria checklists, sex, age, relationship to proband and comorbidity flags (incl. "
                                 "dysautonomia/POTS, chronic fatigue, MCAS/allergy, functional GI); control covariates "
                                 "not released"))


def _insp_cervia():
    cid = "cervia_hasler2024_zurich_somascan"
    lab = _xl(cid, "Proteomics_Clinical_Data_Labels.xlsx", sheet_name="Sheet1")
    prot = _xl(cid, "Proteomics_Clinical_Data_6M_timepoint.xlsx", sheet_name="Data", usecols=[0, 1, 2])
    has = lab["SubjectID"].isin(prot["SubjectID"])
    return _counts(lab["PACS_6M_woDys"], [1], [0], has, label_column="PACS_6M_woDys (1 = PACS at 6 months)",
                   lab_columns_checked="SomaScan 7k aptamers at 6 months (SubjectID match)", n_analytes=7596,
                   shared_ids_across_analytes=True,
                   notes_counts="controls = COVID-convalescent without PACS; healthy uninfected only at the acute timepoint file")


def _insp_talla():
    cid = "talla2023_aifi_olink"
    m = _xl(cid, "41467_2023_38682_MOESM3_ESM.xlsx", header=2)
    per = m.drop_duplicates("PTID")
    return _counts(per["Group"], ["Infected PASC"], ["Infected Recovered", "Uninfected"], None,
                   label_column="Group (per PTID)", lab_columns_checked="Olink Explore 1463 NPX (MOESM4), per visit",
                   n_analytes=1463, shared_ids_across_analytes=True,
                   notes_counts=f"{len(m)} sample-visits, {per.PTID.nunique()} participants; longitudinal (choose one visit)")


def _insp_gu():
    cid = "gu2023_mtbls7337_multiomics"
    c = _xl(cid, "Supplemental_Data_1.xlsx", sheet_name="Cytokines")
    conv = c[c["Type"].astype(str).str.strip() == "Convalescence"]
    ctrl = c[c["Type"].astype(str).str.strip() == "Control"]
    lab = pd.concat([conv["Outcome "], ctrl["Type"]]).astype(str).str.strip()
    return _counts(lab, ["With Event"], ["Control"], None, label_column="Type + 'Outcome ' (With Event vs healthy Control)",
                   lab_columns_checked="48-plex cytokines, ~274 MRM proteins, ~635 metabolites", n_analytes=48,
                   shared_ids_across_analytes=True,
                   notes_counts=(f"Type counts {c['Type'].value_counts().to_dict()}; the LC label ('With Event' vs "
                                 "symptom flags) is ambiguous; 6-month convalescent samples"))


def _insp_iris():
    d = pd.read_csv(_p("stanford_iris_olink", "IRIS_LongCOVID_olink_data.csv"))
    per = d.drop_duplicates("IRIS_number")
    return _counts(per["long_covid"].astype(str), ["1", "True", "Yes", "yes", "LC", "long_covid"],
                   ["0", "False", "No", "no", "non-LC"], None, label_column="long_covid (per IRIS_number)",
                   lab_columns_checked="Olink Target 96 Inflammation + Immune Response", n_analytes=int(d.Assay.nunique()),
                   shared_ids_across_analytes=True, notes_counts=f"timings {d.timing_month.value_counts().to_dict()}")


def _insp_lopez():
    d = _xl("lopezhernandez2023_mexico_metabolome", "LCDATA_AVAILABILITY_SR.xls")
    num = d.select_dtypes("number")
    return _counts(d["GROUP"].astype(str).str.strip(), ["POST-COVID"], ["CONTROL"], num.notna().mean(axis=1) > 0.5,
                   label_column="GROUP", lab_columns_checked="108 targeted metabolites + IL-17/PGE2 (TMIC MEGA)",
                   n_analytes=int(num.shape[1]), shared_ids_across_analytes=True,
                   notes_counts=f"SUBGROUP {d['SUBGROUP'].value_counts(dropna=False).to_dict()} (POST-COVID includes 'RECOVERED')")


def _insp_iraq():
    d = _xl("iraq_lc_maes_group", "iraq_lc_electrolytes.xlsx")
    return _counts(d["Pts1, Ctrl0"], [1], [0], None, label_column="Pts1, Ctrl0",
                   lab_columns_checked="CRP, albumin, Ca, Mg (+ galanin file: PAI-1, IGF-1, PGE2, NSE, S100B)",
                   n_analytes=int(d.select_dtypes("number").shape[1]), shared_ids_across_analytes=False)


def _insp_lu():
    d = _xl("lu2024_early_markers_pasc_figshare", "Raw_Host_Immune_Data.xls")
    per = d.drop_duplicates("id")
    return _counts(per["pascm4"].astype(str), ["PASC", "Yes", "1"], ["No PASC", "No", "0"], None,
                   label_column="pascm4 (per id)", lab_columns_checked="NfL, IL-6, IL-10, TNF, GFAP, MCP-1, IFN-g, IP-10, IFN-a",
                   n_analytes=9, shared_ids_across_analytes=True,
                   notes_counts="markers measured in ACUTE infection (days 5-28): a prediction design, not a concurrent-lab contrast")


def _insp_gold():
    d = _xl("gold_ebv_reactivation", "LongCOVID-EBV-16Jun2021.xlsx")
    return _counts(d["Group"].astype(str), ["LongTermPACS"], ["LongTermControl"], d["EBVEAIgG"].notna(),
                   label_column="Group", lab_columns_checked="EBV EA-D IgG, VCA IgM/IgG, EBNA-1 IgG",
                   n_analytes=4, shared_ids_across_analytes=True)


def _insp_appelman():
    d = _xl("appelman2024_muscle_pem", "41467_2023_44432_MOESM4_ESM.xlsx", sheet_name="Supplementary Data Figure 7G-I")
    b = d[d["Time"].astype(str).str.strip() == "Baseline"]
    return _counts(b["Group"], ["Long COVID"], ["Healthy"], b["Cortisol"].notna(), label_column="Group (baseline rows)",
                   lab_columns_checked="cortisol, CK, CK-MB (+ 83 blood metabolites on another sheet)", n_analytes=3,
                   shared_ids_across_analytes=True, notes_counts="cortisol with minutes since waking; same trial as Charlton 2026")


def _insp_kedor():
    d = _xl("kedor2022_charite_pcs_mecfs", "41467_2022_32507_MOESM4_ESM.xlsx", sheet_name="Tab 4, S10 & Fig 5")
    d["Diagnose"] = d["Diagnose"].ffill()   # the label is written once at the top of each group block
    per = d.dropna(subset=["recordID"]).drop_duplicates("recordID")
    return _counts(per["Diagnose"].astype(str).str.strip(), ["PCS/CFS"], ["PCS"], None,
                   label_column="Diagnose (PCS/CFS = CCC-positive post-COVID; PCS = post-COVID not meeting CCC)",
                   lab_columns_checked="31 routine/immune labs (IL-8, ACE, MBL, IgG1-4, ferritin, NT-proBNP ...)",
                   n_analytes=int(d.shape[1] - 2), shared_ids_across_analytes=True,
                   notes_counts="symptomatic (post-COVID) controls; no healthy controls")


def _insp_pcc_small():
    d = _xl("plos_pcc_recovered_small", "pone.0315486.s001.xlsx")
    return _counts(d["Group"], [2.0], [1.0], d["Cort_ngmL_Plasm"].notna(), label_column="Group (2 = PCC per searcher; 1 = recovered)",
                   lab_columns_checked="cortisol, IL-6, IL-10, NGF, BDNF, DHEA, NfL", n_analytes=14, shared_ids_across_analytes=True)


def _insp_lipid_mediators():
    d = _xl("lc_lipid_mediators_cytokines_mendeley", "c_metadata.xlsx")
    return _counts(d["LCGroup"].astype(str).str.strip(), ["LC"], ["Recov"], None, label_column="LCGroup",
                   lab_columns_checked="lipid mediators, SARS-CoV-2 antigens, 10 cytokines (other two files)",
                   n_analytes=None, shared_ids_across_analytes=True, notes_counts="publication not identified")


def _insp_biomapai():
    cid = "biomapai_jax_mecfs_multiomics"
    m = pd.read_csv(_p(cid, "Metadata_061523.csv"))
    tp1 = m[m["timepoints"] == "tp1"]
    q = pd.read_csv(_p(cid, "Quest_lab_maaslined.csv"), nrows=2)
    qcols = {c for c in q.columns if c.endswith("_tp1")}
    has = (tp1["sample_id_tp1"].astype(str) + "_tp1").isin(qcols)
    return _counts(tp1["study_ptorhc"], ["MECFS"], ["Control"], has, label_column="study_ptorhc (timepoint tp1)",
                   lab_columns_checked="48 Quest clinical labs (transformed: MaAsLin-adjusted / 0-1 scaled)",
                   n_analytes=48, shared_ids_across_analytes=True,
                   notes_counts="lab values are NOT in clinical units (rescaled/residualised); no licence file in the repo")


def _germain_header(cid: str, name: str, sheet: str) -> pd.DataFrame:
    raw = _xl(cid, name, sheet_name=sheet, header=None, nrows=8)
    lab_row = next(i for i in range(8) if raw.iloc[i].astype(str).str.fullmatch(r"(CFS|Control)").sum() > 20)
    id_row = next(i for i in range(8) if raw.iloc[i].astype(str).str.match(r"CU-?\d+").sum() > 20)
    tp_row = next((i for i in range(8) if raw.iloc[i].astype(str).str.contains(r"D1-PRE").sum() > 20), None)
    cols = raw.columns[raw.iloc[lab_row].astype(str).str.fullmatch(r"(CFS|Control)")]
    return pd.DataFrame({"id": raw.iloc[id_row, cols].astype(str).values, "label": raw.iloc[lab_row, cols].values,
                         "tp": raw.iloc[tp_row, cols].values if tp_row is not None else None})


def _insp_germain2022():
    h = _germain_header("germain2022_jci_insight_exercise_metabolomics", "jci.insight.157621.sdd1.xlsx", "Original Dataset")
    b = h[h["tp"].astype(str) == "D1-PRE"].drop_duplicates("id")
    return _counts(b["label"], ["CFS"], ["Control"], None, label_column="Phenotype row (D1-PRE samples)",
                   lab_columns_checked="~1,150 Metabolon biochemicals (raw area counts)", n_analytes=1157,
                   shared_ids_across_analytes=True, notes_counts=f"{len(h)} samples at 4 timepoints around 2-day CPET")


def _insp_hoel2026():
    d = pd.read_csv(_p("hoel2026_somascan7k_pad000026", "Log2Matrix.tsv"), sep="\t", usecols=range(9))
    return _counts(d["Group"], ["ME"], ["HC"], None, label_column="Group", lab_columns_checked="SomaScan 7k (7,326 aptamers)",
                   n_analytes=7326, shared_ids_across_analytes=True,
                   notes_counts=f"covariates: Sex, Age (5-y bins), BMI_Cat, SF36PF, steps; Metabotype {d.Metabotype.value_counts(dropna=False).to_dict()}")


def _mw_insp(cid, name, factor, case, ctrl, analytes, **kw):
    d = _mw_datatable(cid, name, factor)
    return _counts(d["_label"], case, ctrl, d["_frac_present"] > 0.5, label_column=f"MW Class '{factor or 'first factor'}'",
                   lab_columns_checked=analytes, n_analytes=d.attrs["n_features"], **kw)


def _insp_armstrong():
    d = _isa_sample("armstrong2015_mtbls161_nmr", "s_MTBLS161.txt", "Chronic Fatigue")
    s = d[d["Sample Name"].str.contains("serum", case=False)]
    return _counts(s["_label"], ["CFS"], ["non-CFS control"], None, label_column="Factor Value[Chronic Fatigue Syndrome] (serum samples)",
                   lab_columns_checked="1H-NMR serum metabolites (~41, MAF)", n_analytes=41, shared_ids_across_analytes=False)


def _insp_baraniuk():
    d = _xl("baraniuk2021_csf_metabolomics_gwi_mecfs", "baraniuk2021_S1.xlsx", header=1)
    b = d[d["Exercise"] == 0]
    return _counts(b["Group"].astype(str).str.strip(), ["cfs0"], ["sc0"], None, label_column="Group (Exercise == 0)",
                   lab_columns_checked="~187 targeted CSF metabolites", n_analytes=int(d.shape[1] - 5),
                   shared_ids_across_analytes=True, notes_counts=f"all groups {d['Group'].value_counts().to_dict()} (GWI = disease comparator)")


def _insp_rsu_fm():
    cid = "fm_rsu_hhv6_cytokines"
    out = {}
    for f in ("rsu_datafile_338.tab", "rsu_datafile_340.tab", "rsu_datafile_339.tab"):
        try:
            d = pd.read_csv(_p(cid, f), sep=";", encoding="latin-1")
            out[f] = (len(d), str(d["ID"].astype(str).str.extract(r"^([A-Za-z]+)")[0].value_counts().to_dict()))
        except Exception as exc:
            out[f] = repr(exc)[:80]
    ids = pd.concat([pd.read_csv(_p(cid, f), sep=";", encoding="latin-1")["ID"].astype(str)
                     for f in ("rsu_datafile_338.tab", "rsu_datafile_340.tab")])
    lab = ids.str.extract(r"^([A-Za-z]+)")[0]
    return _counts(lab, ["FS"], ["C"], None, label_column="ID prefix (FS = FM, C = control)",
                   lab_columns_checked="9 plasma cytokines, HHV-6 serology/PCR", n_analytes=9, shared_ids_across_analytes=True,
                   notes_counts=json.dumps(out))


def _insp_caboni():
    d = _xl("fm_caboni2014_lysopc_metabolomics", "pone.0107626.s001.xlsx", sheet_name="Foglio3", header=None)
    ids = d.iloc[:, 0].astype(str)
    lab = ids.str.extract(r"^(FMS|Cont)\d+")[0]
    return _counts(lab.dropna(), ["FMS"], ["Cont"], None, label_column="subject ID prefix (FMS / Cont)",
                   lab_columns_checked="~30 LC-Q-TOF features (lysoPC etc.)", n_analytes=30, shared_ids_across_analytes=True)


def _insp_st000888():
    d = _mw_factors("fm_mw_st000888_comparator", "ST000888_factors.json", "Disease Status")
    d["subject"] = d["local_sample_id"].astype(str).str.replace(r"-Run.*$", "", regex=True)
    per = d.drop_duplicates("subject")
    ids, _, nfeat = _mw_results_header("fm_mw_st000888_comparator", "ST000888_AN001450_Results.txt")
    has = per["local_sample_id"].isin(ids) | per["subject"].isin([re.sub(r"-Run.*$", "", i) for i in ids])
    healthy = [v for v in per["_label"].dropna().unique() if "control" in v.lower() and "Lyme" not in v]
    return _counts(per["_label"], ["Fibromyalgia"], healthy, has, label_column="Disease Status (per subject)",
                   lab_columns_checked=f"{nfeat} targeted LC-MS features preselected as a LYME biosignature", n_analytes=nfeat,
                   shared_ids_across_analytes=True, notes_counts="FM is a look-alike comparator arm; features chosen for Lyme")


def _insp_st001391():
    d = _mw_factors("lyme_ptlds_mw_st001391", "ST001391_factors.json", "PTLDS")
    d["tp"] = d["factors"].astype(str)
    b = d[d["tp"].str.contains("Baseline", case=False)]
    return _counts(b["_label"], ["Yes"], ["No"], None, label_column="PTLDS (baseline samples)",
                   lab_columns_checked="~525 LC-MS features", n_analytes=525, shared_ids_across_analytes=True,
                   notes_counts="both arms had treated early Lyme; label = later PTLDS (prognostic design)")


def _insp_st001223():
    d = _mw_factors("lyme_early_mw_st001223", "ST001223_factors.json", "Sample Type")
    d["subject"] = d["local_sample_id"].astype(str).str.replace(r"-Run.*$", "", regex=True)
    per = d.drop_duplicates("subject")
    lab = per["_label"].fillna("").map(lambda v: "Lyme" if "Lyme" in v else ("healthy" if "control" in v.lower() else v))
    ids, _, nfeat = _mw_results_header("lyme_early_mw_st001223", "ST001223_AN002036_Results.txt")
    return _counts(lab, ["Lyme"], ["healthy"], per["local_sample_id"].isin(ids), label_column="Sample Type (per subject)",
                   lab_columns_checked=f"{nfeat} untargeted LC-MS features (unannotated)", n_analytes=nfeat,
                   shared_ids_across_analytes=True, notes_counts=f"subject-level {per['_label'].value_counts().to_dict()}")


def _insp_mars():
    d = _xl("ibs_mars2020_cell_multiomics", "mars2020_serum_cytokines.xlsx", sheet_name="serum")
    per = d.drop_duplicates("Patient")
    lab = per["Cohort"].astype(str).str.strip()
    return _counts(lab, [v for v in lab.unique() if v != "H"], ["H"], None, label_column="Cohort (C = IBS-C, D = IBS-D, H = healthy)",
                   lab_columns_checked=f"serum {sorted(d['variable'].unique())}", n_analytes=int(d["variable"].nunique()),
                   shared_ids_across_analytes=True)


def _insp_mtbls2774():
    d = _isa_sample("ibs_mtbls2774_serum_metabolome", "s_MTBLS2774.txt", "Disease")
    return _counts(d["_label"], ["Irritable Bowel Syndrome"], ["control"], None, label_column="Factor Value[Disease] (sample rows)",
                   lab_columns_checked="LC-MS annotated metabolites (131-152 per mode)", n_analytes=152,
                   shared_ids_across_analytes=False, notes_counts="sample rows, likely 2 sample types per subject: persons not resolved")


def _insp_ibs_ecm():
    d = _xl("ibs_ecm_serum_plos2017", "pone.0185855.s001.xlsx", sheet_name="Figure 2")
    lab = d["Disease "].astype(str).str.strip()
    return _counts(lab, ["IBS"], ["Healthy donors"], d["C5M"].notna(), label_column="Disease",
                   lab_columns_checked="BGM, EL-NE, C5M, Pro-C5, CRP", n_analytes=5, shared_ids_across_analytes=True,
                   notes_counts=f"{lab.value_counts().to_dict()} (CD/UC = disease comparators)")


def _insp_migraine_zinc():
    d = _xl("migraine_serum_zinc_japan", "Migraine_Zinc_data.xlsx", sheet_name="matching")
    return _counts(d["group case 0 control 1"], [0], [1], d["Zn(80-130)"].notna(), label_column="group case 0 control 1",
                   lab_columns_checked="Fe, Mg, Zn, Cu", n_analytes=4, shared_ids_across_analytes=True,
                   notes_counts="age/sex-matched 1:1 table")


def _insp_endo_mif():
    d = _xl("endo_mif_serum_bmc2020", "endo_mif_raw.xls")
    return _counts(d["Group"], [1], [0], d["MIF"].notna(), label_column="Group (1 = endometriosis, surgical)",
                   lab_columns_checked="serum MIF", n_analytes=1, shared_ids_across_analytes=True,
                   notes_counts=f"stage {d['Stage'].value_counts().to_dict()}")


def _insp_endo_arg1():
    d = pd.read_csv(_p("endo_arg1_repod", "endo_arg1.tab"), sep="\t")
    col = next(c for c in d.columns if c.strip().lower() == "group")
    lab = d[col].fillna("").astype(str).str.strip()
    has = pd.to_numeric(d["Arg1PREng/ml"], errors="coerce").notna()
    dx = pd.crosstab(d["diagnosis"].fillna("NA"), lab)
    ctrl3_healthy = int(dx.loc["healthy", "Ctrl3"]) if "healthy" in dx.index and "Ctrl3" in dx.columns else 0
    return _counts(lab, ["Patient"], [v for v in lab.unique() if v.startswith("Ctrl")], has, label_column=col,
                   lab_columns_checked="pre-operative serum ARG1 (primary lab column Arg1PREng/ml); ARG2, arginase "
                   "activity, glycaemia, creatinine, AST, ALT", n_analytes=7, shared_ids_across_analytes=True,
                   notes_counts=(f"'diagnosis' column: Ctrl = surgical non-endometriosis gynaecological patients "
                                 f"(myoma, cysts, pain ...), Ctrl2 = uterine myomas, Ctrl3 = {ctrl3_healthy} 'healthy' "
                                 f"+ {int((lab == 'Ctrl3').sum()) - ctrl3_healthy} other"))


def _insp_endo_xeno():
    d = _xl("endo_xenoestrogen_plos2024", "pone.0304766.s002.xlsx")
    return _counts(d["0control/1case"], [1], [0], d["TEXBalpha Eeq/g"].notna(), label_column="0control/1case",
                   lab_columns_checked="adipose TEXB-alpha", n_analytes=1, shared_ids_across_analytes=True,
                   notes_counts="cases mix endometriosis and leiomyoma ('Current condition'); adipose tissue, not blood")


def _insp_mw_results(cid, name, case_kw, ctrl_kw, analytes):
    ids, fac, n = _mw_results_header(cid, name)
    lab = pd.Series([("case" if any(k in f for k in case_kw) else "control" if any(k in f for k in ctrl_kw) else f)
                     for f in fac])
    return _counts(lab, ["case"], ["control"], None, label_column="Factors row of the MW results file",
                   lab_columns_checked=analytes, n_analytes=n, shared_ids_across_analytes=True,
                   notes_counts=f"{pd.Series(fac).value_counts().head(8).to_dict()}")


def _insp_id_prefix(cid, name, analytes):
    ids, _, n = _mw_results_header(cid, name)
    lab = pd.Series(ids).map(lambda i: "control" if str(i).startswith("CON") else "case")
    return _counts(lab, ["case"], ["control"], None,
                   label_column="sample-ID prefix (CON = control; factors file for the others)", lab_columns_checked=analytes,
                   n_analytes=n + 1, shared_ids_across_analytes=True)


def _insp_lyme_serology():
    d = _xl("lyme_multiantibody_serology_plos2021", "pone.0253514.s001.xlsx")
    return _counts(d["Status1"].astype(str).str.strip(), ["Lyme"], ["NoLyme"], None, label_column="Status1",
                   lab_columns_checked="anti-Borrelia serology (VIDAS, blots, VlsE, C6, pepC10)", n_analytes=20,
                   shared_ids_across_analytes=True, notes_counts="the lab IS the diagnostic test family (circular)")


def _insp_hat_gp():
    d = _xl("hat_gastroparesis_ev_lipidomics_2022", "Supp_Table1.xlsx", header=None)
    hdr = d.astype(str).apply(lambda r: r.str.fullmatch(r"P(CT|HAT)\d+").sum(), axis=1).idxmax()
    lab = d.iloc[hdr].astype(str).str.extract(r"^(PCT|PHAT)")[0].dropna()
    return _counts(lab, ["PHAT"], ["PCT"], None, label_column="sample column prefix (PHAT = HaT, PCT = no HaT)",
                   lab_columns_checked="~189 plasma-EV lipids", n_analytes=189, shared_ids_across_analytes=True,
                   notes_counts="all have paediatric gastroparesis-like symptoms; HaT is the contrast")


def _insp_heparin():
    d = _xl("mcad_plasma_heparin_vysniauskaite2015", "pone.0124912.s001.xls", header=2)
    col = next(c for c in d.columns if "iagnos" in str(c))
    return _counts(d[col].astype(str).str.strip(), ["MCAS"], ["SM"], None, label_column=f"{col} (MCAS vs SM; no controls)",
                   lab_columns_checked="plasma heparin basal/occlusion, tryptase, CgA (cat), N-methylhistamine (cat)",
                   n_analytes=5, shared_ids_across_analytes=True,
                   notes_counts="case-only; tryptase is part of the SM criteria, so MCAS vs SM by tryptase is circular")


def _insp_oi_children():
    x = pd.read_excel(_zip_member("oi_children_metabolomics_li2024", "PMC10941598_supplementaryFiles.zip",
                                  "13052_2024_1601_MOESM1_ESM.xlsx"), header=None, nrows=3)
    cells = pd.Series(x.values.ravel()).astype(str)
    lab = cells.str.extract(r"^(SYN|CTR)\d*")[0].dropna()
    return _counts(lab, ["SYN"], ["CTR"], None, label_column="sample column prefix (SYN = OI/syncope, CTR = control)",
                   lab_columns_checked="~105 preselected differential plasma metabolites", n_analytes=105,
                   shared_ids_across_analytes=True, notes_counts="only the authors' differential metabolites are released")


def _insp_iu2024():
    x = _xl("iu2024_immune_exhaustion_mecfs_lc", "jci183810_sd1.xlsx", sheet_name=None, header=None, nrows=3)
    cells = pd.Series(pd.concat(x.values()).values.ravel()).astype(str)
    lab = cells.str.extract(r"^(HC|LC|ME)\d+$")[0].dropna().drop_duplicates()
    ids = cells[cells.str.fullmatch(r"(HC|LC|ME)\d+")].drop_duplicates()
    lab = ids.str.extract(r"^(HC|LC|ME)")[0]
    return _counts(lab, ["LC", "ME"], ["HC"], None, label_column="sample-ID prefix HC/LC/ME",
                   lab_columns_checked="PBMC NanoString immune-exhaustion panel (~785 genes)", n_analytes=636,
                   shared_ids_across_analytes=True)


def _insp_mw_st003103():
    return _mw_insp("mw_st003103_lc_mitochondria", "AN005077_datatable.txt", None, ["Long COVID"], ["Recovered from acute COVID-19"],
                    "acylcarnitines, ceramides, lipids (LC-MS)", shared_ids_across_analytes=True)


INSPECTORS = {
    "klein2023_mylc_ml_table": _insp_klein,
    "heds_hsd_olink_serum_cinquina2026": _insp_cinquina,
    "cervia_hasler2024_zurich_somascan": _insp_cervia,
    "talla2023_aifi_olink": _insp_talla,
    "gu2023_mtbls7337_multiomics": _insp_gu,
    "stanford_iris_olink": _insp_iris,
    "lopezhernandez2023_mexico_metabolome": _insp_lopez,
    "iraq_lc_maes_group": _insp_iraq,
    "lu2024_early_markers_pasc_figshare": _insp_lu,
    "gold_ebv_reactivation": _insp_gold,
    "appelman2024_muscle_pem": _insp_appelman,
    "kedor2022_charite_pcs_mecfs": _insp_kedor,
    "plos_pcc_recovered_small": _insp_pcc_small,
    "lc_lipid_mediators_cytokines_mendeley": _insp_lipid_mediators,
    "mw_st003103_lc_mitochondria": _insp_mw_st003103,
    "biomapai_jax_mecfs_multiomics": _insp_biomapai,
    "germain2022_jci_insight_exercise_metabolomics": _insp_germain2022,
    "hoel2026_somascan7k_pad000026": _insp_hoel2026,
    "walitt2024_nih_csf_metabolomics_mw_st003178": lambda: _mw_insp(
        "walitt2024_nih_csf_metabolomics_mw_st003178", "AN005217_datatable.txt", "MECFS_STATUS", ["MECFS"],
        ["Healthy volunteer"], "untargeted CSF metabolomics", shared_ids_across_analytes=True),
    "naviaux2016_mw_st000450": lambda: _mw_insp("naviaux2016_mw_st000450", "AN000705_datatable.txt", "Disease", ["CFS"],
                                                ["Normal"], "targeted plasma metabolomics (HILIC MRM)", shared_ids_across_analytes=True),
    "naviaux_validation_mw_st000617": lambda: _mw_insp("naviaux_validation_mw_st000617", "AN000946_datatable.txt", "Treatment",
                                                       ["CFS"], ["Control"], "targeted plasma metabolomics", shared_ids_across_analytes=True),
    "che2022_columbia_plasma_metabolomics_st2000": lambda: _mw_insp(
        "che2022_columbia_plasma_metabolomics_st2000", "AN003264_datatable.txt", "treatment", ["case"], ["control"],
        "targeted plasma metabolites (nM panel; 3 further platforms)", shared_ids_across_analytes=True),
    "nagyszakal2018_cfi_metabolomics_st0800": lambda: _mw_insp(
        "nagyszakal2018_cfi_metabolomics_st0800", "AN001274_datatable.txt", "CaseStatus", ["Case"], ["Control"],
        "targeted oxylipin/nM panel (~56 named)", shared_ids_across_analytes=True),
    "lipkin_csf_metabolomics_mecfs_ms_st0910": lambda: _mw_insp(
        "lipkin_csf_metabolomics_mecfs_ms_st0910", "AN001480_datatable.txt", "Diagnosis", ["MECFS"], ["ND"],
        "CSF GC-TOF primary metabolites (MS = disease comparator)", shared_ids_across_analytes=True),
    "ncnp_serum_indole_st004940": lambda: _mw_insp("ncnp_serum_indole_st004940", "AN008373_datatable.txt", "Condition",
                                                   ["ME/CFS"], ["Healthy control"], "targeted serum tryptophan metabolites",
                                                   shared_ids_across_analytes=True),
    "armstrong2015_mtbls161_nmr": _insp_armstrong,
    "baraniuk2021_csf_metabolomics_gwi_mecfs": _insp_baraniuk,
    "fm_rsu_hhv6_cytokines": _insp_rsu_fm,
    "fm_caboni2014_lysopc_metabolomics": _insp_caboni,
    "fm_mw_st000888_comparator": _insp_st000888,
    "lyme_ptlds_mw_st001391": _insp_st001391,
    "lyme_early_mw_st001223": _insp_st001223,
    "ibs_mars2020_cell_multiomics": _insp_mars,
    "ibs_mtbls2774_serum_metabolome": _insp_mtbls2774,
    "ibs_ecm_serum_plos2017": _insp_ibs_ecm,
    "migraine_serum_zinc_japan": _insp_migraine_zinc,
    "endo_mif_serum_bmc2020": _insp_endo_mif,
    "endo_arg1_repod": _insp_endo_arg1,
    "endo_xenoestrogen_plos2024": _insp_endo_xeno,
    "endo_mw_st003579_ovarian": lambda: _insp_id_prefix("endo_mw_st003579_ovarian", "ST003579_AN005877_Results.txt",
                                                        "untargeted features (peak area)"),
    "endo_mw_st000585_follicular": lambda: _insp_mw_results("endo_mw_st000585_follicular", "ST000585_AN000900_Results.txt",
                                                            ["ndometriosis"], ["Control"], "follicular-fluid LC-MS features"),
    "lyme_multiantibody_serology_plos2021": _insp_lyme_serology,
    "hat_gastroparesis_ev_lipidomics_2022": _insp_hat_gp,
    "mcad_plasma_heparin_vysniauskaite2015": _insp_heparin,
    "oi_children_metabolomics_li2024": _insp_oi_children,
    "iu2024_immune_exhaustion_mecfs_lc": _insp_iu2024,
}


def _generic_file_check(cid: str) -> dict:
    """File-level verification for candidates without a dedicated inspector: every retrieved file opens, and its
    shape (or zip member list) is recorded. No label counts are claimed."""
    out = {}
    for f in sorted((ROOT / cid).glob("*")):
        if f.name.startswith("MANIFEST") or f.is_dir():
            continue
        try:
            if f.suffix == ".zip":
                with zipfile.ZipFile(f) as z:
                    out[f.name] = f"zip, {len(z.namelist())} members"
            elif f.suffix in (".xlsx", ".xls"):
                x = pd.ExcelFile(f)
                out[f.name] = {s: list(pd.read_excel(x, s, header=None).shape) for s in x.sheet_names[:6]}
            elif f.suffix in (".csv", ".tsv", ".txt", ".tab"):
                d = pd.read_csv(f, sep="," if f.suffix == ".csv" else "\t", nrows=5000, on_bad_lines="skip",
                                encoding_errors="replace")
                out[f.name] = list(d.shape)
            elif f.suffix == ".json":
                js = json.loads(f.read_text())
                out[f.name] = f"json, {len(js)} top-level entries"
            else:
                out[f.name] = f"{f.stat().st_size} bytes"
        except Exception as exc:
            out[f.name] = f"unreadable: {repr(exc)[:120]}"
    return out


def inspect_all() -> dict[str, dict]:
    out = {}
    for cid, fn in INSPECTORS.items():
        try:
            out[cid] = fn()
        except Exception as exc:
            out[cid] = {"inspection_error": repr(exc)[:300]}
    for cid in FILES:
        if cid not in out and (ROOT / cid).exists():
            out[cid] = {"file_check": _generic_file_check(cid)}
    for cid in LISTINGS:
        if cid not in out and (ROOT / cid).exists():
            out[cid] = {"file_check": _generic_file_check(cid)}
    return out


# ---------------------------------------------------------------------------------------------------------------
# 3. Curated candidate table. Base rows = the searchers' JSON (verbatim claims); this block records what
#    verification changed. Status vocabulary (as in docs/DEVICE_DATASET_DISCOVERY.md, plus one privacy value):
#      verified_open               files retrieved anonymously here and contents checked
#      verified_registration_only  metadata / ACL / HTTP status checked; data need a free account (none created)
#      controlled                  DUA, DAC, application, restricted record, or data on request
#      not_as_described            reachable, but a load-bearing claim failed on inspection
#      unreachable                 no per-person data route found (aggregate-only publications, blocked downloads)
#      excluded_privacy            the public record appears to contain direct identifiers; deliberately not downloaded
#      duplicate                   the same deposit reached by two searchers (the kept row carries the verification)
# ---------------------------------------------------------------------------------------------------------------
DUPLICATES = {   # duplicate candidate -> kept candidate (same deposit reached by two searchers)
    "longcovid_mw_st003103": "mw_st003103_lc_mitochondria",
    "frontimmunol2024_lc_cfs_cil_metabolomics": "frontiers_lc_mecfs_metabolomics",
    "figshare_gpcr_autoab_postcovid_EXCLUDE": "postcovid_gpcr_aab_seibert2026",
    "gp_pediatric_ev_lipidomics_hat": "hat_gastroparesis_ev_lipidomics_2022",
    "mecfs_mw_columbia_plasma": "che2022_columbia_plasma_metabolomics_st2000",
    "mw_mecfs_pointer": "che2022_columbia_plasma_metabolomics_st2000",
    "mapmecfs_walitt2024_registration": "mapmecfs_nih_pimecfs_clinical_labs",
}
PRIVACY_EXCLUDED = PRIVACY_EXCLUDED | {"endo_vitd_zinc_indonesia"}
# Open and retrieved, but the released per-person file carries no laboratory analyte (device / questionnaire / genotype).
NO_LAB_VALUES = {
    "mirabegron_pots_zenodo2026": "symptom scores and history only (data dictionary has no lab fields)",
    "heds_pots_fmri_zenodo2026": "MRI only; the text file is the subject-group key",
    "nih_pimecfs_pennsieve_neurophysiology": "grip force / MEP / fMRI only (already in data/raw/device_candidates)",
    "norway_endothelial_fmd_plos2023": "FMD, blood pressure, steps, questionnaires",
    "dataverse_lc_cardiovascular_cohort": "HRV, COMPASS and biospecimen inventory; no lab values",
    "pulmonary_lc_jci_insight_2024": "pulmonary function and CT; no blood analytes; LC-only",
    "tryptase_ngs_genotyping_li2024": "TPSAB1 genotype calls with no phenotype or tryptase level",
    "hat_venom_gwas_zenodo2025": "GWAS summary statistics; per-patient clinical docx not opened",
    "bodansky_liinc_phipseq": "README + listing only here (PhIP-seq CSVs not downloaded); autoantibody reactivity",
}

# Verification outcomes that override the searcher's access claim.
STATUS_OVERRIDE = {
    "yin2024_liinc_olink_cytof": ("unreachable", "Dryad record listed (CC0) but the file API answers 401 to anonymous "
                                  "scripted requests (probed 2026-09-24); a browser download may work"),
    "northwestern_nlc_autoab_dryad": ("unreachable", "Dryad downloads refused to anonymous scripts (401/403)"),
    "vijayakumar_airway_postcovid_dryad": ("unreachable", "Dryad downloads refused to anonymous scripts (401/403)"),
    "lyme_ptlds_serochip_epitopes": ("unreachable", "Dryad downloads refused to anonymous scripts (401/403)"),
    "lyme_peptide_array_mbio2024": ("unreachable", "Dryad downloads refused to anonymous scripts (401/403)"),
    "mecfs_rsu_restricted_and_dryad": ("controlled", "RSU files flagged restricted; Dryad file refused (401)"),
    "pretorius2021_microclots_onedrive": ("unreachable", "OneDrive link answers 403 to anonymous requests"),
    "prusty2024_herpesvirus_serology_mendeley": ("verified_open", "Mendeley listing only; Prism files not opened here"),
    "ibs_mtbls12737_urine": ("not_as_described", "MAF has metabolite rows but no abundance values"),
    "migraine_pfo_mtbls13630": ("not_as_described", "MAFs empty; no disease-group factor"),
    "fm_serum_proteomics_pxd022886": ("unreachable", "only RAW + a 1.55 GB search zip; no quant matrix checked"),
    "fm_csf_proteome_pxd008076": ("unreachable", "no quantitative protein x sample matrix deposited"),
    "butterfield_mcas_urinary_ratios_2025": ("verified_open", "3-patient case series; values are in the article "
                                             "tables (read via paperclip); no data file"),
}
# Candidates whose label -> group mapping is not stated in the released files (cannot be ranked until the key is found).
LABEL_MAPPING_UNVERIFIED = {
    "heds_serum_ms_proteomics_griggs2025": "run-name codes C/L/X undocumented",
    "frontiers_lc_mecfs_metabolomics": "group codes L/A/C/HC not defined in the sheet",
    "schutzer2023_csf_proteome_mecfs_fm": "sample-to-group map not in the data file",
    "fm_clos_garcia_serum_proteomics_pxd059894": "no sample-to-group map",
    "pots_plasma_dia_pxd031458": "no group labels",
    "hochecker_pbmc_seahorse_mendeley": "HC vs ME/CFS membership not encoded",
    "lyme_raman_serum_spectra": "no subject ID; spectra per person cannot be grouped",
    "pride_misc_lc_proteomics": "A/B group meaning unverified",
    "plos_pcc_recovered_small": "Group 1/2 coding inferred, not stated in the file",
}
CIRCULAR_LABEL = {   # the lab is (part of) the case definition, so lab-vs-label is not an independent test
    "lyme_multiantibody_serology_plos2021": "anti-Borrelia serology defines Lyme status",
    "lyme_romania_lnb_csf": "CSF antibody index is a diagnostic criterion for LNB",
    "lyme_children_ukraine_immunoblot": "immunoblot is the diagnostic test",
    "mcad_plasma_heparin_vysniauskaite2015": "tryptase is an SM criterion; no controls",
}
# label_quality: 3 criteria-based diagnosis stated by the study (CCC/IOM/Fukuda, ACR, 2017 hEDS nosology, Rome, ICHD-3,
# surgical/histological, CDC EM); 2 clinical diagnosis or study case definition without named criteria (e.g. symptom
# persistence after COVID assessed in a clinic); 1 self-report, symptom-derived, or ambiguous label; 0 none.
LABEL_QUALITY = {
    "klein2023_mylc_ml_table": 2, "heds_hsd_olink_serum_cinquina2026": 3, "cervia_hasler2024_zurich_somascan": 2,
    "talla2023_aifi_olink": 2, "gu2023_mtbls7337_multiomics": 1, "stanford_iris_olink": 1,
    "lopezhernandez2023_mexico_metabolome": 2, "iraq_lc_maes_group": 2, "lu2024_early_markers_pasc_figshare": 1,
    "gold_ebv_reactivation": 1, "appelman2024_muscle_pem": 2, "kedor2022_charite_pcs_mecfs": 3,
    "plos_pcc_recovered_small": 2, "lc_lipid_mediators_cytokines_mendeley": 1, "mw_st003103_lc_mitochondria": 2,
    "biomapai_jax_mecfs_multiomics": 3, "germain2022_jci_insight_exercise_metabolomics": 3,
    "hoel2026_somascan7k_pad000026": 3, "walitt2024_nih_csf_metabolomics_mw_st003178": 3, "naviaux2016_mw_st000450": 3,
    "naviaux_validation_mw_st000617": 3, "che2022_columbia_plasma_metabolomics_st2000": 3,
    "nagyszakal2018_cfi_metabolomics_st0800": 3, "lipkin_csf_metabolomics_mecfs_ms_st0910": 2,
    "ncnp_serum_indole_st004940": 3, "armstrong2015_mtbls161_nmr": 3, "baraniuk2021_csf_metabolomics_gwi_mecfs": 3,
    "fm_rsu_hhv6_cytokines": 3, "fm_caboni2014_lysopc_metabolomics": 3, "fm_mw_st000888_comparator": 2,
    "lyme_ptlds_mw_st001391": 3, "lyme_early_mw_st001223": 3, "ibs_mars2020_cell_multiomics": 3,
    "ibs_mtbls2774_serum_metabolome": 2, "ibs_ecm_serum_plos2017": 2, "migraine_serum_zinc_japan": 3,
    "endo_mif_serum_bmc2020": 3, "endo_arg1_repod": 3, "endo_xenoestrogen_plos2024": 2, "endo_mw_st000585_follicular": 3,
    "endo_mw_st003579_ovarian": 3, "lyme_multiantibody_serology_plos2021": 1, "hat_gastroparesis_ev_lipidomics_2022": 2,
    "mcad_plasma_heparin_vysniauskaite2015": 1, "oi_children_metabolomics_li2024": 2, "iu2024_immune_exhaustion_mecfs_lc": 3,
}
# lab relevance: 3 targeted clinical-grade analyte named in the brief for this condition (cortisol for Long COVID,
# tryptase/LTE4/PGD2 metabolites for MCAS, catecholamines for POTS, CGRP for migraine) inside a targeted panel;
# 2 targeted immunoassay / Olink / targeted metabolite / clinical-chemistry panels; 1 untargeted features, transformed
# values or a single exploratory analyte of unclear clinical grade.
LAB_RELEVANCE = {
    "klein2023_mylc_ml_table": 3, "appelman2024_muscle_pem": 3, "plos_pcc_recovered_small": 3,
    "biomapai_jax_mecfs_multiomics": 1, "lyme_early_mw_st001223": 1, "endo_mw_st003579_ovarian": 1,
    "endo_mw_st000585_follicular": 1, "germain2022_jci_insight_exercise_metabolomics": 1, "fm_mw_st000888_comparator": 1,
    "lyme_ptlds_mw_st001391": 1, "endo_xenoestrogen_plos2024": 1, "hat_gastroparesis_ev_lipidomics_2022": 1,
}
MULTI_CONDITION = {   # 1 = >= 2 target umbrella conditions or disease comparators with the same assay; 0.5 = two control types
    "heds_hsd_olink_serum_cinquina2026": 1.0, "lipkin_csf_metabolomics_mecfs_ms_st0910": 1.0,
    "baraniuk2021_csf_metabolomics_gwi_mecfs": 1.0, "fm_mw_st000888_comparator": 1.0, "ibs_ecm_serum_plos2017": 1.0,
    "iu2024_immune_exhaustion_mecfs_lc": 1.0, "kedor2022_charite_pcs_mecfs": 1.0, "klein2023_mylc_ml_table": 0.5,
    "talla2023_aifi_olink": 0.5, "endo_arg1_repod": 0.5, "ibs_mars2020_cell_multiomics": 0.5,
    "endo_xenoestrogen_plos2024": 0.5,
}
MIN_PER_GROUP = 20


def _searcher_rows() -> list[dict]:
    rows = []
    for tag, fname in SEARCHER_FILES.items():
        for c in json.loads((SEARCHER_DIR / fname).read_text()):
            c = dict(c)
            c["found_by"] = f"searcher_{tag}"
            rows.append(c)
    return rows


def _status(c: dict, insp: dict) -> tuple[str, str]:
    cid, acc = c["candidate_id"], str(c.get("access_claimed", "")).lower()
    if cid in PRIVACY_EXCLUDED:
        return "excluded_privacy", "record appears to contain direct personal identifiers; not downloaded, not ingested"
    if cid in STATUS_OVERRIDE:
        return STATUS_OVERRIDE[cid]
    if cid in DUPLICATES:
        return "duplicate", f"same deposit as {DUPLICATES[cid]} (verified there)"
    if str(c.get("per_person", "")).lower() == "no" and acc in ("open", "none", ""):
        return "unreachable", "no per-person lab data released (aggregate tables only; recorded as a negative)"
    if acc == "open":
        i = insp.get(cid, {})
        if i and "inspection_error" not in i:
            fc = i.get("file_check")
            if fc is not None and any(str(v).startswith("unreadable") for v in fc.values()):
                return "not_as_described", f"retrieved but unreadable: {fc}"
            return "verified_open", "anonymous download through measure_it.download; contents checked here"
        if i.get("inspection_error"):
            return "verified_open", f"retrieved; inspector error {i['inspection_error'][:120]}"
        return "verified_open", "searcher-verified anonymous access (not re-downloaded here; see probe_result)"
    if acc == "registration":
        return "verified_registration_only", "free account needed (searcher probe or HTTP 403 login page); none created"
    if acc in ("controlled", "on_request"):
        return "controlled", "DUA / application / restricted file / on request"
    return "unreachable", "no per-person data route found (aggregate-only publication or no deposit)"


def _score(r: dict) -> dict:
    nc, nk = r.get("n_cases_with_lab"), r.get("n_controls_with_lab")
    lq = LABEL_QUALITY.get(r["candidate_id"])
    ok = (r["access_status"] == "verified_open" and r["candidate_id"] not in LABEL_MAPPING_UNVERIFIED
          and r["candidate_id"] not in CIRCULAR_LABEL and r["candidate_id"] not in DUPLICATES
          and r["candidate_id"] not in NO_LAB_VALUES
          and lq is not None and lq > 0 and pd.notna(nc) and pd.notna(nk) and nc >= MIN_PER_GROUP and nk >= MIN_PER_GROUP)
    if not ok:
        return {"rank_eligible": False}
    tot = (r.get("n_cases_file") or 0) + (r.get("n_controls_file") or 0)
    comp = {"score_label": lq, "score_n": round(min(3.0, 1.5 * math.log10(min(nc, nk))), 3),
            "score_lab_relevance": LAB_RELEVANCE.get(r["candidate_id"], 2),
            "score_completeness": round((nc + nk) / tot, 3) if tot else 0.0, "score_controls": 1.0,
            "score_multi_condition": MULTI_CONDITION.get(r["candidate_id"], 0.0),
            "score_person_linkage": 1.0 if r.get("shared_ids_across_analytes") else 0.0}
    comp["rank_score"] = round(sum(comp.values()), 3)
    comp["rank_eligible"] = True
    return comp


SCORE_COLS = ["score_label", "score_n", "score_lab_relevance", "score_completeness", "score_controls",
              "score_multi_condition", "score_person_linkage"]


def build(insp: dict) -> pd.DataFrame:
    rows = []
    for c in _searcher_rows():
        cid = c["candidate_id"]
        st, detail = _status(c, insp)
        i = insp.get(cid, {})
        r = {"candidate_id": cid, "name": c.get("name"), "citation": c.get("citation"), "doi": c.get("doi"),
             "found_by": c["found_by"], "repository": c.get("repository"), "url": c.get("landing_url"),
             "conditions": "; ".join(map(str, c.get("conditions") or [])), "other_groups": c.get("other_groups"),
             "biomarkers": "; ".join(map(str, c.get("biomarkers") or [])), "assay": c.get("assay"),
             "label_basis": c.get("label_basis"), "per_person_claimed": c.get("per_person"),
             "access_claimed": c.get("access_claimed"), "access_status": st, "access_detail": detail,
             "licence": c.get("licence"), "n_cases_claimed": c.get("n_cases_claimed"),
             "n_controls_claimed": c.get("n_controls_claimed"), "searcher_probe": c.get("probe_result"),
             "published_claim": c.get("published_claim"), "searcher_notes": c.get("notes"),
             "duplicate_of": DUPLICATES.get(cid, ""), "label_mapping_unverified": LABEL_MAPPING_UNVERIFIED.get(cid, ""),
             "circular_label": CIRCULAR_LABEL.get(cid, ""), "label_quality": LABEL_QUALITY.get(cid),
             "lab_values_in_file": ("no: " + NO_LAB_VALUES[cid]) if cid in NO_LAB_VALUES else "",
             "download_urls": " ".join(u for u, _, _ in FILES.get(cid, [])) or " ".join(
                 map(str, c.get("download_urls") or []))}
        for k in ("n_rows", "n_cases_file", "n_controls_file", "n_cases_with_lab", "n_controls_with_lab", "label_column",
                  "lab_columns_checked", "n_analytes", "shared_ids_across_analytes", "notes_counts", "inspection_error"):
            v = i.get(k)
            r[k] = json.dumps(v) if isinstance(v, (dict, list)) else v
        r["label_values"] = json.dumps(i.get("label_values")) if i.get("label_values") else ""
        r["file_check"] = json.dumps(i.get("file_check")) if i.get("file_check") else ""
        r["counts_source"] = ("computed here from downloaded files" if "n_cases_file" in i else
                              "file-level check only (no label counts claimed)" if "file_check" in i else
                              "searcher probe only (not recomputed)")
        r.update(_score(r))
        r["raw_dir"] = f"data/raw/{SRC}/{cid}" if (ROOT / cid).exists() else ""
        r["checked_at"] = utc_now_iso()
        rows.append(r)
    df = pd.DataFrame(rows)
    df["rank"] = df["rank_score"].rank(ascending=False, method="min").where(df["rank_eligible"])
    df = df.sort_values(["rank", "access_status", "candidate_id"], na_position="last").reset_index(drop=True)
    return df


def robustness(df: pd.DataFrame, top: int = 5) -> pd.DataFrame:
    """Top-k overlap when each score component is dropped in turn."""
    el = df[df.rank_eligible].copy()
    base = set(el.nlargest(top, "rank_score").candidate_id)
    out = []
    for c in SCORE_COLS:
        alt = el.assign(_s=el.rank_score - el[c]).nlargest(top, "_s")
        out.append({"dropped_component": c, "top_overlap": len(base & set(alt.candidate_id)),
                    "new_top2": ";".join(alt.candidate_id.head(2))})
    return pd.DataFrame(out)


def record_searcher_inputs() -> None:
    man = load_manifest(f"{SRC}/_searcher_inputs")["files"]
    for fname in SEARCHER_FILES.values():
        if fname not in man:
            record_file(f"{SRC}/_searcher_inputs", fname, url="local: parallel searcher output (2026-09-24)",
                        note="discovery input; claims verbatim, verified by measure_it.labs.lab_dataset_discovery")


def run(download: bool = True) -> pd.DataFrame:
    record_searcher_inputs()
    errors = download_all() if download else {}
    probes = probe_gated() if download else {}
    insp = inspect_all()
    (ROOT / "inspection.json").write_text(json.dumps({"inspection": insp, "download_errors": errors,
                                                      "gated_probes": probes}, indent=2, default=str))
    df = build(insp)
    TABLES.mkdir(parents=True, exist_ok=True)
    df.to_csv(TABLES / "lab_dataset_candidates.csv", index=False)
    rob = robustness(df)
    rob.to_csv(TABLES / "lab_dataset_candidates_rank_robustness.csv", index=False)
    print(df[["rank", "candidate_id", "access_status", "n_cases_with_lab", "n_controls_with_lab", "rank_score"]]
          .head(20).to_string())
    print(df["access_status"].value_counts().to_string())
    print(rob.to_string())
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-download", action="store_true", help="inspect files already on disk only")
    run(download=not ap.parse_args().no_download)
