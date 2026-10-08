# GitHub release plan

Prepared 2026-09-28 for publishing this repository on GitHub.

**Update 2026-10-08.** D1 is decided: option (c), one fresh squashed commit on branch `main` (section 3.3); the
development history stays on the owner's machine (`master`, `maestro-stages-b-f`). D5 is done (the two Stanford
per-person tables and the Cinquina tables are git-ignored) and D7 is done except the registry licence strings in the source modules (`docs/DATA_LICENSES.md` section 4 has the
correct ones; e.g. Cinquina still reads CC BY). Still open: D2 (author e-mail
of the squashed commit), D3 (contacts) and D6 (name). No remote exists and nothing has been pushed. The audit figures
below (12 commits, 808 unit tests) are as of 2026-09-28; at the squash CI runs 1,007 unit tests.
D4 was revisited the same day: git now keeps 410 files (code, docs, write-ups, embedded figures, provenance
digests, hand-made records); the regenerated tables, maps, run logs, `.meta.json` sidecars and per-download sidecars
(902 files) moved to the release asset `measure-it-reference-results-<date>.zip` (`results/README.md`).

## 1. Decisions

### Made by the owner (2026-09-28)

| decision | choice | consequence |
|---|---|---|
| visibility | **private** repository on the owner's personal account; made public later, after review | GitHub Free private repos: 2,000 Actions minutes/month (CI uses about 2-3 min per push); branch-protection rules and private vulnerability reporting need a paid plan or a public repo |
| code licence | **MIT**, copyright 2026 <OWNER NAME> (`LICENSE`) | covers code and project-written docs only; data and derived tables stay under their sources' terms (README "Licence", `docs/DATA_LICENSES.md`) |
| data distribution | **code only; no data bundle release asset** | users build the data themselves with `uv run measure-it pipeline --workers 16` (network, about 40 GB of downloads, about 45 GB of disk, untimed and never run end to end from an empty checkout; README "Building the data"). `docs/DATA_LICENSES.md` records what a future bundle could contain |

Code-licence options that were considered: MIT (shortest, maximally permissive, no patent grant), Apache-2.0
(permissive plus an explicit patent licence and NOTICE handling; longer), BSD-3-Clause (MIT-like plus a
no-endorsement clause).

### Still open

| # | decision | options | recommendation |
|---|---|---|---|
| D1 | **git history** (section 3) | (a) push the 12 commits as they are; (b) rewrite them with `git filter-repo`; (c) publish one fresh squashed commit | (b) or (c) **before the first push**. The history holds expired pre-signed download URLs (provider AWS key ids, signatures, a GitHub release JWT, a Pennsieve STS token), absolute paths with the local user name, and the personal e-mail as author. None is a live credential, but rewriting after the repository has been public is much harder than before the first push |
| D2 | **commit author e-mail** | keep `<personal e-mail>`; or rewrite to the GitHub no-reply address (`<id>+<user>@users.noreply.github.com`, shown under GitHub Settings > Emails) | decide before D1; with (a) the personal address becomes public when the repo does |
| D3 | **contacts** in `SECURITY.md` (`SECURITY_CONTACT_EMAIL`) and `CODE_OF_CONDUCT.md` (`CONDUCT_CONTACT_EMAIL`), and `OWNER/REPO` in `.github/ISSUE_TEMPLATE/config.yml` | a dedicated address, the no-reply address plus GitHub private vulnerability reporting (public repos), or a personal address | fill in before going public |
| D4 | **results in git** (`results/tables` 85 MB, `results/figures` 42 MB, `results/maps` 17 MB in HEAD) | keep (the README and FINAL_REPORT cite them); drop `results/maps` and the two largest CSVs (`test3_reactome_ora.csv` 21 MB, `test3_gene_sets.csv` 12 MB) and let the pipeline regenerate them | keep: every file is under 25 MB and the whole push is about 85 MB |
| D5 | **four per-person derived tables in `results/tables`** (section 4.3) | keep all four; or remove the two built from the Stanford Mishra 2020 / Alavi 2022 downloads | **remove `stanford_acute_participants.csv` and `stanford_acute_detection_windows.csv` before going public** (the downloads carry no licence and the Mishra article is not open access; `docs/DATA_LICENSES.md`), or get the Snyder lab's permission; the Uwakwe (ODC-By) and fibromyalgia thermography (CC BY) tables can stay with attribution. The pipeline regenerates them locally |
| D5b | **Cinquina 2026 (hEDS/HSD Olink) is CC BY-NC-ND 4.0**, not CC BY as the registry says | keep the six group-level `results/tables/heds_hsd_olink_serum_cinquina2026_*` tables and the write-up (statistics computed from the file; whether they are "adapted material" is uncertain); or drop the per-analyte table (`..._protein_effects.csv`, 458 analytes with group medians) and keep only the headline statistics | decide before going public; the private phase is unaffected |
| D6 | repository name | `measure-it-public` (package name) or another | any; the README uses relative links only |
| D7 | small follow-ups outside this change set | set `license = "MIT"` in `pyproject.toml`; change the Dockerfile label `org.opencontainers.image.licenses` from "see README.md" to "MIT"; make `download.download_file` strip signed query strings from `final_url` (as `ingestion/mapmecfs.py:_strip_query` already does) so a cold build does not write new signed URLs into `MANIFEST.json`; correct the registry licence strings listed in `docs/DATA_LICENSES.md` section 4 (Cinquina CC BY-NC-ND, Reactome CC0, RepOD CC BY 4.0, NUCC no-derivatives) in the modules that write them | do them; CI's hygiene job fails if a signed URL is committed again |

## 2. Audit findings (2026-09-28)

### 2.1 Secrets and privacy

Scanned: every tracked file at HEAD and every blob reachable from any commit (2,234 blobs), for cloud keys, tokens,
JWTs, private keys, passwords, e-mail addresses, the local user and host names, absolute local paths, the two
privacy-excluded records, and per-person rows.

| finding | where | status |
|---|---|---|
| pre-signed download URLs (AWS `X-Amz-Credential` with provider key ids such as `AKIAJ45O...`, `AKIAZT3G...`, `AKIA3OGA...`; Google `X-Goog-Credential` / `GoogleAccessId` + `Signature`; GitHub release-asset `sig` + JWT; Pennsieve `X-Amz-Security-Token`) | `final_url` (and in one case `url`) of 54 `data/raw/**/MANIFEST.json` files, e.g. `data/raw/fm_thermography/MANIFEST.json:7`, `data/raw/mondo_ontology/MANIFEST.json:17`, `data/raw/device_candidates/pennsieve_356_mecfs_emg/MANIFEST.json:7,11` | **fixed in the working tree** (query strings stripped; the path part, sha256 and original `url` kept). Still in history: introduced or changed in commits `2478d0a`, `5a0de28`, `209f666` and `4bb55f4`, and present in every later commit's tree. All are the providers' short-lived download signatures, expired days after retrieval; none is a key of this project |
| absolute path `<local checkout path>/...` (reveals the local user name) | `results/pipeline_runs/20260924T215649Z_vs_pre_integration.csv:172`, `results/pipeline_runs/20260924_integration_vs_pre_integration.csv:172` | **fixed** (made relative). In history also in `results/pipeline_runs/20260924T0*.jsonl` and `20260924T1*.jsonl` ("log" fields; already relative at HEAD) |
| author and committer `<OWNER NAME> <<personal e-mail>>` | metadata of all 12 commits; `.git/config` `user.email` | decision D2 |
| e-mail addresses in content | `mpsnyder@stanford.edu` (`docs/DATA_ACCESS_PLAN.md`), `mapMECFS@rti.org` (`src/measure_it/measurements/device_dataset_discovery.py`, `results/tables/device_dataset_candidates.csv`), `phosp@leicester.ac.uk` (`results/tables/lab_dataset_candidates.csv`), `x@example.org` (test), `feross@feross.org` (licence header inside `results/maps/plotly.min.js`) | public institutional or library addresses; left as they are |
| API keys, tokens, passwords, private keys | none found (the only `AKIA...` matches are the provider ids inside the signed URLs above) | none |
| host name (`<host name>`), user name outside the paths above | none | none |
| privacy-excluded records figshare 30127723 and Mendeley vzydrgtymt | ids, titles, column *headers* and non-null counts only: `docs/LAB_DATASET_DISCOVERY.md:104-107`, `results/FINAL_REPORT.md:1414-1416`, `results/tables/lab_dataset_candidates.csv:64-66`, `src/measure_it/labs/lab_dataset_discovery.py:22` | acceptable (no values; no file from either record was ever committed) |
| per-person rows from gated sources | none: mapMECFS login-gated files were never downloaded; no NHANES per-person rows in tracked files (SEQN ids appear only as examples in docs) | none |
| data payloads ever committed | only the Census geography tables `data/processed/{county,state}_boundaries_2024.geoparquet`, `geographies.parquet`, `zcta_centroids.parquet`, `zcta_county_crosswalk.parquet` in early commits (public domain, removed since) | harmless |
| reflog-only commit `fd0e41a` (the pre-amend version of "Ingest 25 public sources") with processed Parquet tables, including Stanford per-participant daily wearable data (largest blob 132 MB) | local `.git` only; not reachable from any branch, so **it is not pushed** | clean locally (section 3.4) |

### 2.2 Size

* The local pack is 264 MB, but 271 MB of it are the reflog-only objects above. A fresh clone of the branch packs to
  **85 MB**; a single squashed commit would be about **55 MB**.
* Largest blobs anywhere in history: `results/tables/test3_reactome_ora.csv` 21.3 MB (4 versions),
  `results/tables/test3_gene_sets.csv` 12.0 MB, `results/tables/linked_omics_univariate.csv` 8.9 MB,
  `results/tables/labs_dx_nhanes_effects.csv` 5.7 MB, `results/maps/fig07_infrastructure_overlay.html` 5.7 MB,
  `results/figures/static_snapshot_desktop.png` 5.0 MB, `results/maps/plotly.min.js` 4.8 MB.
* **No file or blob exceeds 50 MB** (GitHub's warning) **or 100 MB** (GitHub's limit); none is needed in Git LFS. CI
  fails if a tracked file over 50 MB is added.

### 2.3 Tests on a clean clone

`git clone` of HEAD into an empty directory (no `data/raw` payloads, no `data/processed` tables), `uv sync --frozen`,
`uv run pytest -m "not network and not data"`: **808 passed, 620 deselected** in 16 s, after one marker fix
(`tests/test_fm_thermography.py::test_condition_id_matches_ontology_normalizer` reads
`data/processed/condition_registry.parquet` and was unmarked; now `@pytest.mark.data`). The lint subset
(`ruff --select E9,F63,F7,F82`) passes. `docker build` of the code-only target succeeds and the API answers `/health`
with `degraded` (no data), as intended. The lock resolves to wheels on Linux x86_64 (what GitHub runners use) and
macOS arm64 14+.

## 3. History options (D1)

Work on a copy; keep the original until the pushed result has been checked.

### 3.1 Option (a): push as is

No action. The signed URLs are expired and belong to the data providers. GitHub push protection may still flag the
AWS key ids / session token; if a push is blocked, use (b) or (c) rather than bypassing.

### 3.2 Option (b): rewrite with git filter-repo

```bash
uv tool install git-filter-repo                        # or: pipx install git-filter-repo
git clone --no-local ~/measure-it-public ~/measure-it-public-rewrite
cd ~/measure-it-public-rewrite

cat > ../mi-replacements.txt <<'EOF'
<local checkout path>/==>
regex:("(?:final_url|url)": "[^"?]+)\?[^"]*?(?:X-Amz-|X-Goog-|sig=|Signature=|GoogleAccess[I]d=)[^"]*"==>\1"
EOF
# optional (D2): map the author e-mail to the GitHub no-reply address
cat > ../mi-mailmap.txt <<'EOF'
<OWNER NAME> <NOREPLY_ID+USER@users.noreply.github.com> <<personal e-mail>>
EOF

git filter-repo --replace-text ../mi-replacements.txt --mailmap ../mi-mailmap.txt

# check: nothing left, same tree as the cleaned working tree
git log --all -p | grep -c -E 'X-Amz-Credential|X-Goog-Credential|GoogleAccessId|/home/<user>' # expect 0
git log --format='%an <%ae> | %cn <%ce>' | sort -u
```

The regex is the one used to clean HEAD: applied to the 54 manifests at HEAD it reproduces the cleaned working tree
exactly. Commit hashes change; the commit messages are kept (they contain no paths or secrets; they do contain
`Co-Authored-By: Claude ... <noreply@anthropic.com>` lines).

### 3.3 Option (c): one fresh commit

```bash
cd ~/measure-it-public
git checkout --orphan public-main                     # after committing the publication changes on master
git commit -m "measure-it-public 0.1.0: public-data invisible-illness measurement engine"
git branch -M public-main main                        # push this branch; keep master locally as the full history
```

Simplest and smallest (about 55 MB), and it drops every historical blob; the development history stays only on the
owner's machine. Set the author first if D2 says so (`git config user.email ...`).

### 3.4 Local cleanup (any option)

The 132 MB processed-table blobs of the amended commit `fd0e41a` are kept alive only by the reflog. They are never
pushed, but to drop them locally:

```bash
git reflog expire --expire=now --all && git gc --prune=now
```

## 4. Data redistribution summary (for later)

See `docs/DATA_LICENSES.md`. For this code-only release:

1. **Nothing in the repository is a raw data file.** `data/` holds only provenance sidecars (DATA_AUDIT.md,
   MANIFEST.json, registry fragments, `.meta.json`), enforced by CI.
2. **`results/` holds derived aggregate tables, figures and write-ups** of the reference build: model metrics,
   rankings of counties and facilities (from public-domain federal data), molecular enrichment (Reactome CC BY 4.0,
   Open Targets CC0, GWAS Catalog), and short quotes from published articles. These are covered by the attribution in
   `docs/DATA_LICENSES.md`.
3. **Four tracked tables are per person** (participant id + derived values, no demographics or locations):
   `results/tables/stanford_uwakwe_oof_predictions.csv` (2,268 rows; Uwakwe 2025, ODC-By 1.0),
   `results/tables/fm_thermography_oof_predictions.csv` (177; PLOS ONE CC BY 4.0),
   `results/tables/stanford_acute_participants.csv` (107) and `results/tables/stanford_acute_detection_windows.csv`
   (314; Mishra 2020 / Alavi 2022, no stated data licence) — decision D5.
4. **Licence surprises found by the audit** (full table in `docs/DATA_LICENSES.md`): Cinquina 2026 is CC BY-NC-ND
   (decision D5b); the mapMECFS data-use agreement forbids redistribution (nothing from mapMECFS itself was
   downloaded; the Walitt tables come from the CC BY article, Pennsieve and GEO); HPO and the NUCC taxonomy forbid
   altering their content; Reactome data are CC0. None of this affects the private, code-only phase.

## 5. Publishing steps (private repository, code only)

1. **Review and commit the publication changes** (this change set plus the bring-your-own-data work):

   ```bash
   cd ~/measure-it-public
   git status                                          # LICENSE, CONTRIBUTING.md, CODE_OF_CONDUCT.md, SECURITY.md,
                                                       # .github/, docs/DATA_LICENSES.md, docs/GITHUB_RELEASE_PLAN.md,
                                                       # README.md, 54 MANIFEST.json, 2 run CSVs, 1 test marker, BYOD files
   uv run pytest -m "not network and not data"         # 808+ passed
   git add -A && git commit -m "Prepare for GitHub: MIT licence, CI, contributing/security docs, data licences; strip signed URLs"
   ```

2. **History** (D1): run section 3.2 or 3.3 now if chosen; push from the rewritten copy.
3. **Create the repository on github.com**: New repository > owner = personal account > name (D6) > **Private** >
   do *not* add a README, .gitignore or licence (the repository must be empty).
4. **Push**:

   ```bash
   git branch -M master main                           # optional: GitHub's default branch name
   git remote add origin git@github.com:<account>/<repo>.git     # or https://github.com/<account>/<repo>.git
   git push -u origin main
   ```

5. **On github.com**: Actions tab, check that "CI" (and "Docker image", when its paths changed) pass. Settings >
   General: disable Wiki/Projects if unused; Settings > Code security: enable Dependabot alerts and secret scanning
   where the plan allows.
6. **Optional tag** (no assets): `git tag -a v0.1.0 -m "Alpha, reference build 2026-09-25" && git push origin v0.1.0`;
   then Releases > Draft a new release from `v0.1.0` with notes that say it is code only and point to README "Building
   the data". GitHub adds the source archives automatically; do not upload a data bundle.
7. **Invite reviewers** (Settings > Collaborators) for the pre-public review.

### Going public later

1. Settle D3 (contacts), D5 (per-person tables) and D5b (Cinquina tables); replace `OWNER/REPO` in `.github/ISSUE_TEMPLATE/config.yml`.
2. Re-run the scans of section 2.1 on the then-current history (`git log --all -p | grep -E ...`), and the CI hygiene
   job.
3. Settings > General > Danger zone > Change visibility > Public. Then enable private vulnerability reporting
   (Settings > Code security) and, optionally, a branch ruleset on `main` requiring the CI checks.
4. Add a citation (a `CITATION.cff`, or a Zenodo DOI via the GitHub-Zenodo integration on the next release) and
   update README "Citation".

## 6. Pre-publish checklist

- [ ] D1 history decision executed and checked (no signed URLs, no `/home/<user>` paths in `git log --all -p`).
- [ ] D2 author e-mail as intended (`git log --format='%an <%ae>' | sort -u`).
- [ ] `git status` clean; the bring-your-own-data files are committed and linked docs exist
      (`docs/WALKTHROUGH.md`, `docs/BRING_YOUR_OWN_DATA.md`).
- [ ] `uv run pytest -m "not network and not data"` passes on a fresh clone (`git clone <local path> /tmp/x && cd /tmp/x && uv sync --frozen && uv run pytest -m "not network and not data"`).
- [ ] CI hygiene checks pass locally (the three `run:` blocks of `.github/workflows/ci.yml` job `hygiene`).
- [ ] No file over 50 MB: `git ls-files -z | xargs -0 stat -c '%s %n' | sort -rn | head`.
- [ ] `data/` contains only sidecars: `git ls-files data | grep -v -E '(DATA_AUDIT.md|MANIFEST.json|registry_entry.yaml|QUERY_LOG.json|queries.json|\.meta\.json|\.gitkeep)$'` prints nothing.
- [ ] `dist/` (the 959 MB bundle and the static snapshot) is not tracked (`.gitignore` covers it) and is not uploaded anywhere.
- [ ] README renders on GitHub: relative links resolve (docs/, results/, LICENSE).
- [ ] `LICENSE` present (MIT, 2026, <OWNER NAME>); D7 follow-ups done or noted.
- [ ] Before going public: D3, D5, D5b, `OWNER/REPO`, and the privacy re-scan.
