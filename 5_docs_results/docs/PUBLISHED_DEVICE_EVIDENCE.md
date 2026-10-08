# Published device evidence for invisible illnesses

**Every number on this page is a published claim — not reproduced by this project.** It is copied from the verbatim
quote stored with the row in `data/processed/published_device_evidence.parquet` (66 claims, 61 publications,
2026-09-24). The build re-fetches each abstract or open full text from Europe PMC and checks every quote verbatim:
66 of 66 quotes verified, 66 of 66 rows with their AUROC / sensitivity / specificity / stated counts found inside
the quote, 66 of 66 DOIs matching Europe PMC. Rows are cited here as `PDE-###`; the object id is
`published_evidence:PDE-###`. Module: `measure_it.measurements.published_evidence`; audit and search log:
`data/raw/published_device_evidence/DATA_AUDIT.md`.

This is a targeted curation, not a systematic review. A missing study is not evidence that none exists.

## The short answer

The question was whether any device other than nailfold capillaroscopy can be shown, from public evidence, to be a
reasonable way to uncover an invisible illness. What the literature shows, in order of strength:

1. **Orthostatic heart rate with posture (active stand / NASA lean test, and wearable HR + accelerometer) for POTS
   and POTS inside Long COVID.** This is the best-supported device use. It is principled, because POTS is *defined*
   by a device-measured heart-rate rise on standing, and partly circular for the same reason (next section). Its
   practical value is that in symptomatic Long COVID cohorts an objective orthostatic tachycardia is found in 20-31%
   (PDE-045 31% of 467; PDE-046 20% of 221), and **symptoms do not pick out who has it**: symptom distributions did
   not differ between POTS and non-POTS Long COVID (PDE-045), and 83% of lean-test-negative patients still reported
   orthostatic symptoms (PDE-046). A measurement adds information that symptoms do not. These rows do not show that
   finding the subgroup changes treatment or outcome; POTS patients did have lower physical activity and shorter
   6-minute walk distance than the other Long COVID patients (PDE-045). The one free-living
   wearable study found (chest ECG + accelerometer, HR change at detected posture changes) reports cross-validated AUC 0.900
   (PDE-006), but in 66 POTS and 20 control *measurements*, with train/test splits at session level (the same person
   may be on both sides) and no external cohort. A 24-h ECG rule in Long COVID reports sensitivity 91.11% and
   specificity 96.36-99% (PDE-007), but with cut-offs chosen on the same data and healthy controls from another
   country.
2. **Imaging for endometriosis** has real diagnostic-accuracy evidence against an independent reference standard
   (surgery): transvaginal ultrasound for deep endometriosis, sensitivity 0.79 / specificity 0.94 over 934
   participants (Cochrane, PDE-059), endometrioma 0.93 / 0.96 (PDE-060). It is a clinic imaging test, and no imaging
   test was accurate enough to replace surgery for pelvic endometriosis overall.
3. **Gastric emptying devices for gastroparesis** (wireless capsule AUROC 0.83, PDE-056; breath test 89% sensitive at
   80% specificity, PDE-057). Like POTS, gastroparesis is defined by a device measurement, so this is agreement with
   the defining test.
4. **ME/CFS: 2-day CPET and tilt with cerebral blood flow** show group differences, not validated person-level
   diagnosis. The largest positive 2-day CPET study (84 vs 71, PDE-013) is contradicted by an independent null
   replication (58 vs 25, PDE-015); the only classification figure (95.1%, PDE-014) is in-sample with 10 controls.
   Tilt Doppler cerebral blood flow is abnormal in 90% of 429 ME/CFS patients against a cut-off built from 44 healthy
   controls (PDE-022). Heart rate on standing does **not** separate ME/CFS from matched healthy controls (33% of
   controls met the POTS HR criterion, PDE-021; no tilt difference in the NIH cohort, PDE-018).
5. **Consumer wearables in Long COVID** show consistent group differences in activity and HRV in the largest study
   (RECOVER Fitbit, 1,475 people: steps -1,624/day, HRV -4.4 ms, resting HR +1.5 bpm, PDE-039) but no person-level
   accuracy, a survey-index label and a small resting-HR effect. The pooled HRV literature is null (PDE-038). The
   one published wearable Long COVID classifier (AUROC 95.1% with symptoms, PDE-040) could not be reproduced from its
   public release: the HR window it needs starts at a diagnosis date the release lacks. A different, pre-specified
   whole-record resting-HR analysis by this project was null (AUROC 0.504).
6. **Fibromyalgia**: accelerometry separates groups (1,881 fewer steps/day, 413 vs 188, PDE-029); small-fibre pathology
   is found in about half of patients (49%, PDE-033 and PDE-034) but is not specific to fibromyalgia.
7. **Newer modalities** (VR pupillometry in Long COVID AUC 0.8711, PDE-049; bowel-sound belt for IBS 87% / 87% on a
   30-person test set, PDE-062) are small case-control proof-of-concept studies.

**None of the 17 diagnostic-accuracy rows reports validation in an external cohort** (13 primary studies, none
externally validated; 4 pooled reviews, where it does not apply), and the accuracy figures come from
case-control designs against healthy people, which inflate accuracy relative to the real clinical question
(this illness versus other causes of the same symptoms).

So the defensible claim is narrower than "a device diagnoses an invisible illness": **a wearable heart-rate +
posture measurement is a principled, published-evidence-supported way to uncover the objective orthostatic-tachycardia
subgroup (POTS) inside Long COVID and related illnesses, a subgroup that symptoms cannot identify.** This project has
not demonstrated it: no public dataset processed here has heart rate linked to posture (see
`results/PHASE3_MEASUREMENT_EVIDENCE.md`), and the only public Long COVID wearable label is null for resting heart
rate. The published evidence also says what a demonstration would need: symptomatic controls rather than healthy
volunteers, person-level splits, an age-appropriate threshold, and an external cohort.

## POTS is defined by a device: why that makes wearable HR + posture principled and partly circular

POTS is defined as chronic orthostatic intolerance with a heart-rate rise of at least 30 bpm on standing and no
orthostatic hypotension (PDE-001); the adult criteria do not suit adolescents, whose healthy 95th percentile on tilt was 42.9 bpm (PDE-005). The diagnosis is therefore a
device measurement (HR while posture changes) plus symptoms.

* **Principled.** A wearable that measures heart rate and detects posture measures the defining quantity itself, in
  daily life and repeatedly, instead of once in a clinic. Among the target conditions, no other has its defining
  criterion within reach of a consumer-grade sensor.
* **Partly circular.** Any study that scores a HR-based device against a POTS diagnosis compares the device with a
  reference standard built from the same quantity (incorporation bias). High accuracy then shows agreement with the
  criterion, not that the criterion captures an illness. It also inherits the criterion's weaknesses:
  * the 30 bpm threshold depends on test mode: 74% of diagnosed POTS patients met it on active stand vs 98% on tilt
    (PDE-004); tilt raises HR more and is less specific (40% specificity at 10 min vs 67% for standing, PDE-002/-003);
  * healthy young people cross it: 42% of healthy 8-19-year-olds on tilt (PDE-005), 33% of healthy adult controls on
    a NASA lean test (HR-only criterion, PDE-021), but 0 of 112 healthy adults when symptoms during the test are
    also required (PDE-066). This suggests the symptom requirement carries much of the specificity, but it compares
    two different studies, populations and protocols;
  * a wearable sees HR rises from any cause (exercise, stairs, anxiety, dehydration), so posture detection and
    context are part of the measurement, not an add-on (PDE-006 builds on detected posture-change events).
* **What escapes the circularity.** Two things a device adds that the clinic definition does not: (1) identifying who,
  among people with the same symptoms, has objective orthostatic tachycardia (PDE-045, PDE-046); (2) measuring it
  repeatedly in daily life. Both are measurement-deployment claims that fit this engine; neither is a claim that
  the device diagnoses Long COVID or ME/CFS.

## How the rows map to the measurement classes in `configs/measurements.yaml`

| class | rows | best published result | person-level accuracy published? |
|---|---|---|---|
| `autonomic_testing` | PDE-001-005, 009-011, 018, 020-024, 034, 035, 044-046, 048, 054, 055, 065, 066 | POTS stand test sensitivity 87% / specificity 67% vs healthy (n 15/15) | only for POTS, small, vs healthy |
| `posture_detection` | PDE-001, 002, 004, 006, 021, 028, 046, 055, 066 | wearable ECG + accelerometer, AUC 0.900 (PDE-006) | yes, one small feasibility study |
| `wearable_heart_rate` | PDE-001, 002, 006, 007, 020, 021, 039-042, 046, 052, 055, 066 | RECOVER RHR +1.5 bpm (group) | POTS only (PDE-006/007); Long COVID claim not reproduced (PDE-040) |
| `ecg_ambulatory` | PDE-006, 007, 008, 019, 045 | 24-h ECG POTS rule se 91.11% (PDE-007) | POTS only, in-sample |
| `hrv` | PDE-007, 008, 019, 020, 031, 032, 034, 038, 039, 041-043, 052, 054 | ME/CFS HF-HRV SMD -0.34 (PDE-020); Long COVID pooled null (PDE-038) | no (one null ML study, PDE-043) |
| `accelerometry` | PDE-006, 025-030, 039, 052 | FM -1,881 steps/day (PDE-029); LC -1,624 steps/day (PDE-039) | no |
| `sleep_objective` | PDE-027, 052 | ME/CFS lower sleep efficiency (PDE-027) | no |
| `ppg` | PDE-039, 052 | consumer HR accuracy best supported (PDE-052) | no |
| `continuous_spo2` | PDE-052 | measurement validity only | no |
| `cpet` | PDE-012-017, 030, 037, 064 | 2-day CPET group decline (PDE-012, -013) vs null replication (PDE-015) | in-sample only (PDE-014) |
| `small_fiber_testing` | PDE-011, 023, 033, 034, 047, 048, 065 | FM small-fibre pathology 49% (PDE-033) | no |
| `vascular_imaging` | PDE-022, 023 | ME/CFS tilt CBF abnormal in 90% (PDE-022) | sensitivity only, cut-off from the same controls |
| `digital_gait` | PDE-028, 036 | FM slower gait (PDE-036) | no |
| `digital_cognitive` | PDE-050, 053 | eye-movement correlations 0.210-0.371 (PDE-050) | no |
| no class (gap) | pupillometry PDE-049; voice PDE-051; gastric motility PDE-056-058; pelvic imaging PDE-059-061; bowel sound PDE-062, 063 | endometriosis TVUS 0.79 / 0.94 (PDE-059) | yes for endometriosis and gastroparesis |

The registry has no class for pupillometry, voice, gastric emptying / gastric electrophysiology, pelvic imaging or
bowel-sound analysis; `measurement_class_gap` records this per row. Adding classes is a curated-config change for the
owner of `configs/measurements.yaml`.

## Device x condition sections

Each section: what the literature reports, the quality of that evidence, and the measurement classes. "Healthy
controls" means volunteers without the illness; the clinically relevant comparison (the illness versus other causes of
the same symptoms) is almost never studied.

### 1. Active stand / tilt heart-rate increment x POTS

Classes: `autonomic_testing`, `posture_detection`, `wearable_heart_rate`.

* Reported: 30 bpm criterion at 10 min, active stand sensitivity 87% / specificity 67%, tilt 93% / 40%, in 15 POTS
  vs 15 healthy controls (PDE-002, PDE-003). In 60 diagnosed POTS patients the criterion was met on 74% of active
  stands vs 98% of tilts, 83% with a 27 bpm stand threshold (PDE-004). 42% of healthy children and adolescents exceed
  30 bpm on a 5-min tilt, 95th percentile 42.9 bpm (PDE-005). With symptoms required during a 10-min lean test, 0 of 112
  healthy adults met POTS or orthostatic-hypotension criteria (mean HR rise 9.89 ± 8.15 bpm, PDE-066).
* Quality: small case-control or normative studies; cases diagnosed with the same criterion (incorporation bias); no
  study against symptomatic non-POTS patients. Strong on physiology, weak as diagnostic-accuracy evidence.

### 2. Wearable / ambulatory HR (+ posture) x POTS, including POTS inside Long COVID

Classes: `ecg_ambulatory`, `wearable_heart_rate`, `accelerometry`, `posture_detection`, `hrv`.

* Reported: chest ECG + accelerometer worn >= 24 h in daily life, HR change at detected posture changes, deep
  learning: 10-fold CV AUC 0.900 ± 0.103, sensitivity 0.950, specificity 0.850 (PDE-006). 24-h ECG in Long COVID
  (HR spikes, HR after awakening, RMSSD): sensitivity 91.11%, specificity 99% vs healthy and 96.36% vs Long COVID
  without POTS; RMSSD alone had AUC > 0.50 only (PDE-007). Holter HRV index in adolescents: AUROC 0.874 (19 vs 44,
  PDE-008).
* Quality: all case-control, all internal validation only. PDE-006 counts measurements, splits by session, has 20
  control measurements, keeps the best-performing model per fold (metrics recorded every epoch), and reports 80%
  specificity on its hold-out test set (its ± values are called both SDs and 95% CIs in the paper); PDE-007 picks cut-offs on the same data and uses a foreign healthy database; PDE-008 reports the
  best of many indices. These are feasibility results that motivate, but do not establish, a wearable POTS screen.

### 3. Orthostatic tests x ME/CFS

Classes: `autonomic_testing`, `posture_detection`, `wearable_heart_rate`, `vascular_imaging`.

* NASA lean test, 150 ME/CFS vs 75 matched healthy controls: by 10 min 33% of controls, 42% of ME/CFS < 4 years and
  38% of ME/CFS > 10 years met POTS HR criteria; the difference was in symptoms and, for recent-onset cases, a
  narrowed pulse pressure (PDE-021). **Null for HR.**
* NIH post-infectious ME/CFS: tilt up to 40 min showed no group differences in orthostatic hypotension, orthostatic
  tachycardia or test-stopping symptoms (PDE-018). **Null.**
* Tilt with extracranial Doppler, 429 ME/CFS vs 44 healthy controls: end-tilt cerebral blood flow fell 26% vs 7%;
  90% of patients abnormal against a 13% cut-off (2 SD of the controls) (PDE-022). **Positive, single clinic,
  threshold from the same controls, no disease comparators.**
* Referral autonomic lab: reduced orthostatic CBF velocity in 88% of ME/CFS, POTS in 19% (PDE-023); Spanish cohort POTS
  in 31% of ME/CFS (PDE-024). Control rates for these flags are not in the abstracts.
* Quality: the one large positive result is cerebral blood flow, not heart rate. Heart rate on standing is not a
  discriminating ME/CFS measurement in these studies.

### 4. Orthostatic tests x Long COVID, PTLDS and EDS/HSD

Classes: `autonomic_testing`, `ecg_ambulatory`, `posture_detection`, `wearable_heart_rate`.

* Long COVID: POTS in 143 of 467 (31%) highly symptomatic patients with no symptom difference between POTS and
  non-POTS (PDE-045); positive lean test in 45 of 221 (20%), 83% of test-negative patients still symptomatic,
  algorithm vs expert agreement 92% (PDE-046); tilt HR 8 bpm higher than recovered comparators (95% CI 1.1 to 14.4),
  but only 3 of 16 (19%) clinically abnormal (PDE-044); POTS 22% in a referral autonomic lab (PDE-023), 13.8% in a
  Spanish cohort (PDE-024).
* PTLDS: 9 of 210 (4.29%) had orthostatic tachycardia on a 10-min active stand, not significantly different from
  healthy controls (PDE-055). **Null.**
* EDS/HSD: 48.6% of 35 hypermobility-type patients had postural tachycardia on tilt (PDE-054); no tilt controls.
* Quality: prevalence studies without controls or with small ones. The useful message is the Long COVID one: the
  test finds an objective subgroup that symptoms do not.

### 5. 2-day CPET (day-2 decline) and single CPET x ME/CFS

Class: `cpet`.

* 2-day CPET meta-analysis: 5 studies, 98 ME/CFS and 51 controls; patient-control workload difference at the
  ventilatory threshold -10.8 on test 1 vs -33.0 on test 2 (PDE-012). Largest study, 84 ME/CFS vs 71 sedentary
  controls: ME/CFS failed to reproduce CPET-1 on CPET-2, also in 55 fitness-matched pairs (PDE-013). Classification
  accuracy 95.1% in 51 CFS vs 10 controls, in-sample (PDE-014). **Independent replication, 58 vs 25: no day-2
  change in either group (peak VO2 22.3 -> 22.5 in ME/CFS); the authors reject 2-day CPET for defining PEM (PDE-015).**
* Single CPET: pooled peak VO2 5.2 mL/kg/min lower than healthy controls over 32 effects, 95% prediction interval
  -1.9 to 12.2 (PDE-016); in 214 vs 189 (CDC MCAM) many differences vanish after fitness matching, ventilatory
  inefficiency and perceived exertion remain (PDE-017); NIH cohort peak VO2 lower in 8 vs 9 (PDE-064).
* Quality: group differences only; the evidence for day-2 decline is contested; no study against other fatiguing
  illnesses or deconditioned patients; burdensome (two maximal tests, post-exertional risk).

### 6. CPET x Long COVID

Class: `cpet`. Meta-analysis of 9 studies (464 with symptoms vs 359 without): peak VO2 -4.9 mL/kg/min (95% CI -6.4 to
-3.4), rated low certainty with high risk of bias (PDE-037). Group difference; deconditioning common.

### 7. Actigraphy / step counts x ME/CFS

Classes: `accelerometry`, `sleep_objective`, `posture_detection`.

* 277 CFS vs 47 controls: shorter, less intense activity peaks; about one-fourth "pervasively passive" (PDE-025). FM
  and/or CFS 38 vs 27: lower peak activity (8,654 vs 12,913 units) but **similar average activity** (1,525 vs 1,602,
  P = 0.47) (PDE-026). 38 vs 38: lower sleep efficiency and more night-to-night variability; total sleep time not
  different (PDE-027). IMU "UpTime" separates 5 severe, 5 moderate and 5 healthy people (PDE-028). Systematic review:
  less daily activity (PDE-030).
* Quality: group differences, mixed on average activity, no discrimination statistics. Activity measures track
  severity and pacing better than they diagnose.

### 8. Actigraphy / step counts x fibromyalgia

Class: `accelerometry`. 413 women with fibromyalgia vs 188 controls: 1,881 ± 262 fewer steps/day and 39 ± 8 more
sedentary minutes/day (PDE-029); FM/CFS peak-only difference (PDE-026). Group-level separation only: the abstract
gives no person-level accuracy and no between-person spread of steps, so whether a ~1,900-step gap separates
individuals is not reported.

### 9. Consumer wearable HR / steps x Long COVID

Classes: `wearable_heart_rate`, `hrv`, `accelerometry`, `ppg`.

* RECOVER-Adult Fitbit, 1,475 people a median 21 months after infection, 498 with high symptom burden (LCRI):
  HRV -4.4 ms (-6.5 to -2.4), resting HR +1.5 bpm (+0.7 to +2.4), steps -1,624/day (-1,952 to -1,296) (PDE-039).
  Label is a symptom index among infected people; self-selected device sharers; no person-level accuracy.
* Stanford (Uwakwe 2025): combined HR + symptom model ROC-AUC 95.1%, about 5% above symptoms alone (PDE-040). **This
  project could not reproduce it:** the paper's HR features come from a post-diagnosis window, and the public
  release has no diagnosis date. This project's own, different, pre-specified whole-record resting-HR analysis was
  null: AUROC 0.504 (0.404-0.605), -0.024 over the demographic baseline. That null is not a failed replication of
  the paper's model.
* US military: persistent physiological change in 9.4% of 663 COVID-positive people vs 2,513 negative controls
  (PDE-041); the label is infection, not Long COVID.
* Quality: one large positive group-difference study; the only classifier claim could not be reproduced here, and
  this project's own resting-HR test on the same release was null.

### 10. HRV x Long COVID

Class: `hrv`. Meta-analysis, 11 studies, 593 vs 565: SDNN SMD 0.26 (-0.03 to 0.56), rMSSD 0.11 (-0.15 to 0.36), 9 of 11
at high risk of bias (PDE-038). **Null pooled result.** Wearable HRV lower in 121 Long COVID vs 21 controls (p = 0.027,
PDE-042); HRV-only machine learning did not separate post-COVID from healthy (F1 56%, PDE-043); RECOVER -4.4 ms
(PDE-039).

### 11. HRV / heart rate x ME/CFS

Classes: `hrv`, `ecg_ambulatory`, `wearable_heart_rate`. Meta-analysis of 64 articles: resting HR 4.14 ± 1.38 bpm
higher, tilt HR SMD 0.92 ± 0.24, resting HF-HRV SMD -0.34 ± 0.22 (PDE-020). NIH 24-h ECG: lower time-domain HRV in 14
vs 19 (PDE-019). Systematic review: in CFS, HRV reduced only during sleep (PDE-032). Group effects of this size overlap
heavily between individuals.

### 12. HRV x fibromyalgia

Class: `hrv`. Chronic-pain meta-analysis: consistent moderate-to-large decrease in high-frequency HRV, heavily
influenced by fibromyalgia studies (PDE-031; pooled number not in the open abstract); systematic review: lower HRV
in most of 10 FM case-control studies (PDE-032). Autonomic small-fibre impairment prevalence 45% (PDE-034).

### 13. Sudomotor and small-fibre testing x POTS

Classes: `autonomic_testing` (QSART, sweat tests, Sudoscan), `small_fiber_testing`. Half of 152 Mayo POTS patients had
sudomotor abnormalities (PDE-009); abnormal QSART in 17 of 30 (56%, PDE-010); neuropathic POTS (abnormal skin biopsy
IENFD + small-fibre function) in 9 of 24 (PDE-011). Prevalence among cases; these tests define a neuropathic subgroup,
not POTS.

### 14. Small-fibre testing x fibromyalgia

Classes: `small_fiber_testing`, `autonomic_testing`. Pooled prevalence of small-fibre pathology 49% (38-60%) in 222
patients: skin biopsy 45%, corneal confocal microscopy 59% (PDE-033); 20 studies, 903 patients: somatic 49%, autonomic
45% (PDE-034). Sudoscan hand conductance 71.4 vs 74.4 µS in 50 vs 50 (P = 0.003), no difference on feet (PDE-035).
Small-fibre pathology also occurs in other conditions, so prevalence is not specificity.

### 15. Small-fibre and sudomotor testing x Long COVID and ME/CFS

Classes: `small_fiber_testing`, `autonomic_testing`. Referral autonomic lab: small-fibre neuropathy 67% of Long COVID,
53% of ME/CFS (PDE-023); pathological Sudoscan palms 34% of ME/CFS and 19.5% of post-COVID (PDE-024); corneal nerve
fibre density lower in post-COVID patients with neurological symptoms (40 vs 30, p = 0.032, PDE-047); skin biopsy
abnormal in 10 of 16 referred Long COVID patients (PDE-048); 3 of 9 on a non-invasive protocol (PDE-065). Selected
referral samples, few controls.

### 16. Pupillometry x Long COVID

No class (`measurement_class_gap`). VR infrared pupillometry, 112 Long COVID, 44 post-COVID without Long COVID, 29
controls: constriction time alone AUC 0.8711 vs controls and 0.8140 vs post-COVID; selected features AUC 0.9000
(PDE-049). Feature and model selection in one small cohort; 29 controls; no external validation.

### 17. Eye tracking x ME/CFS and Long COVID

Class: `digital_cognitive`. ME/CFS: similar saccade latencies, impaired antisaccade accuracy and smooth pursuit
(PDE-053; n not in abstract). Post-COVID: eye-movement metrics correlate 0.210-0.371 with cognitive tests in 103
patients, no controls (PDE-050).

### 18. Voice x Long COVID (and ME/CFS)

No class. Smartphone voice in 64 post-COVID vs 70 controls: jitter, shimmer and HNR differed, speech timing did not
(PDE-051). For ME/CFS only a patient survey was found (not tabulated; DATA_AUDIT.md).

### 19. Gait x fibromyalgia (and ME/CFS)

Class: `digital_gait`. 36 studies, 3,369 FM vs 709 without: slower gait, shorter stride, lower cadence (PDE-036); the
pooled number is for the 6-minute walk test (stopwatch, not a device). ME/CFS UpTime IMU (PDE-028).

### 20. Gastric emptying devices x gastroparesis

No class. Wireless motility capsule, 61 gastroparesis vs 87 healthy: AUROC 0.83 (scintigraphy 4 h 0.82); 300-min
cut-off sensitivity 0.65, specificity 0.87 (PDE-056). Breath test vs simultaneous scintigraphy, 129 patients + 38
healthy: 89% sensitive for delayed emptying at 80% specificity (PDE-057). Body-surface gastric mapping found more
abnormalities than an emptying test (33.3% vs 22.7%, 42.7% combined in 75 patients) with no reference standard
(PDE-058). The condition is defined by a device measurement, so these are agreement figures.

### 21. Imaging x endometriosis

No class. Cochrane (surgical reference): deep endometriosis TVUS sensitivity 0.79 / specificity 0.94 (934
participants), endometrioma TVUS 0.93 / 0.96 and MRI 0.95 / 0.91; no imaging test met criteria to replace surgery
for pelvic endometriosis overall (PDE-059, PDE-060). Paired TVUS vs MRI, rectosigmoid deep disease: 0.85 / 0.96 vs
0.85 / 0.95 (PDE-061). The strongest diagnostic-accuracy evidence in the table, because the reference standard is
independent of the test; it applies to visible (deep, ovarian) disease.

### 22. Bowel-sound analysis x IBS

No class. Acoustic belt: 90% / 92% leave-one-out in 31 IBS + 37 healthy, 87% / 87% on the next 15 + 15 (PDE-062).
Meta-analysis of 4 studies: sensitivity 0.94, specificity 0.89, AUC 0.97, with the authors doubting validity
(PDE-063). Healthy controls only; the clinical problem is IBS versus other gut disease.

### 23. No device evidence found: MCAS, migraine, acute Lyme disease

Searches (DATA_AUDIT.md) found no device study with a diagnostic comparison for MCAS; for migraine, wearable studies
forecast attacks in patients already diagnosed. Acute Lyme disease was not searched for device studies (only
post-treatment Lyme, PDE-055). `get_published_device_evidence('MCAS')` returns UNKNOWN / NOT AVAILABLE.

## What this page does not show

* It does not show that any device diagnoses Long COVID, ME/CFS or fibromyalgia. For these the literature reports
  group differences, prevalence of abnormal tests among cases, or small in-sample classifiers.
* It does not replace a systematic review: searches were targeted and the inclusion rule was written during the
  search (DATA_AUDIT.md).
* Nothing here was reproduced by this project. PDE-040 was attempted and could not be reproduced from its public
  release.

## Query

```python
from measure_it.measurements.published_evidence import get_published_device_evidence
get_published_device_evidence("POTS", "wearable")      # rows whose classes intersect the wearable bundle
get_published_device_evidence("Long COVID", "pupillometry")   # no class: matched on device words
get_published_device_evidence("MCAS")                  # {'status': 'UNKNOWN / NOT AVAILABLE', ...}
```

Every returned row carries `claim_label` ("published claim — not reproduced by this project"), the design, reference
standard, risk-of-bias notes and the verbatim quote.

## Appendix: all rows

`?` = not reported in the text read (or a count summed by the curator; see `count_derivation`). Numbers are as published, not reproduced.

| id | condition(s) | device / test | claim | direction | n cases / controls | headline (as published) | ext. valid. | source |
|---|---|---|---|---|---|---|---|---|
| PDE-001 | pots | Orthostatic heart-rate increment (active stand or head-up tilt; any device that measures HR and posture) | defining_criterion | not_applicable | ? / ? | >= 30 bpm on standing, chronic (> 3 months) orthostatic intolerance, no orthostatic hypotension | not_applicable | Choi H 2026, PMID 41758178 |
| PDE-002 | pots | 10-min active stand, 30 bpm criterion | diagnostic_accuracy | positive | 15 / 15 | STAND10 sensitivity 87%, specificity 67% | no | Plash WB 2013, PMID 22931296 |
| PDE-003 | pots | 10-min head-up tilt, 30 bpm criterion | diagnostic_accuracy | mixed | 15 / 15 | TILT10 sensitivity 93%, specificity 40% | no | Plash WB 2013, PMID 22931296 |
| PDE-004 | pots | 10-min active stand vs head-up tilt (beat-to-beat haemodynamics), same day | test_comparison | mixed | 60 / 0 | active stand 74% vs tilt 98%; 83% vs 98% with a 27 bpm stand threshold | no | Uppal J 2025, PMID 40273723 |
| PDE-005 | pots | 5-min 70-degree head-up tilt, orthostatic HR increment | normative_false_positive_rate | null | 654 / 106 | 42% of normal controls had an HR increment >= 30 bpm; 95th percentile 42.9 bpm | not_applicable | Singer W 2012, PMID 21996154 |
| PDE-006 | pots | Chest-worn wearable ECG + accelerometer (Faros 180), free-living >= 24 h, HR change at detected posture changes, deep learning | diagnostic_accuracy | positive | 66 / 20 | AUC 0.900 ± 0.103; sensitivity 0.950 ± 0.062; specificity 0.850 ± 0.142 | no | Choi H 2026, PMID 41758178 |
| PDE-007 | pots, long_covid | 24-h ambulatory ECG: HR spikes > 30 bpm, HR rise after awakening, RMSSD (combined rule) | diagnostic_accuracy | positive | 45 / 155 | se 91.11%, sp 99% vs healthy; sp 96.36% vs LC without POTS | no | Hupin D 2025, PMID 40734755 |
| PDE-008 | pots | Holter ECG HRV (Kubios PNS index) | diagnostic_accuracy | positive | 19 / 44 | AUROC 0.874; sensitivity 77.8%; specificity 88.6% at Youden threshold | no | Kakavand B 2026, PMID 42380672 |
| PDE-009 | pots | QSART + thermoregulatory sweat test (sudomotor) | prevalence_in_cases | positive | 152 / 0 | half | not_applicable | Thieben MJ 2007, PMID 17352367 |
| PDE-010 | pots | QSART | prevalence_in_cases | positive | 30 / 0 | 17 of 30 (56%) | not_applicable | Peltier AC 2010, PMID 20035362 |
| PDE-011 | pots | Skin biopsy IENFD + QST + sudomotor testing | prevalence_in_cases | positive | 24 / 10 | 9 of 24 | not_applicable | Gibbons CH 2013, PMID 24386408 |
| PDE-012 | me_cfs | 2-day CPET (day-2 decline) | group_difference | positive | 98 / 51 | -10.8 at test 1 vs -33.0 at test 2 (p < 0.05) | not_applicable | Lim EJ 2020, PMID 33327624 |
| PDE-013 | me_cfs | 2-day CPET | group_difference | positive | 84 / 71 | ME/CFS failed to reproduce CPET-1 on CPET-2 (work, VO2, HR and others); controls reproduced; results similar in 55 fitness-matched pairs | not_applicable | Keller B 2024, PMID 38965566 |
| PDE-014 | me_cfs | 2-day CPET, multivariate classification | diagnostic_accuracy | positive | 51 / 10 | 95.1% overall accuracy; no group difference on test 1 | no | Snell CR 2013, PMID 23813081 |
| PDE-015 | me_cfs | 2-day CPET | group_difference | null | 58 / 25 | ME/CFS 22.3 -> 22.5; controls 23.4 -> 22.8; NS | not_applicable | Mancini DM 2026, PMID 42212259 |
| PDE-016 | me_cfs | Single maximal CPET (peak VO2) | group_difference | positive | ? / ? | 5.2 lower | not_applicable | Franklin JD 2019, PMID 30557887 |
| PDE-017 | me_cfs | Single ramped cycle CPET (fitness-matched comparison) | group_difference | mixed | 214 / 189 | breathing frequency, ventilatory equivalents and RPE still differed; HR and CI differences did not survive fitness matching | not_applicable | Cook DB 2022, PMID 35290404 |
| PDE-018 | me_cfs | Head-up tilt up to 40 min | group_difference | null | ? / ? | no group differences | not_applicable | Walitt B 2024, PMID 38383456 |
| PDE-019 | me_cfs | 24-h ambulatory ECG HRV | group_difference | positive | 14 / 19 | lower in PI-ME/CFS (rMSSD p = 0.019, unadjusted) | not_applicable | Walitt B 2024, PMID 38383456 |
| PDE-020 | me_cfs | HR and HRV (resting, tilt, orthostatic, exercise; ECG-based) | group_difference | positive | ? / ? | resting HR MD 4.14 ± 1.38 bpm higher; HR on tilt SMD 0.92 ± 0.24; resting HF-HRV SMD -0.34 ± 0.22 | not_applicable | Nelson MJ 2019, PMID 31651868 |
| PDE-021 | me_cfs | 10-min NASA Lean Test (BP + HR each minute) | group_difference | null | 150 / 75 | 33% of healthy controls, 42% of ME/CFS < 4 y, 38% of ME/CFS > 10 y; narrowed pulse pressure only in the < 4 y group | not_applicable | Lee J 2020, PMID 32799889 |
| PDE-022 | me_cfs | Head-up tilt with extracranial Doppler cerebral blood flow | diagnostic_accuracy | positive | 429 / 44 | 90% of ME/CFS abnormal; mean end-tilt reduction 26% vs 7% in controls | no | van Campen CLMC 2020, PMID 32140630 |
| PDE-023 | me_cfs, long_covid | Autonomic battery (tilt with transcranial Doppler CBFv, Valsalva, deep breathing, sudomotor) + skin biopsy | prevalence_in_cases | positive | 313 / 73 | reduced orthostatic CBFv 92%/88%; SFN 67%/53%; POTS 22%/19% | not_applicable | Novak P 2026, PMID 41576003 |
| PDE-024 | me_cfs, long_covid | Sudoscan electrochemical skin conductance; tilt; heat/cold evoked potentials | prevalence_in_cases | positive | 137 / 50 | POTS 31% / 13.8%; Sudoscan palms 34% / 19.5% | not_applicable | Azcue N 2023, PMID 37968647 |
| PDE-025 | me_cfs | 12-day actometer (ankle movement sensor) | group_difference | positive | 277 / 47 | less intense, shorter peaks; approximately one-fourth pervasively passive | not_applicable | van der Werf SP 2000, PMID 11164063 |
| PDE-026 | me_cfs, fibromyalgia | 5-day ambulatory activity monitor | group_difference | mixed | 38 / 27 | peak activity 8,654 vs 12,913 units (P = 0.003); average activity 1,525 vs 1,602 (P = 0.47) | not_applicable | Kop WJ 2005, PMID 15641057 |
| PDE-027 | me_cfs | 7-day wrist accelerometer sleep | group_difference | positive | 38 / 38 | poorer sleep efficiency and greater night-to-night variability; total sleep time did not differ | not_applicable | Saurel M 2026, PMID 42129014 |
| PDE-028 | me_cfs | Leg-worn IMU 'UpTime' (time with feet on the floor) | severity_stratification | positive | 10 / 5 | severe < 20%; moderate 20-30%; healthy > 30% | no | Palombo T 2020, PMID 33168001 |
| PDE-029 | fibromyalgia | Triaxial accelerometry (steps, sedentary time, MVPA) | group_difference | positive | 413 / 188 | -1,881 ± 262 steps/day; +39 ± 8 sedentary min/day | not_applicable | Segura-Jiménez V 2015, PMID 26108350 |
| PDE-030 | me_cfs | Objective physical activity, exercise capacity, strength (systematic review) | group_difference | positive | ? / ? | less daily physical activity; conflicting exercise-capacity data | not_applicable | Nijs J 2011, PMID 21166613 |
| PDE-031 | fibromyalgia | HRV (ECG) | group_difference | positive | ? / ? | consistent moderate-to-large decrease (numeric pooled value not in the abstract) | not_applicable | Tracy LM 2016, PMID 26431423 |
| PDE-032 | fibromyalgia, me_cfs | HRV (ECG) | group_difference | mixed | ? / ? | FM: lower HRV in most studies; CFS: HRV reduced only during sleep | not_applicable | Meeus M 2013, PMID 23838093 |
| PDE-033 | fibromyalgia | Skin biopsy IENFD and corneal confocal microscopy (small fibre pathology) | prevalence_in_cases | positive | 222 / ? | 49% overall; skin biopsy 45%; corneal confocal microscopy 59% | not_applicable | Grayston R 2019, PMID 30314675 |
| PDE-034 | fibromyalgia | Somatic (skin biopsy, CCM, microneurography) and autonomic (HRV, SSR, skin conductance, tilt) small-fibre tests | prevalence_in_cases | positive | 903 / ? | somatic 49%; autonomic 45% | not_applicable | Galosi E 2022, PMID 35626288 |
| PDE-035 | fibromyalgia | Sudoscan electrochemical skin conductance | group_difference | mixed | 50 / 50 | 71.4 ± 11.2 vs 74.4 ± 10.3 µS (P = 0.003); no difference on feet | not_applicable | Pickering G 2020, PMID 31705738 |
| PDE-036 | fibromyalgia | Gait analysis and walk tests | group_difference | positive | 3369 / 709 | slower gait, shorter stride, lower cadence; shorter 6MWT distance | not_applicable | Carrasco-Vega E 2022, PMID 35341042 |
| PDE-037 | long_covid | CPET peak VO2 | group_difference | positive | 464 / 359 | -4.9 | not_applicable | Durstenfeld MS 2022, PMID 36223120 |
| PDE-038 | long_covid | HRV (SDNN, rMSSD, LF/HF) | group_difference | null | 593 / 565 | SDNN 0.26 (p = 0.09); rMSSD 0.11 (p = 0.41); LF/HF -0.271 (p = 0.12) | not_applicable | Schoene D 2026, PMID 41814525 |
| PDE-039 | long_covid | Consumer Fitbit (resting HR, HRV, steps, activity, SpO2), 6-month window | group_difference | positive | 498 / 977 | HRV -4.4 ms; resting HR +1.5 bpm; steps -1,624/day | not_applicable | Vogel JM 2026, PMID 42330737 |
| PDE-040 | long_covid | Consumer smartwatch HR features + symptom features (ML) | diagnostic_accuracy | positive | 31 / 95 | 95.1% (about 5% above symptoms-only) | no | Uwakwe CK 2025, PMID 41264615 |
| PDE-041 | long_covid | Wearable (nightly HR, HRV) | prevalence_in_cases | positive | 663 / 2513 | 9.4% | no | Borhani S 2025, PMID 40790065 |
| PDE-042 | long_covid | Continuous multi-day wearable HRV | group_difference | positive | 121 / 21 | lower in Long COVID (p = 0.027) | not_applicable | Ruijgt TM 2026, PMID 42501245 |
| PDE-043 | long_covid | Wearable HRV indices (SDNN, RMSSD, LF%, HF%) + ML | diagnostic_accuracy | null | 20 / 20 | differences not significant; F1 56% for post-COVID detection from HRV alone | no | Sanches CA 2025, PMID 40553043 |
| PDE-044 | long_covid | 30-min 70-degree head-up tilt | group_difference | mixed | 16 / 10 | HR 8 bpm higher throughout tilt; 3 of 16 (19%) clinically abnormal | not_applicable | Durstenfeld MS 2025, PMID 41134786 |
| PDE-045 | long_covid, pots | 48-h ECG + head-up tilt + active stand in clinically suspected POTS | prevalence_in_cases | positive | 467 / 0 | 143 of 467 (31%) | not_applicable | Björnson M 2025, PMID 41025260 |
| PDE-046 | long_covid, pots | 10-min NASA Lean Test with automated algorithm | prevalence_in_cases | positive | 221 / 0 | 45 of 221 (20%) positive by algorithm; 83% of NLT-negative patients still reported orthostatic symptoms; algorithm vs experts 92% | not_applicable | Siebler M 2026, PMID 42387884 |
| PDE-047 | long_covid | Corneal confocal microscopy (corneal nerve fibre density) | group_difference | positive | 40 / 30 | lower (p = 0.032) | not_applicable | Bitirgen G 2022, PMID 34312122 |
| PDE-048 | long_covid | Skin biopsy IENFD, electrodiagnostics, autonomic function tests | prevalence_in_cases | positive | 17 / 0 | 10 of 16 (63%) | not_applicable | Oaklander AL 2022, PMID 35232750 |
| PDE-049 | long_covid | Virtual-reality infrared pupillometry (pupillary light response) | diagnostic_accuracy | positive | 112 / 73 | AUC 0.8711; selected features AUC 0.9000 | no | Tang CH 2025, PMID 39631631 |
| PDE-050 | long_covid | Eye tracker (fixation, saccades, smooth pursuit, pupil) vs neuropsychological tests | within_patient_correlation | mixed | 103 / 0 | 0.210 to 0.371 | not_applicable | Goset J 2026, PMID 41851162 |
| PDE-051 | long_covid | Smartphone voice recording, acoustic analysis (jitter, shimmer, HNR) | group_difference | mixed | 64 / 70 | jitter, shimmer, HNR differed; speech parameters did not | not_applicable | Berti LC 2025, PMID 39791670 |
| PDE-052 | long_covid, post_infectious_syndrome | Consumer wearables (17 biometrics): accuracy and validity | measurement_validity | mixed | ? / ? | heart rate measurement and atrial fibrillation detection best supported; limited evidence for clinical utility | not_applicable | Kaplan DM 2026, PMID 42258488 |
| PDE-053 | me_cfs | Eye tracking (saccades, antisaccades, smooth pursuit) | group_difference | mixed | ? / ? | saccade latencies similar; antisaccade accuracy and smooth pursuit impaired | not_applicable | Badham SP 2013, PMID 23918092 |
| PDE-054 | eds_hsd, pots | Cardiovascular reflex battery incl. head-up tilt; spectral HR/BP variability | prevalence_in_cases | positive | 35 / ? | 48.6% postural orthostatic tachycardia; 31.4% orthostatic intolerance; 20% normal | not_applicable | Celletti C 2017, PMID 28286774 |
| PDE-055 | ptlds | 10-min active stand test | group_difference | null | 210 / ? | 9 of 210 (4.29%); not significantly different from healthy controls | not_applicable | Adler BL 2025, PMID 41427440 |
| PDE-056 | gastroparesis | Wireless motility capsule (SmartPill) gastric emptying time | diagnostic_accuracy | positive | 61 / 87 | GET AUROC 0.83 (GES-4h 0.82); 300-min cut-off sensitivity 0.65, specificity 0.87 | no | Kuo B 2008, PMID 17973643 |
| PDE-057 | gastroparesis | [13C]-Spirulina gastric emptying breath test | diagnostic_accuracy | positive | 129 / 38 | 89% (150- and 180-min samples combined) | no | Szarka LA 2008, PMID 18406670 |
| PDE-058 | gastroparesis | Body-surface gastric mapping (Gastric Alimetry) vs gastric emptying test | test_comparison | mixed | 75 / 0 | GET 22.7%; Gastric Alimetry spectral 33.3%; combined 42.7% | not_applicable | Wang WJ 2024, PMID 37782524 |
| PDE-059 | endometriosis | Transvaginal ultrasound | diagnostic_accuracy | positive | ? / ? | sensitivity 0.79, specificity 0.94 | not_applicable | Nisenblat V 2016, PMID 26919512 |
| PDE-060 | endometriosis | Transvaginal ultrasound; MRI | diagnostic_accuracy | positive | ? / ? | TVUS 0.93 / 0.96; MRI 0.95 / 0.91 | not_applicable | Nisenblat V 2016, PMID 26919512 |
| PDE-061 | endometriosis | Transvaginal ultrasound vs MRI in the same patients | diagnostic_accuracy | positive | ? / ? | TVS 0.85 / 0.96; MRI 0.85 / 0.95 | not_applicable | Guerriero S 2018, PMID 29154402 |
| PDE-062 | ibs | Acoustic belt bowel-sound analysis | diagnostic_accuracy | positive | 46 / 52 | 87% / 87% (15 IBS + 15 healthy); LOOCV 90% / 92% | no | Du X 2019, PMID 30908308 |
| PDE-063 | ibs | Computerized bowel sound analysis (non-invasive devices) | diagnostic_accuracy | positive | ? / ? | 0.94 / 0.89 / AUC 0.97 | not_applicable | Yan XX 2024, PMID 38594814 |
| PDE-064 | me_cfs | Single maximal CPET | group_difference | positive | 8 / 9 | lower in PI-ME/CFS (p = 0.004), about 3.3 METs; peak HR p = 0.07 | not_applicable | Walitt B 2024, PMID 38383456 |
| PDE-065 | long_covid | Non-invasive small-fibre protocol (sympathetic skin response, cutaneous silent period, thermal thresholds, Sudoscan ESC) | prevalence_in_cases | positive | 9 / 0 | 3 of 9 (33%) | not_applicable | Khoo A 2026, PMID 41635577 |
| PDE-066 | pots | 10-min lean test (NASA-type), BP and HR every minute | normative_false_positive_rate | positive | 0 / 112 | none of 112; HR rise 9.89 ± 8.15 bpm | not_applicable | Iftekhar N 2025, PMID 41111962 |
