#!/usr/bin/env Rscript
# ===========================================================================
# Synthetic health-administrative data generator for the Quan Summit workshop
# "Defining Study Populations and Variables".
#
# Produces two linkable files, for teaching purposes only. No real patient
# data is involved and no row corresponds to any real person.
#
#   1. A simplified CIHI Discharge Abstract Database (DAD) extract,
#      diagnoses coded in ICD-10-CA.
#   2. A simplified physician-claims extract, diagnoses coded in ICD-9
#      (up to three diagnosis fields per claim, three digits each, as in
#      Alberta practitioner claims), linkable to the DAD on PATIENT_ID.
#
# DAD columnr
# -----------
#   PATIENT_ID       synthetic person-level key (links the two files, and
#                    lets you find repeat hospitalizations)
#   ENCOUNTER_ID     synthetic encounter-level key
#   AGE              age in years at admission
#   SEX              M / F / O / U  (DAD convention)
#   ADMIT_DTTM       admission date + time
#   DISCH_DTTM       discharge date + time
#   ADMIT_CATEGORY   U urgent, L elective, N newborn
#   INST_ID          facility identifier (transfers move between facilities)
#   DISCH_DISP       discharge disposition (see below)
#   DXCODE1..25      ICD-10-CA diagnosis codes, left-filled. DXCODE1 is
#                    always populated and is always the most responsible
#                    diagnosis. There are no DXTYPE columns: this extract
#                    does not carry diagnosis type.
#
# Claims columns
# --------------
#   CLAIM_ID, PATIENT_ID, SERVICE_DATE, SERVICE_LOCATION, PROVIDER_TYPE,
#   DXCODE1..3       ICD-9 diagnosis codes, left-filled, DXCODE1 populated
#
# Usage
# -----
#   Rscript generate_icd10ca_data.R --patients 20000
#   Rscript generate_icd10ca_data.R --patients 2000 --truth --report
#   Rscript generate_icd10ca_data.R --patients 5000 --hf-rate 0.25 --seed 7
#
# Base R only, no packages. Tested on R 4.3.
#
# Comorbidity coverage
# --------------------
# The condition registry below is built so that abstracts and claims carry
# codes the published algorithms will actually find:
#
#   * Quan H et al. Coding algorithms for defining comorbidities in
#     ICD-9-CM and ICD-10 administrative data. Med Care 2005;43:1130-9
#     -- the enhanced Charlson (17 conditions) and Elixhauser (31).
#   * Tonelli M et al. Methods for identifying 30 chronic conditions:
#     application to administrative data. BMC Med Inform Decis Mak
#     2015;15:31, with the 2019 correction.
#
# Every condition lists the algorithm categories it feeds (`algos`), and
# the emitted codes are drawn from those published lists, so a student who
# implements Quan's or Tonelli's code lists will pick them up. A few
# Elixhauser categories are acute rather than chronic (fluid and
# electrolyte disorders, blood loss and deficiency anaemia, coagulopathy,
# weight loss) and are fed from the acute code pool instead; they are
# marked there.
#
# Two honest limitations, both worth saying out loud in the workshop:
#   * Claims carry three-digit ICD-9, which is what Alberta practitioner
#     claims actually hold. Algorithms specifying 4- or 5-digit ICD-9 codes
#     therefore cannot be expressed exactly (Tonelli's peptic ulcer 531.7,
#     chronic pain 338.0, atrial fibrillation 427.31, and so on). Three
#     digits is what the data give you; the widening is real.
#   * Tonelli treats some conditions as remitting (cancers after 5 years,
#     depression and peptic ulcer after 2). Conditions here are permanent
#     once assigned, so a remission rule applied to this data will not
#     change anything.
#
# Notes on ICD-10-CA vs ICD-10-CM -- these matter for teaching:
#   * ICD-10-CA has no I50.2x / I50.3x / I50.4x systolic-vs-diastolic
#     subdivisions. Heart failure is I50.0, I50.1, I50.9 only.
#   * There are no 7th-character extensions and no "x" placeholders.
#   * Codes are 3 to 6 characters; Canadian detail is added at the 5th/6th
#     character rather than by restructuring the category.
#   * Laterality is generally NOT coded, unlike ICD-10-CM.
# Codes are stored WITH the decimal for readability and stripped on output
# unless --decimal is passed, because the DAD carries codes without it.
#
# What is deliberately built into the data (the point of the workshop)
# --------------------------------------------------------------------
#   * Comorbidity is a latent patient-level truth. Each condition is then
#     recorded in the DAD with a condition-specific sensitivity, and in the
#     claims with a different one. Hypertension is given a low DAD
#     sensitivity on purpose, so a DAD-only definition misses most cases
#     and the claims-based definition (Quan 2009: 2 claims in 2 years)
#     recovers them. --truth writes the latent flags out so you can score
#     whatever algorithm the workshop builds.
#   * Patients have repeat hospitalizations, and the gap to the next
#     admission depends on age and comorbidity, so a survival model fit to
#     time-to-readmission has something real to find.
#   * Acute-care transfers appear as a discharge with disposition 01
#     followed by an admission at a different facility a few hours later.
#     Counting those as readmissions is the classic error the episode-of-
#     care step is meant to catch.
#   * In-hospital deaths end a patient's record and are a competing risk.
#   * Admissions stop at --end, so follow-up is administratively censored.
#
# Discharge disposition (DISCH_DISP): 01 transfer to acute care, 02
# transfer to continuing care, 03 transfer to other, 04 home with support
# services, 05 home, 06 signed out (against medical advice), 07 died.
# NOTE: this is the pre-2018-19 CIHI code set. A real 2022-2024 extract
# uses the newer scheme (10/20/30/40/90 transfers, 61/62/65 left, 72/73/74
# died) -- check the abstracting manual for the year of data you hold.
#
# Claims SERVICE_LOCATION: OFFICE, ED, INPT (billed during an inpatient
# stay), LTC, HOME.
# ===========================================================================

Sys.setenv(TZ = "UTC")
MAX_DX <- 25
MAX_CLAIM_DX <- 3
MAX_ENC_PER_PATIENT <- 14
SECS_PER_DAY <- 86400

inv_logit <- function(x) 1 / (1 + exp(-x))

# ---------------------------------------------------------------------------
# Heart failure as the most responsible diagnosis
# ---------------------------------------------------------------------------
# Drawn from the Quan 2005 Charlson congestive heart failure list (I09.9,
# I11.0, I13.0, I13.2, I25.5, I42.0, I42.5-I42.9, I43, I50, P29.0). I50.x
# dominates, but I42.0 and I25.5 appear too, so "I50 only" and "the Quan
# list" give different cohorts -- which is the point of the code-list
# sensitivity exercise. I11.0 and I13.0 are used only for patients who
# really have hypertension (and, for I13.0, kidney disease): a combination
# code should not appear for someone without the condition it combines.
HF_CODES        <- c("I50.0", "I50.9", "I50.1", "I42.0", "I25.5")
HF_CODE_W       <- c(0.44,    0.36,    0.13,    0.05,    0.02)
HF_HTN_CODE     <- "I11.0"   # hypertensive heart disease with heart failure
HF_HTN_CKD_CODE <- "I13.0"   # hypertensive heart and renal disease with HF

# ---------------------------------------------------------------------------
# The chronic condition registry
# ---------------------------------------------------------------------------
# Prevalence is inv_logit(logit65 + slope * (age - 65) + hf_boost * is_hf),
# i.e. logit65 is the log-odds of having the condition at age 65.
#
#   algos       which published algorithm(s) this condition feeds
#   known_p     chance the condition has actually been diagnosed. An
#               undiagnosed condition is invisible to every administrative
#               source, which puts a ceiling on every algorithm's sensitivity
#   dad_sens    chance a diagnosed condition is coded on a given abstract
#   dad_fp      chance it is coded for someone who does not have it
#   claims_per_year      expected claims carrying the condition
#   claims_fp_per_year   same, for someone who does not have it
#   readmit_beta         log-hazard contribution to the next admission
#   sex         restrict the condition to one DAD sex code
#
# ICD-9 codes are three digits, matching what the claims file carries.
cond <- function(label, algos = character(0), logit65, slope, hf_boost = 0,
                 min_age = 18, icd10, icd9, icd10_w = NULL, icd9_w = NULL,
                 known_p = 0.80, dad_sens = 0.55, dad_fp = NULL,
                 claims_per_year = 1.4, claims_fp_per_year = NULL,
                 provider = c(GP = 0.80, SPEC = 0.20), sex = NA,
                 readmit_beta = 0) {
  if (is.null(icd10_w)) icd10_w <- rep(1, length(icd10))
  if (is.null(icd9_w))  icd9_w  <- rep(1, length(icd9))
  # Miscoding scales with how common the condition is. A flat false-positive
  # rate would give a rare condition more false cases than true ones, which
  # is the opposite of how administrative data behave: specificity for
  # chronic conditions is consistently above 99%.
  p65 <- 1 / (1 + exp(-logit65))
  if (is.null(dad_fp)) dad_fp <- 0.020 * p65
  if (is.null(claims_fp_per_year)) claims_fp_per_year <- 0.12 * p65
  list(label = label, algos = algos, logit65 = logit65, slope = slope,
       hf_boost = hf_boost, min_age = min_age, icd10 = icd10, icd10_w = icd10_w,
       icd9 = icd9, icd9_w = icd9_w, known_p = known_p, dad_sens = dad_sens,
       dad_fp = dad_fp, claims_per_year = claims_per_year,
       claims_fp_per_year = claims_fp_per_year, provider = provider,
       sex = sex, readmit_beta = readmit_beta)
}

CHRONIC <- list(

  # --- circulatory -------------------------------------------------------
  HTN = cond(
    "Hypertension", c("elixhauser", "tonelli"),
    logit65 = 0.30, slope = 0.050, hf_boost = 1.00, min_age = 25,
    # Low DAD sensitivity on purpose: hypertension is rarely the reason for
    # admission and rarely drives resource use, so coders often leave it off.
    known_p = 0.86, dad_sens = 0.34, dad_fp = 0.010,
    icd10 = c("I10", "I11.9", "I12.9"), icd10_w = c(8, 1.5, 0.5),
    icd9  = c("401", "402", "403"),     icd9_w  = c(8, 1.5, 0.5),
    claims_per_year = 2.4, claims_fp_per_year = 0.12, readmit_beta = 0.45),
  IHD = cond(
    "Ischaemic heart disease", character(0),
    logit65 = -2.20, slope = 0.045, hf_boost = 1.20, min_age = 35,
    known_p = 0.82, dad_sens = 0.74,
    icd10 = c("I25.1", "I20.0"), icd10_w = c(8, 1),
    icd9  = c("414", "413"),     icd9_w  = c(8, 2),
    claims_per_year = 1.4, provider = c(GP = 0.55, SPEC = 0.45)),
  MI_PRIOR = cond(
    # Charlson counts old MI (I25.2) as well as the acute codes; Tonelli
    # requires a most-responsible hospitalization for I21/I22, which comes
    # from the acute pool instead.
    "Previous myocardial infarction", c("charlson", "tonelli"),
    logit65 = -2.90, slope = 0.045, hf_boost = 1.10, min_age = 35,
    known_p = 0.95, dad_sens = 0.62,
    icd10 = c("I25.2"), icd9 = c("412"),
    claims_per_year = 0.9, provider = c(GP = 0.6, SPEC = 0.4)),
  AFIB = cond(
    "Atrial fibrillation", c("elixhauser", "tonelli"),
    logit65 = -3.00, slope = 0.055, hf_boost = 1.30, min_age = 40,
    known_p = 0.85, dad_sens = 0.80,
    # Tonelli's atrial fibrillation list is I48.0 / ICD-9 427.3, so I48.0
    # carries most of the weight here.
    icd10 = c("I48.0", "I48.1", "I48.2", "I48.9"), icd10_w = c(5, 1.5, 2, 2),
    icd9  = c("427"),
    claims_per_year = 1.6, provider = c(GP = 0.5, SPEC = 0.5),
    readmit_beta = 0.20),
  VALVULAR = cond(
    "Valvular disease", "elixhauser",
    logit65 = -3.50, slope = 0.055, hf_boost = 1.20, min_age = 40,
    known_p = 0.85, dad_sens = 0.66,
    icd10 = c("I35.0", "I34.0", "I05.0"), icd10_w = c(5, 4, 1),
    icd9  = c("424"),
    claims_per_year = 1.1, provider = c(GP = 0.45, SPEC = 0.55)),
  PULM_CIRC = cond(
    "Pulmonary circulation disorders", "elixhauser",
    logit65 = -4.40, slope = 0.040, hf_boost = 1.40, min_age = 30,
    known_p = 0.80, dad_sens = 0.60,
    # I27.8/I27.9 are in Quan's chronic pulmonary disease list too, so the
    # codes used here are the ones that are specific to this category.
    icd10 = c("I27.0", "I28.8"), icd10_w = c(4, 1),
    icd9  = c("416"),
    claims_per_year = 1.0, provider = c(GP = 0.4, SPEC = 0.6)),
  PVD = cond(
    "Peripheral vascular disease", c("charlson", "elixhauser", "tonelli"),
    logit65 = -3.50, slope = 0.050, hf_boost = 0.60, min_age = 40,
    known_p = 0.78, dad_sens = 0.68,
    # Tonelli's peripheral vascular disease is exactly I70.2 / ICD-9 440.2.
    icd10 = c("I70.2", "I73.9", "I71.4"), icd10_w = c(7, 2, 1),
    icd9  = c("440"),
    claims_per_year = 1.1, provider = c(GP = 0.6, SPEC = 0.4)),
  STROKE_TIA = cond(
    # Charlson's cerebrovascular category spans I60-I69, so sequelae codes
    # count. Tonelli requires an acute stroke code (I60/I61/I63/I64/G45),
    # which comes from the acute pool, or a physician claim.
    "Previous stroke or TIA", c("charlson", "tonelli"),
    logit65 = -2.90, slope = 0.055, hf_boost = 0.50, min_age = 40,
    known_p = 0.92, dad_sens = 0.64,
    icd10 = c("I69.3", "I69.4"), icd10_w = c(4, 1),
    icd9  = c("434", "435"), icd9_w = c(3, 2),
    claims_per_year = 1.0, provider = c(GP = 0.6, SPEC = 0.4)),

  # --- respiratory -------------------------------------------------------
  COPD = cond(
    "COPD", c("charlson", "elixhauser", "tonelli"),
    logit65 = -2.60, slope = 0.040, hf_boost = 0.90, min_age = 40,
    known_p = 0.72, dad_sens = 0.72, dad_fp = 0.012,
    icd10 = c("J44.0", "J44.1", "J44.9"), icd10_w = c(1.2, 1.8, 3.0),
    icd9  = c("491", "492", "496"),       icd9_w  = c(3, 1, 4),
    claims_per_year = 1.8, readmit_beta = 0.45),
  ASTHMA = cond(
    # J45 is in Charlson's J40-J47 chronic pulmonary range but Tonelli
    # deliberately removes it from chronic pulmonary disease, because
    # asthma has its own algorithm. Same patient, two answers.
    "Asthma", c("charlson", "elixhauser", "tonelli"),
    logit65 = -2.40, slope = -0.010, hf_boost = 0.20, min_age = 5,
    known_p = 0.85, dad_sens = 0.55,
    icd10 = c("J45.9", "J45.0"), icd10_w = c(4, 1),
    icd9  = c("493"),
    claims_per_year = 1.3),

  # --- endocrine / metabolic ---------------------------------------------
  DIABETES = cond(
    # Charlson and Elixhauser both split uncomplicated from complicated;
    # the split here is in the codes, not in the latent truth, so students
    # see the distinction come out of the data.
    "Diabetes mellitus", c("charlson", "elixhauser", "tonelli"),
    logit65 = -1.60, slope = 0.030, hf_boost = 0.70, min_age = 20,
    known_p = 0.90, dad_sens = 0.85, dad_fp = 0.008,
    icd10 = c("E11.9", "E10.9", "E11.2", "E11.4", "E11.5", "E10.2"),
    icd10_w = c(7, 0.8, 1.6, 1.0, 1.0, 0.4),
    icd9  = c("250"),
    claims_per_year = 2.4, readmit_beta = 0.20),
  HYPOTHYROID = cond(
    "Hypothyroidism", c("elixhauser", "tonelli"),
    logit65 = -2.10, slope = 0.020, min_age = 18,
    known_p = 0.88, dad_sens = 0.52,
    icd10 = c("E03.9", "E03.8"), icd10_w = c(6, 1),
    icd9  = c("244"),
    claims_per_year = 1.2),
  OBESITY = cond(
    "Obesity", "elixhauser",
    logit65 = -1.10, slope = -0.005, hf_boost = 0.45, min_age = 18,
    known_p = 0.55, dad_sens = 0.30, dad_fp = 0.006,
    icd10 = c("E66.9"), icd9 = c("278"),
    claims_per_year = 0.8),
  DYSLIPID = cond(
    "Dyslipidaemia", character(0),
    logit65 = -0.60, slope = 0.020, hf_boost = 0.30, min_age = 25,
    known_p = 0.80, dad_sens = 0.42,
    icd10 = c("E78.5"), icd9 = c("272"),
    claims_per_year = 1.4, claims_fp_per_year = 0.09),

  # --- renal -------------------------------------------------------------
  CKD = cond(
    # Charlson/Elixhauser renal failure use N18-N19 and dialysis codes;
    # Tonelli's administrative fallback (Ronksley 2012) is the whole of
    # N00-N23, which also sweeps in acute kidney injury from the acute pool.
    "Chronic kidney disease", c("charlson", "elixhauser", "tonelli"),
    logit65 = -2.80, slope = 0.060, hf_boost = 1.10, min_age = 30,
    known_p = 0.78, dad_sens = 0.66,
    icd10 = c("N18.3", "N18.4", "N18.5", "N18.9", "Z99.2"),
    icd10_w = c(6, 2.5, 1.2, 1.5, 0.5),
    icd9  = c("585"),
    claims_per_year = 1.6,
    provider = c(GP = 0.55, SPEC = 0.45), readmit_beta = 0.55),

  # --- digestive / liver -------------------------------------------------
  PEPTIC_ULCER = cond(
    "Peptic ulcer disease", c("charlson", "elixhauser", "tonelli"),
    logit65 = -3.90, slope = 0.030, min_age = 25,
    known_p = 0.80, dad_sens = 0.55,
    icd10 = c("K25.7", "K25.9", "K27.9"), icd10_w = c(3, 3, 1),
    icd9  = c("531", "533"), icd9_w = c(3, 1),
    claims_per_year = 0.9),
  LIVER_MILD = cond(
    "Chronic liver disease, mild", c("charlson", "elixhauser"),
    logit65 = -3.90, slope = 0.010, min_age = 20,
    known_p = 0.75, dad_sens = 0.60,
    icd10 = c("K70.0", "K73.9", "K76.0"), icd10_w = c(2, 2, 3),
    icd9  = c("571"),
    claims_per_year = 1.0, provider = c(GP = 0.6, SPEC = 0.4)),
  CIRRHOSIS = cond(
    "Cirrhosis", c("charlson", "elixhauser", "tonelli"),
    logit65 = -4.80, slope = 0.010, min_age = 25,
    known_p = 0.85, dad_sens = 0.72,
    icd10 = c("K70.3", "K74.6"), icd10_w = c(2, 3),
    icd9  = c("571"),
    claims_per_year = 1.4, provider = c(GP = 0.5, SPEC = 0.5),
    readmit_beta = 0.25),
  LIVER_SEVERE = cond(
    # Charlson's moderate/severe liver disease: decompensation on top of
    # cirrhosis. Tonelli asks for cirrhosis AND decompensation codes.
    "Decompensated liver disease", c("charlson", "tonelli"),
    logit65 = -5.60, slope = 0.005, min_age = 25,
    known_p = 0.90, dad_sens = 0.78,
    icd10 = c("K76.7", "I85.0", "K72.9"), icd10_w = c(3, 2, 1),
    icd9  = c("572", "456"), icd9_w = c(3, 1),
    claims_per_year = 1.2, provider = c(GP = 0.4, SPEC = 0.6),
    readmit_beta = 0.30),
  HEP_B = cond(
    "Chronic viral hepatitis B", "tonelli",
    logit65 = -5.50, slope = -0.005, min_age = 18,
    known_p = 0.70, dad_sens = 0.55,
    icd10 = c("B18.1"), icd9 = c("070"),
    claims_per_year = 1.2, provider = c(GP = 0.5, SPEC = 0.5)),
  IBD = cond(
    "Inflammatory bowel disease", "tonelli",
    logit65 = -4.40, slope = -0.005, min_age = 15,
    known_p = 0.88, dad_sens = 0.70,
    icd10 = c("K50.9", "K51.9"), icd10_w = c(1, 1),
    icd9  = c("555", "556"),
    claims_per_year = 1.6, provider = c(GP = 0.45, SPEC = 0.55)),
  IBS = cond(
    "Irritable bowel syndrome", "tonelli",
    logit65 = -3.50, slope = -0.010, min_age = 15,
    known_p = 0.70, dad_sens = 0.40,
    icd10 = c("K58.9", "K58.0"), icd10_w = c(3, 1),
    icd9  = c("564"),
    claims_per_year = 1.2),
  CONSTIPATION = cond(
    # Tonelli's severe constipation. Note the ICD-9 collision: at three
    # digits, 564 is shared with irritable bowel syndrome, so the claims
    # side of these two algorithms cannot be separated here. 560 can.
    "Severe constipation", "tonelli",
    logit65 = -3.90, slope = 0.045, min_age = 18,
    known_p = 0.65, dad_sens = 0.45,
    icd10 = c("K59.0", "K56.0"), icd10_w = c(6, 1),
    icd9  = c("560", "564"), icd9_w = c(3, 1),
    claims_per_year = 1.0),

  # --- neurological ------------------------------------------------------
  DEMENTIA = cond(
    "Dementia", c("charlson", "tonelli"),
    logit65 = -3.60, slope = 0.110, hf_boost = 0.20, min_age = 55,
    known_p = 0.70, dad_sens = 0.62,
    icd10 = c("F03", "G30.9", "F05.1"), icd10_w = c(5, 3, 1),
    icd9  = c("290", "331"), icd9_w = c(6, 4),
    claims_per_year = 1.4, provider = c(GP = 0.75, SPEC = 0.25),
    readmit_beta = 0.25),
  PARKINSON = cond(
    "Parkinson's disease", c("elixhauser", "tonelli"),
    logit65 = -4.80, slope = 0.070, min_age = 45,
    known_p = 0.90, dad_sens = 0.70,
    icd10 = c("G20"), icd9 = c("332"),
    claims_per_year = 1.8, provider = c(GP = 0.5, SPEC = 0.5)),
  MS = cond(
    "Multiple sclerosis", c("elixhauser", "tonelli"),
    logit65 = -5.50, slope = -0.020, min_age = 18,
    known_p = 0.92, dad_sens = 0.72,
    icd10 = c("G35"), icd9 = c("340"),
    claims_per_year = 2.0, provider = c(GP = 0.4, SPEC = 0.6)),
  EPILEPSY = cond(
    "Epilepsy", c("elixhauser", "tonelli"),
    logit65 = -4.20, slope = 0.010, min_age = 5,
    known_p = 0.85, dad_sens = 0.66,
    icd10 = c("G40.9", "G40.2"), icd10_w = c(3, 1),
    icd9  = c("345"),
    claims_per_year = 1.6, provider = c(GP = 0.5, SPEC = 0.5)),
  HEMIPLEGIA = cond(
    "Hemiplegia or paraplegia", c("charlson", "elixhauser"),
    logit65 = -4.40, slope = 0.035, min_age = 18,
    known_p = 0.92, dad_sens = 0.74,
    icd10 = c("G81.9", "G82.2", "G82.1"), icd10_w = c(5, 2, 1),
    icd9  = c("342", "344"), icd9_w = c(3, 2),
    claims_per_year = 1.4, provider = c(GP = 0.55, SPEC = 0.45)),
  CHRONIC_PAIN = cond(
    "Chronic pain", "tonelli",
    logit65 = -2.00, slope = 0.010, min_age = 18,
    known_p = 0.70, dad_sens = 0.30,
    icd10 = c("M54.5", "M51.1", "M79.7", "M47.8"), icd10_w = c(5, 2, 2, 1),
    icd9  = c("724", "722", "729"), icd9_w = c(4, 2, 3),
    claims_per_year = 2.2, claims_fp_per_year = 0.12),

  # --- mental health / substance use -------------------------------------
  DEPRESSION = cond(
    "Depression", c("elixhauser", "tonelli"),
    logit65 = -1.90, slope = -0.005, hf_boost = 0.35, min_age = 14,
    known_p = 0.65, dad_sens = 0.38,
    icd10 = c("F32.9", "F33.9", "F41.2"), icd10_w = c(5, 3, 1),
    icd9  = c("311", "296"), icd9_w = c(8, 2),
    claims_per_year = 1.8, claims_fp_per_year = 0.10),
  SCHIZOPHRENIA = cond(
    "Schizophrenia", c("elixhauser", "tonelli"),
    logit65 = -4.40, slope = -0.020, min_age = 16,
    known_p = 0.90, dad_sens = 0.72,
    icd10 = c("F20.9", "F25.9"), icd10_w = c(4, 1),
    icd9  = c("295"),
    claims_per_year = 2.4, provider = c(GP = 0.5, SPEC = 0.5)),
  ALCOHOL = cond(
    "Alcohol misuse", c("elixhauser", "tonelli"),
    logit65 = -3.00, slope = -0.030, min_age = 18,
    known_p = 0.60, dad_sens = 0.58,
    icd10 = c("F10.2", "F10.1", "K29.2"), icd10_w = c(5, 2, 1),
    icd9  = c("303", "305"), icd9_w = c(3, 1),
    claims_per_year = 1.2, readmit_beta = 0.20),
  DRUG_ABUSE = cond(
    "Drug abuse", "elixhauser",
    logit65 = -3.90, slope = -0.045, min_age = 16,
    known_p = 0.60, dad_sens = 0.58,
    icd10 = c("F19.2", "F11.2"), icd10_w = c(3, 2),
    icd9  = c("304", "305"), icd9_w = c(3, 1),
    claims_per_year = 1.2),
  SMOKING = cond(
    "Tobacco use", character(0),
    logit65 = -1.40, slope = -0.020, hf_boost = 0.30, min_age = 18,
    known_p = 0.70, dad_sens = 0.45,
    icd10 = c("F17.2", "Z72.0"), icd10_w = c(1.2, 1),
    # 305 at three digits is shared with alcohol and drug misuse, so
    # smoking claims use the personal-history code instead.
    icd9  = c("V15"),
    claims_per_year = 0.6),

  # --- musculoskeletal / skin --------------------------------------------
  RHEUM = cond(
    "Rheumatoid arthritis / connective tissue disease",
    c("charlson", "elixhauser", "tonelli"),
    logit65 = -3.90, slope = 0.015, min_age = 20,
    known_p = 0.88, dad_sens = 0.62,
    icd10 = c("M05.9", "M06.9", "M32.9"), icd10_w = c(3, 4, 1),
    icd9  = c("714", "710"), icd9_w = c(4, 1),
    claims_per_year = 1.8, provider = c(GP = 0.45, SPEC = 0.55)),
  OSTEOARTH = cond(
    "Osteoarthritis", character(0),
    logit65 = -1.30, slope = 0.030, min_age = 40,
    known_p = 0.75, dad_sens = 0.40,
    icd10 = c("M17.1", "M16.1"), icd10_w = c(6, 4),
    icd9  = c("715"),
    claims_per_year = 1.2, claims_fp_per_year = 0.09),
  PSORIASIS = cond(
    "Psoriasis", "tonelli",
    logit65 = -3.90, slope = -0.005, min_age = 15,
    known_p = 0.80, dad_sens = 0.35,
    icd10 = c("L40.0", "L40.9"), icd10_w = c(3, 2),
    icd9  = c("696"),
    claims_per_year = 1.2, provider = c(GP = 0.5, SPEC = 0.5)),

  # --- cancer ------------------------------------------------------------
  # Tonelli treats these as remitting after 5 years; here they are
  # permanent once assigned.
  CANCER_LUNG = cond(
    "Cancer, lung", c("charlson", "elixhauser", "tonelli"),
    logit65 = -4.80, slope = 0.030, min_age = 35,
    known_p = 0.95, dad_sens = 0.80,
    icd10 = c("C34.9", "C34.1"), icd10_w = c(3, 1),
    icd9  = c("162"),
    claims_per_year = 2.4, provider = c(GP = 0.35, SPEC = 0.65),
    readmit_beta = 0.30),
  CANCER_COLORECTAL = cond(
    "Cancer, colorectal", c("charlson", "elixhauser", "tonelli"),
    logit65 = -4.60, slope = 0.030, min_age = 35,
    known_p = 0.95, dad_sens = 0.78,
    icd10 = c("C18.9", "C20"), icd10_w = c(3, 2),
    icd9  = c("153", "154"), icd9_w = c(3, 2),
    claims_per_year = 2.0, provider = c(GP = 0.35, SPEC = 0.65)),
  CANCER_BREAST = cond(
    "Cancer, breast", c("charlson", "elixhauser", "tonelli"),
    logit65 = -3.50, slope = 0.020, min_age = 30, sex = "F",
    known_p = 0.95, dad_sens = 0.72,
    icd10 = c("C50.9"), icd9 = c("174"),
    claims_per_year = 1.8, provider = c(GP = 0.4, SPEC = 0.6)),
  CANCER_PROSTATE = cond(
    "Cancer, prostate", c("charlson", "elixhauser", "tonelli"),
    logit65 = -3.20, slope = 0.045, min_age = 45, sex = "M",
    known_p = 0.92, dad_sens = 0.70,
    icd10 = c("C61"), icd9 = c("185"),
    claims_per_year = 1.6, provider = c(GP = 0.4, SPEC = 0.6)),
  CANCER_METASTATIC = cond(
    "Cancer, metastatic", c("charlson", "elixhauser", "tonelli"),
    logit65 = -4.60, slope = 0.030, min_age = 30,
    known_p = 0.97, dad_sens = 0.84,
    icd10 = c("C78.0", "C79.5", "C78.7"), icd10_w = c(3, 2, 2),
    icd9  = c("197", "198"), icd9_w = c(3, 2),
    claims_per_year = 3.0, provider = c(GP = 0.3, SPEC = 0.7),
    readmit_beta = 0.50),
  LYMPHOMA = cond(
    "Cancer, lymphoma or myeloma", c("charlson", "elixhauser", "tonelli"),
    logit65 = -5.10, slope = 0.030, min_age = 20,
    known_p = 0.95, dad_sens = 0.80,
    icd10 = c("C83.3", "C85.9", "C90.0"), icd10_w = c(2, 2, 2),
    icd9  = c("200", "202", "203"), icd9_w = c(2, 2, 2),
    claims_per_year = 2.4, provider = c(GP = 0.3, SPEC = 0.7),
    readmit_beta = 0.25),
  HIV = cond(
    "HIV / AIDS", c("charlson", "elixhauser"),
    logit65 = -5.50, slope = -0.035, min_age = 18,
    known_p = 0.80, dad_sens = 0.70,
    icd10 = c("B24", "B20.9"), icd10_w = c(2, 1),
    icd9  = c("042"),
    claims_per_year = 3.0, provider = c(GP = 0.4, SPEC = 0.6))
)

CONDITION_NAMES <- names(CHRONIC)
ALGORITHMS <- c("charlson", "elixhauser", "tonelli")

# ---------------------------------------------------------------------------
# Acute / episode-specific ICD-10-CA pool
# ---------------------------------------------------------------------------
# Chronic-condition codes are NOT in this pool: those are emitted only
# through the latent-truth mechanism above, so a hypertension code never
# appears for a patient without hypertension except at the stated false
# positive rate. The `feeds` column notes the published comorbidity
# categories an acute code will nonetheless trigger -- several Elixhauser
# categories are acute by nature and are populated only from here.
acute <- function(code, desc, min_age = 0, max_age = 105, sex = NA, weight = 1,
                  feeds = "") {
  data.frame(code = code, desc = desc, min_age = min_age, max_age = max_age,
             sex = sex, weight = weight, feeds = feeds,
             stringsAsFactors = FALSE)
}

ACUTE_POOL <- do.call(rbind, list(
  # --- circulatory ---------------------------------------------------
  acute("I21.0", "Acute transmural MI of anterior wall",            35, 105, NA, 0.9,
        "Charlson/Tonelli myocardial infarction"),
  acute("I21.1", "Acute transmural MI of inferior wall",            35, 105, NA, 0.9,
        "Charlson/Tonelli myocardial infarction"),
  acute("I21.4", "Acute subendocardial myocardial infarction",      35, 105, NA, 1.8,
        "Charlson/Tonelli myocardial infarction"),
  acute("I21.9", "Acute myocardial infarction, unspecified",        35, 105, NA, 0.9,
        "Charlson/Tonelli myocardial infarction"),
  acute("I26.9", "Pulmonary embolism",                              20, 105, NA, 1.0,
        "Elixhauser pulmonary circulation disorders"),
  acute("I63.9", "Cerebral infarction, unspecified",                40, 105, NA, 1.2,
        "Charlson cerebrovascular; Tonelli stroke/TIA"),
  acute("I64",   "Stroke, not specified as haemorrhage or infarct", 40, 105, NA, 0.5,
        "Charlson cerebrovascular; Tonelli stroke/TIA"),
  acute("I61.9", "Intracerebral haemorrhage, unspecified",          40, 105, NA, 0.4,
        "Charlson cerebrovascular; Tonelli stroke/TIA"),
  acute("G45.9", "Transient cerebral ischaemic attack",             40, 105, NA, 0.8,
        "Charlson cerebrovascular; Tonelli stroke/TIA"),
  acute("I65.2", "Occlusion and stenosis of carotid artery",        45, 105, NA, 0.6),
  acute("I80.2", "Phlebitis and thrombophlebitis, deep vessels",    20, 105, NA, 0.7),
  acute("I95.1", "Orthostatic hypotension",                         50, 105, NA, 0.6),
  # --- respiratory ---------------------------------------------------
  acute("J18.9", "Pneumonia, unspecified",                           0, 105, NA, 4.0),
  acute("J96.0", "Acute respiratory failure",                        0, 105, NA, 1.8),
  acute("J90",   "Pleural effusion, not elsewhere classified",      20, 105, NA, 1.2),
  acute("J69.0", "Pneumonitis due to food and vomit",               50, 105, NA, 0.8),
  acute("U07.1", "COVID-19, virus identified",                       0, 105, NA, 1.2),
  # --- endocrine / metabolic -----------------------------------------
  acute("E87.1", "Hypo-osmolality and hyponatraemia",                0, 105, NA, 2.2,
        "Elixhauser fluid and electrolyte disorders"),
  acute("E87.6", "Hypokalaemia",                                     0, 105, NA, 1.8,
        "Elixhauser fluid and electrolyte disorders"),
  acute("E86",   "Volume depletion",                                 0, 105, NA, 2.0,
        "Elixhauser fluid and electrolyte disorders"),
  acute("E46",   "Unspecified protein-energy malnutrition",         50, 105, NA, 1.0,
        "Elixhauser weight loss"),
  acute("R63.4", "Abnormal weight loss",                            40, 105, NA, 0.8,
        "Elixhauser weight loss"),
  acute("E05.9", "Thyrotoxicosis, unspecified",                     18, 105, NA, 0.5),
  # --- renal ---------------------------------------------------------
  acute("N17.9", "Acute kidney failure, unspecified",                0, 105, NA, 2.5,
        "Tonelli CKD (N00-N23 includes AKI)"),
  acute("N39.0", "Urinary tract infection, site not specified",      0, 105, NA, 3.0),
  acute("N40",   "Hyperplasia of prostate",                         50, 105, "M", 1.8),
  # --- digestive ------------------------------------------------------
  acute("K92.2", "Gastrointestinal haemorrhage, unspecified",       18, 105, NA, 1.2),
  acute("K80.2", "Calculus of gallbladder",                         18, 105, NA, 1.2),
  acute("K57.3", "Diverticular disease of large intestine",         40, 105, NA, 1.0),
  acute("K21.9", "Gastro-oesophageal reflux disease",               18, 105, NA, 2.0),
  # --- infection ------------------------------------------------------
  acute("A41.9", "Sepsis, unspecified",                              0, 105, NA, 2.0),
  acute("A04.7", "Enterocolitis due to Clostridium difficile",      18, 105, NA, 0.7),
  acute("B95.6", "Staphylococcus aureus as cause of disease",        0, 105, NA, 0.6),
  # --- musculoskeletal / injury ---------------------------------------
  acute("S72.0", "Fracture of neck of femur",                       50, 105, NA, 1.8),
  acute("S06.0", "Concussion",                                       5, 105, NA, 0.8),
  acute("M81.9", "Osteoporosis, unspecified",                       50, 105, NA, 1.6),
  acute("W19",   "Unspecified fall",                                40, 105, NA, 2.0),
  acute("W01",   "Fall on same level from slipping or tripping",    40, 105, NA, 1.5),
  # --- neuro / mental health ------------------------------------------
  acute("F05.9", "Delirium, unspecified",                           50, 105, NA, 1.8),
  acute("F41.9", "Anxiety disorder, unspecified",                   12, 105, NA, 1.6),
  acute("G47.3", "Sleep apnoea",                                    25, 105, NA, 1.5),
  # --- blood ----------------------------------------------------------
  acute("D62",   "Acute posthaemorrhagic anaemia",                   0, 105, NA, 1.0),
  acute("D50.0", "Iron deficiency anaemia secondary to blood loss",  0, 105, NA, 1.2,
        "Elixhauser blood loss anaemia"),
  acute("D50.9", "Iron deficiency anaemia, unspecified",             0, 105, NA, 1.8,
        "Elixhauser deficiency anaemia"),
  acute("D64.9", "Anaemia, unspecified",                             0, 105, NA, 2.0),
  acute("D68.3", "Haemorrhagic disorder due to anticoagulants",     40, 105, NA, 0.8,
        "Elixhauser coagulopathy"),
  acute("D69.6", "Thrombocytopenia, unspecified",                    0, 105, NA, 0.8,
        "Elixhauser coagulopathy"),
  # --- symptoms / signs -------------------------------------------------
  acute("R07.4", "Chest pain, unspecified",                         18, 105, NA, 1.5),
  acute("R06.0", "Dyspnoea",                                        18, 105, NA, 1.8),
  acute("R55",   "Syncope and collapse",                            18, 105, NA, 1.2),
  acute("R50.9", "Fever, unspecified",                               0, 105, NA, 1.2),
  acute("R33",   "Retention of urine",                              40, 105, NA, 1.0),
  # --- obstetric / newborn ----------------------------------------------
  acute("O80",   "Single spontaneous delivery",                     15,  50, "F", 6.0),
  acute("O82",   "Single delivery by caesarean section",            15,  50, "F", 3.0),
  acute("O70.1", "Second degree perineal laceration",               15,  50, "F", 1.5),
  acute("O24.4", "Diabetes mellitus arising in pregnancy",          15,  50, "F", 1.2),
  acute("O14.0", "Moderate pre-eclampsia",                          15,  50, "F", 0.8),
  acute("Z37.0", "Single live birth",                               15,  50, "F", 4.0),
  acute("Z38.0", "Singleton, born in hospital",                      0,   0, NA, 6.0),
  acute("P07.3", "Other preterm infants",                            0,   0, NA, 1.5),
  acute("P59.9", "Neonatal jaundice, unspecified",                   0,   0, NA, 2.0),
  acute("P92.5", "Neonatal difficulty in feeding at breast",         0,   0, NA, 1.0),
  # --- factors influencing health status ----------------------------------
  acute("Z51.5", "Palliative care",                                 40, 105, NA, 1.0),
  acute("Z95.0", "Presence of cardiac pacemaker",                   50, 105, NA, 0.8,
        "Elixhauser cardiac arrhythmias"),
  acute("Z92.1", "Personal history of long-term anticoagulant use", 40, 105, NA, 1.5)
))

# Codes that get an extra push as secondary diagnoses when the most
# responsible diagnosis is heart failure.
HF_SECONDARY_BOOST <- c("E87.1" = 2.5, "E87.6" = 2.0, "N17.9" = 2.5,
                        "D64.9" = 2.0, "Z92.1" = 2.0, "R06.0" = 2.5,
                        "Z95.0" = 2.5, "G47.3" = 1.8, "J90" = 2.5,
                        "J96.0" = 1.8, "Z51.5" = 1.5, "D50.9" = 1.5)

# Reasons an admission happens for someone with heart failure when the most
# responsible diagnosis is not heart failure itself.
HF_ALT_MRDX <- c("J18.9", "J44.1", "N17.9", "I21.4", "A41.9", "E87.1",
                 "R06.0", "S72.0", "U07.1", "N39.0", "I63.9", "K92.2")

# ---------------------------------------------------------------------------
# ICD-9 pool for the claims file (three digits, no decimal)
# ---------------------------------------------------------------------------
# Acute / non-chronic reasons for an office or ED visit.
ICD9_ACUTE <- data.frame(
  code = c("460", "465", "466", "486", "599", "780", "786", "789", "784",
           "719", "692", "477", "530", "535", "558", "787", "781",
           "V70", "V58", "V72", "V06"),
  weight = c(2.0, 2.5, 1.8, 1.5, 2.2, 3.0, 2.6, 1.6, 1.8,
             1.5, 1.4, 1.2, 1.3, 1.0, 1.2, 1.4, 1.0,
             4.0, 2.5, 1.5, 1.2),
  stringsAsFactors = FALSE)

# Heart failure in ICD-9 (Tonelli chronic heart failure: 428 among others).
ICD9_HF <- c("428", "425")
ICD9_HF_W <- c(9, 1)

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

wsample1 <- function(x, w) {
  if (length(x) == 1L) return(x)
  x[sample.int(length(x), 1L, prob = w)]
}

draw_age <- function() {
  # Age mixture roughly shaped like an acute-care inpatient population.
  u <- runif(1)
  if (u < 0.06) return(0L)
  if (u < 0.11) return(sample(1:17, 1))
  if (u < 0.30) return(sample(18:49, 1))
  if (u < 0.58) return(sample(50:69, 1))
  if (u < 0.88) return(sample(70:84, 1))
  sample(85:101, 1)
}

draw_sex <- function() {
  # DAD sex codes: M, F, O (other), U (unknown). O and U are rare but real,
  # and an analysis that silently drops them is making a decision it should
  # be making out loud.
  sample(c("M", "F", "O", "U"), 1, prob = c(48.5, 49.5, 1.0, 1.0))
}

draw_dx_count <- function(age) {
  mean_extra <- if (age == 0) 1.2 else if (age < 18) 1.8 else if (age < 50) 2.8
                else if (age < 70) 4.0 else if (age < 85) 5.0 else 6.0
  n <- 1L + as.integer(rexp(1, rate = 1 / mean_extra))
  max(1L, min(MAX_DX, n))
}

draw_los_hours <- function(age, is_hf, n_dx) {
  los_days <- rlnorm(1, meanlog = 0.85, sdlog = 0.85)
  if (is_hf) los_days <- los_days * 1.6
  if (age >= 80) los_days <- los_days * 1.35 else if (age >= 65) los_days <- los_days * 1.15
  los_days <- los_days * (1 + 0.02 * n_dx)
  los_days <- min(los_days, 180)
  max(4, los_days * 24) + runif(1, -3, 3)
}

# ---------------------------------------------------------------------------
# Phase 1: patients (latent truth)
# ---------------------------------------------------------------------------

make_patients <- function(n_patients, hf_patient_rate, start) {
  age0 <- vapply(seq_len(n_patients), function(i) draw_age(), integer(1))
  sex  <- vapply(seq_len(n_patients), function(i) draw_sex(), character(1))

  # Birth date, kept internally only. The workshop notes that birthdate is
  # identifiable and normally is not released -- the analyst computes age and
  # releases that instead, which is exactly what happens here.
  birth <- as.Date(start) - (age0 * 365.25 + runif(n_patients, 0, 365))

  # Heart failure is a patient attribute, not an encounter attribute: once a
  # patient has it they keep it. Strongly age dependent.
  hf_p <- inv_logit(-1.2 + 0.055 * (age0 - 70))
  hf_p[age0 < 40] <- 0
  hf_p <- pmin(hf_p * (hf_patient_rate / max(mean(hf_p), 1e-9)), 0.97)
  is_hf <- runif(n_patients) < hf_p

  conds <- matrix(FALSE, nrow = n_patients, ncol = length(CONDITION_NAMES),
                  dimnames = list(NULL, CONDITION_NAMES))
  known <- conds
  for (nm in CONDITION_NAMES) {
    cc <- CHRONIC[[nm]]
    p <- inv_logit(cc$logit65 + cc$slope * (age0 - 65) + cc$hf_boost * is_hf)
    p[age0 < cc$min_age] <- 0
    if (!is.na(cc$sex)) p[sex != cc$sex] <- 0
    conds[, nm] <- runif(n_patients) < p
    # Being diagnosed is a second lottery. A condition nobody has diagnosed
    # cannot appear in any administrative source, so it puts a ceiling on
    # the sensitivity of every algorithm the workshop will write -- while
    # still affecting the patient's risk of readmission.
    known[, nm] <- conds[, nm] & (runif(n_patients) < cc$known_p)
  }
  # Patients with heart failure are under close follow-up, so more of their
  # comorbidity has been picked up -- but not all of it.
  for (nm in c("HTN", "CKD", "AFIB", "IHD", "COPD")) {
    catchup <- is_hf & conds[, nm] & !known[, nm] & (runif(n_patients) < 0.6)
    known[catchup, nm] <- TRUE
  }

  list(patient_id = sprintf("P%06d", seq_len(n_patients)),
       age0 = age0, sex = sex, birth = birth, is_hf = is_hf,
       conds = conds, known = known)
}

# ---------------------------------------------------------------------------
# Phase 2: encounters (dates, transfers, deaths)
# ---------------------------------------------------------------------------

READMIT_BETA <- vapply(CHRONIC, function(cc) cc$readmit_beta, numeric(1))

gen_patient_encounters <- function(i, pts, start_dt, end_dt, base_gap_days,
                                   transfer_rate) {
  is_hf <- pts$is_hf[i]
  cvec <- pts$conds[i, ]
  # Readmission risk follows the underlying conditions, not the coded ones,
  # which is why the fitted hazard ratios come out attenuated.
  # Additive in the conditions that carry a readmission effect. Only the
  # conditions the workshop is likely to model carry one, so the sum stays
  # bounded and each one's marginal effect is still there to be estimated.
  lh_comorb <- sum(READMIT_BETA[cvec])
  has_ckd <- isTRUE(cvec[["CKD"]])
  has_dementia <- isTRUE(cvec[["DEMENTIA"]])
  has_mets <- isTRUE(cvec[["CANCER_METASTATIC"]])

  admit <- start_dt + runif(1, 0, as.numeric(difftime(end_dt, start_dt, units = "secs")))
  out <- list()
  n <- 0L
  is_transfer_in <- FALSE

  repeat {
    n <- n + 1L
    age <- max(0L, as.integer(floor(as.numeric(as.Date(admit) - pts$birth[i]) / 365.25)))

    n_dx_hint <- draw_dx_count(age)
    los <- draw_los_hours(age, is_hf, n_dx_hint)
    if (is_transfer_in) los <- max(12, los * 0.7)   # the receiving stay
    disch <- admit + los * 3600

    admit_cat <- if (age == 0 && n == 1L) "N"
                 else if (is_transfer_in) "U"
                 else if (runif(1) < 0.13 && !is_hf) "L" else "U"

    # In-hospital death: a competing risk for readmission, and the reason a
    # readmission denominator has to exclude patients discharged dead.
    p_death <- inv_logit(-3.6 + 0.05 * (age - 70) + 0.6 * is_hf +
                         0.4 * has_ckd + 0.3 * has_dementia + 0.9 * has_mets)
    died <- runif(1) < p_death

    transferred <- !died && runif(1) < transfer_rate * (1 + 0.02 * max(0, age - 60))

    disp <- if (died) "07"
            else if (transferred) "01"
            else sample(c("05", "04", "02", "03", "06"), 1,
                        prob = c(0.55, 0.25,
                                 if (age >= 75) 0.20 else 0.06,
                                 0.05, 0.03))

    out[[n]] <- list(admit = admit, disch = disch, age = age,
                     admit_cat = admit_cat, disp = disp,
                     is_transfer_in = is_transfer_in, n_dx_hint = n_dx_hint)

    if (died) return(list(enc = out, death = disch))
    if (n >= MAX_ENC_PER_PATIENT) return(list(enc = out, death = NA))

    if (transferred) {
      # Acute-to-acute transfer: a new abstract at a different facility,
      # hours after the first discharge. Looks like a readmission unless you
      # build episodes of care.
      admit <- disch + runif(1, 0.5, 8) * 3600
      is_transfer_in <- TRUE
    } else {
      is_transfer_in <- FALSE
      lh <- 0.025 * (age - 70) + 0.60 * is_hf + lh_comorb
      # Two-component gap: an early-readmission component that a 30- or
      # 90-day outcome will pick up, and a long background component.
      p_early <- inv_logit(-2.7 + 0.7 * lh)
      gap <- if (runif(1) < p_early) max(0.6, rexp(1, rate = 1 / 45))
             else rexp(1, rate = exp(lh) / base_gap_days)
      admit <- disch + gap * SECS_PER_DAY
    }
    if (admit > end_dt) return(list(enc = out, death = NA))
  }
}

# ---------------------------------------------------------------------------
# Phase 3: diagnoses
# ---------------------------------------------------------------------------

eligible_acute <- function(age, sex) {
  ok <- ACUTE_POOL$min_age <= age & age <= ACUTE_POOL$max_age &
        (is.na(ACUTE_POOL$sex) | ACUTE_POOL$sex == sex)
  ACUTE_POOL[ok, , drop = FALSE]
}

build_diagnoses <- function(age, sex, hf_mrdx, kvec, p_code, n_total,
                            is_hf_patient) {
  pool <- eligible_acute(age, sex)
  if (nrow(pool) == 0)
    pool <- ACUTE_POOL[is.na(ACUTE_POOL$sex) & ACUTE_POOL$min_age <= age, , drop = FALSE]

  # --- most responsible diagnosis ---
  if (hf_mrdx) {
    mrdx <- if (isTRUE(kvec[["HTN"]]) && isTRUE(kvec[["CKD"]]) && runif(1) < 0.10)
              HF_HTN_CKD_CODE
            else if (isTRUE(kvec[["HTN"]]) && runif(1) < 0.14) HF_HTN_CODE
            else wsample1(HF_CODES, HF_CODE_W)
  } else if (is_hf_patient && age >= 40 && runif(1) < 0.55) {
    alt <- intersect(HF_ALT_MRDX, pool$code)
    mrdx <- if (length(alt)) sample(alt, 1) else wsample1(pool$code, pool$weight)
  } else {
    mrdx <- wsample1(pool$code, pool$weight)
  }
  codes <- mrdx

  # A patient with heart failure who is admitted for something else usually
  # still has it coded, in a secondary field. This is what separates
  # "admitted FOR heart failure" from "admitted WITH heart failure", and it
  # is why the choice of diagnosis field changes the cohort so much.
  if (is_hf_patient && !hf_mrdx && age >= 40 && runif(1) < 0.58)
    codes <- c(codes, wsample1(HF_CODES, HF_CODE_W))

  # --- chronic comorbidities, recorded with condition-specific sensitivity ---
  # p_code is precomputed per patient: sensitivity where the condition is
  # present and diagnosed, the false-positive rate where it is not.
  p <- p_code
  if (hf_mrdx) p <- pmin(0.98, p * 1.15)
  hit <- which(runif(length(p)) < p)
  for (k in hit) {
    nm <- CONDITION_NAMES[k]
    # The hypertension combination codes already carry these diagnoses.
    if (nm == "HTN" && mrdx %in% c(HF_HTN_CODE, HF_HTN_CKD_CODE)) next
    if (nm == "CKD" && mrdx == HF_HTN_CKD_CODE) next
    cc <- CHRONIC[[k]]
    pick <- wsample1(cc$icd10, cc$icd10_w)
    if (!(pick %in% codes)) codes <- c(codes, pick)
  }

  # --- fill the remaining slots with acute / episode codes ---
  n_total <- min(max(n_total, length(codes)), MAX_DX)
  cand <- pool$code
  w <- pool$weight
  if (hf_mrdx) {
    boost <- HF_SECONDARY_BOOST[cand]
    boost[is.na(boost)] <- 1
    w <- w * as.numeric(boost)
  }
  guard <- 0L
  while (length(codes) < n_total && guard < 400L) {
    guard <- guard + 1L
    pick <- wsample1(cand, w)
    if (!(pick %in% codes)) codes <- c(codes, pick)
  }
  if (length(codes) > MAX_DX) codes <- c(codes[1], sample(codes[-1], MAX_DX - 1))

  # Secondary diagnoses are not filed in any meaningful order; shuffle them
  # so students cannot cheat by reading position instead of searching all 25.
  if (length(codes) > 2) codes <- c(codes[1], sample(codes[-1]))
  codes
}

# ---------------------------------------------------------------------------
# Phase 4: physician claims (ICD-9), built vectorised over patients
# ---------------------------------------------------------------------------

gen_claims <- function(pts, win_lo, win_hi, years, stay_admit, stay_disch,
                       stay_pt) {
  n <- length(pts$patient_id)
  active <- years > 0.02

  pt_idx <- integer(0); dx1 <- character(0)
  prov_gp <- numeric(0)

  add_stream <- function(lambda, codes, code_w, p_gp) {
    cnt <- rpois(n, pmax(0, lambda))
    cnt[!active] <- 0L
    if (sum(cnt) == 0) return(invisible(NULL))
    idx <- rep.int(seq_len(n), cnt)
    pt_idx <<- c(pt_idx, idx)
    dx1 <<- c(dx1, sample(codes, length(idx), replace = TRUE, prob = code_w))
    prov_gp <<- c(prov_gp, rep.int(p_gp, length(idx)))
    invisible(NULL)
  }

  # 1. Condition-specific follow-up claims. This is the stream that makes
  #    "k claims within a window" definitions work.
  for (nm in CONDITION_NAMES) {
    cc <- CHRONIC[[nm]]
    has <- pts$known[, nm]
    lambda <- ifelse(has, cc$claims_per_year, cc$claims_fp_per_year) * years
    add_stream(lambda, cc$icd9, cc$icd9_w, unname(cc$provider["GP"]))
  }

  # 2. Heart-failure follow-up claims.
  add_stream(ifelse(pts$is_hf, 2.6, 0.04) * years, ICD9_HF, ICD9_HF_W, 0.55)

  # 3. Everything else: acute visits, checkups, minor complaints.
  base_rate <- 1.6 + 0.035 * pmax(0, pts$age0 - 40) + 0.25 * rowSums(pts$known)
  add_stream(base_rate * years, ICD9_ACUTE$code, ICD9_ACUTE$weight, 0.85)

  m <- length(pt_idx)
  if (m == 0) return(NULL)

  # Service dates: uniform inside each patient's observable window.
  t <- win_lo[pt_idx] + runif(m) * (win_hi[pt_idx] - win_lo[pt_idx])

  # Secondary diagnoses come from the patient's own conditions, or from the
  # acute pool. A flattened per-patient pool keeps this vectorised.
  pool_list <- lapply(seq_len(n), function(i) {
    cs <- CONDITION_NAMES[pts$known[i, ]]
    p <- unlist(lapply(cs, function(nm) CHRONIC[[nm]]$icd9), use.names = FALSE)
    if (pts$is_hf[i]) p <- c(p, "428")
    p
  })
  pool_len <- vapply(pool_list, length, integer(1))
  pool_flat <- unlist(pool_list, use.names = FALSE)
  pool_start <- cumsum(c(0L, pool_len))[seq_len(n)]

  draw_extra <- function(sel) {
    out <- character(length(sel))
    len <- pool_len[pt_idx[sel]]
    from_pool <- len > 0 & runif(length(sel)) < 0.6
    if (any(from_pool)) {
      j <- which(from_pool)
      out[j] <- pool_flat[pool_start[pt_idx[sel[j]]] +
                          1L + floor(runif(length(j)) * len[j])]
    }
    if (any(!from_pool)) {
      j <- which(!from_pool)
      out[j] <- sample(ICD9_ACUTE$code, length(j), replace = TRUE,
                       prob = ICD9_ACUTE$weight)
    }
    out
  }

  dx2 <- rep("", m)
  has2 <- runif(m) < 0.45
  if (any(has2)) dx2[has2] <- draw_extra(which(has2))
  dx2[dx2 == dx1] <- ""

  dx3 <- rep("", m)
  has3 <- has2 & dx2 != "" & runif(m) < 0.30
  if (any(has3)) dx3[has3] <- draw_extra(which(has3))
  dx3[dx3 == dx1 | dx3 == dx2] <- ""
  # Keep the diagnosis fields left-filled.
  move <- dx2 == "" & dx3 != ""
  dx2[move] <- dx3[move]; dx3[move] <- ""

  # A claim dated inside one of the patient's inpatient stays is billed
  # from hospital.
  loc <- sample(c("OFFICE", "ED", "LTC", "HOME"), m, replace = TRUE,
                prob = c(0.78, 0.15, 0.05, 0.02))
  if (length(stay_admit)) {
    ord <- order(stay_admit)
    sa <- stay_admit[ord]; sd <- stay_disch[ord]; sp <- stay_pt[ord]
    k <- findInterval(t, sa)
    inpt <- k > 0
    inpt[inpt] <- sp[k[inpt]] == pt_idx[inpt] & t[inpt] <= sd[k[inpt]]
    loc[inpt] <- "INPT"
  }

  prov <- ifelse(runif(m) < prov_gp, "GP", "SPEC")

  o <- order(t)
  data.frame(pt = pt_idx[o],
             date = format(as.POSIXct(t[o], origin = "1970-01-01", tz = "UTC"),
                           "%Y-%m-%d"),
             loc = loc[o], prov = prov[o],
             dx1 = dx1[o], dx2 = dx2[o], dx3 = dx3[o],
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Validation -- the invariants the teaching data has to satisfy
# ---------------------------------------------------------------------------

left_filled_problem <- function(M) {
  filled <- M != ""
  n_filled <- rowSums(filled)
  last_filled <- max.col(filled * col(filled), ties.method = "last")
  last_filled[n_filled == 0] <- 0
  any(n_filled != last_filled)
}

validate_dad <- function(dad, hf_prefixes) {
  problems <- character(0)
  M <- as.matrix(dad[, paste0("DXCODE", seq_len(MAX_DX)), drop = FALSE])

  if (any(M[, 1] == "")) problems <- c(problems, "DXCODE1 empty on some rows")
  if (anyDuplicated(dad$ENCOUNTER_ID)) problems <- c(problems, "duplicate ENCOUNTER_ID")
  if (left_filled_problem(M)) problems <- c(problems, "gap in DXCODE sequence")
  dup <- apply(M, 1, function(r) { r <- r[r != ""]; anyDuplicated(r) > 0 })
  if (any(dup)) problems <- c(problems, "duplicate diagnosis code within an abstract")
  if (any(dad$DISCH_DTTM <= dad$ADMIT_DTTM))
    problems <- c(problems, "discharge not after admit")

  obs <- mean(dad$DXCODE1 %in% hf_prefixes)
  if (obs < 0.10)
    problems <- c(problems, sprintf(
      "heart-failure MRDx rate %.1f%% is below the 10%% floor", 100 * obs))
  list(hf_rate = obs, problems = problems)
}

validate_claims <- function(claims, patient_ids, start, end) {
  problems <- character(0)
  M <- as.matrix(claims[, paste0("DXCODE", seq_len(MAX_CLAIM_DX)), drop = FALSE])
  if (any(M[, 1] == "")) problems <- c(problems, "claims DXCODE1 empty on some rows")
  if (left_filled_problem(M)) problems <- c(problems, "gap in claims DXCODE sequence")
  dup <- apply(M, 1, function(r) { r <- r[r != ""]; anyDuplicated(r) > 0 })
  if (any(dup)) problems <- c(problems, "duplicate diagnosis code within a claim")
  if (!all(claims$PATIENT_ID %in% patient_ids))
    problems <- c(problems, "claim for a PATIENT_ID that is not in the DAD file")
  if (any(claims$SERVICE_DATE < start) || any(claims$SERVICE_DATE > end))
    problems <- c(problems, "claim service date outside the study window")
  if (anyDuplicated(claims$CLAIM_ID)) problems <- c(problems, "duplicate CLAIM_ID")
  problems
}

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

usage <- function() paste0(
  "Usage: Rscript generate_icd10ca_data.R [options]\n\n",
  "  --patients N          number of patients to generate (default 20000)\n",
  "  --dad-out PATH        DAD output CSV (default dad_synthetic.csv)\n",
  "  --claims-out PATH     claims output CSV (default claims_synthetic.csv)\n",
  "  --truth               also write patient-level latent comorbidity flags\n",
  "  --truth-out PATH      where to write them (default patient_truth.csv)\n",
  "  --report              print the per-condition capture table\n",
  "  --seed N              RNG seed (default 42)\n",
  "  --hf-rate P           target share of abstracts whose DXCODE1 is a\n",
  "                        heart-failure code; must be >= 0.10 (default 0.15)\n",
  "  --hf-patient-rate P   share of patients who have heart failure (0.14)\n",
  "  --start YYYY-MM-DD    earliest admission date (default 2022-01-01)\n",
  "  --end YYYY-MM-DD      latest admission date; follow-up is censored\n",
  "                        here (default 2024-12-31)\n",
  "  --base-gap-days D     mean days to the next admission for a reference\n",
  "                        patient: age 70, no comorbidity, no HF (1000)\n",
  "  --transfer-rate P     share of stays ending in an acute-care transfer\n",
  "  --decimal             write ICD-10-CA codes with the decimal (I50.0).\n",
  "                        Default is the DAD convention without it (I500).\n")

parse_args <- function(argv) {
  opts <- list(patients = 20000L, dad_out = "../data/dad_synthetic.csv",
               claims_out = "../data/claims_synthetic.csv",
               truth_out = "patient_truth.csv", seed = 42L,
               hf_rate = 0.15, hf_patient_rate = 0.14,
               start = "2022-01-01", end = "2024-12-31",
               base_gap_days = 1000, transfer_rate = 0.045,
               decimal = FALSE, truth = FALSE, report = FALSE)
  keys <- c("--patients" = "patients", "--dad-out" = "dad_out",
            "--claims-out" = "claims_out", "--truth-out" = "truth_out",
            "--seed" = "seed", "--hf-rate" = "hf_rate",
            "--hf-patient-rate" = "hf_patient_rate", "--start" = "start",
            "--end" = "end", "--base-gap-days" = "base_gap_days",
            "--transfer-rate" = "transfer_rate")
  argv <- unlist(lapply(argv, function(a)
    if (grepl("^--[a-z-]+=", a)) c(sub("=.*$", "", a), sub("^[^=]*=", "", a)) else a))
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("--help", "-h")) { cat(usage()); quit(status = 0) }
    if (a == "--decimal") { opts$decimal <- TRUE; i <- i + 1L; next }
    if (a == "--truth")   { opts$truth   <- TRUE; i <- i + 1L; next }
    if (a == "--report")  { opts$report  <- TRUE; i <- i + 1L; next }
    if (is.na(keys[a])) stop("unknown argument: ", a)
    if (i == length(argv)) stop("missing value for ", a)
    key <- unname(keys[a]); val <- argv[i + 1L]
    opts[[key]] <- if (key %in% c("patients", "seed")) as.integer(val)
                   else if (key %in% c("hf_rate", "hf_patient_rate",
                                       "base_gap_days", "transfer_rate")) as.numeric(val)
                   else val
    i <- i + 2L
  }
  opts
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_args(argv)
  if (opts$hf_rate < 0.10) stop("--hf-rate must be at least 0.10")
  if (opts$patients < 1) stop("--patients must be positive")

  start_dt <- as.POSIXct(paste(opts$start, "00:00:00"), tz = "UTC")
  end_dt   <- as.POSIXct(paste(opts$end, "23:59:00"), tz = "UTC")
  if (end_dt <= start_dt) stop("--end must be after --start")

  set.seed(opts$seed)
  fmt <- if (opts$decimal) identity else function(x) gsub(".", "", x, fixed = TRUE)
  n_pt <- opts$patients

  # -- phase 1 --------------------------------------------------------------
  pts <- make_patients(n_pt, opts$hf_patient_rate, opts$start)

  # -- phase 2 --------------------------------------------------------------
  enc_by_pt <- vector("list", n_pt)
  death_dt <- rep(NA_real_, n_pt)
  for (i in seq_len(n_pt)) {
    res <- gen_patient_encounters(i, pts, start_dt, end_dt,
                                  opts$base_gap_days, opts$transfer_rate)
    enc_by_pt[[i]] <- res$enc
    if (!is.na(res$death)) death_dt[i] <- as.numeric(res$death)
  }
  n_enc_by_pt <- vapply(enc_by_pt, length, integer(1))
  n_enc <- sum(n_enc_by_pt)

  # Decide which abstracts get heart failure as the most responsible
  # diagnosis, so the encounter-level rate lands on --hf-rate.
  eligible_hf_enc <- sum(vapply(seq_len(n_pt), function(i) {
    if (!pts$is_hf[i]) return(0L)
    sum(vapply(enc_by_pt[[i]], function(e) as.integer(e$age >= 40), integer(1)))
  }, integer(1)))
  p_hf_mrdx <- if (eligible_hf_enc > 0) opts$hf_rate * n_enc / eligible_hf_enc else 0
  hf_capped <- p_hf_mrdx > 0.92 || p_hf_mrdx < 0.05
  p_hf_mrdx <- min(max(p_hf_mrdx, 0.05), 0.92)

  # -- phase 3 --------------------------------------------------------------
  dad_cols <- c("PATIENT_ID", "ENCOUNTER_ID", "AGE", "SEX", "ADMIT_DTTM",
                "DISCH_DTTM", "ADMIT_CATEGORY", "INST_ID", "DISCH_DISP",
                paste0("DXCODE", seq_len(MAX_DX)))
  dad <- matrix("", nrow = n_enc, ncol = length(dad_cols),
                dimnames = list(NULL, dad_cols))
  facilities <- sprintf("F%03d", 1:12)
  sens_vec <- vapply(CHRONIC, function(cc) cc$dad_sens, numeric(1))
  fp_vec   <- vapply(CHRONIC, function(cc) cc$dad_fp, numeric(1))

  stay_admit <- numeric(n_enc); stay_disch <- numeric(n_enc)
  stay_pt <- integer(n_enc)
  row_i <- 0L

  for (i in seq_len(n_pt)) {
    kvec <- pts$known[i, ]
    p_code <- ifelse(kvec, sens_vec, fp_vec)
    home_fac <- sample(facilities, 1)
    for (e in enc_by_pt[[i]]) {
      row_i <- row_i + 1L
      hf_mrdx <- pts$is_hf[i] && e$age >= 40 && runif(1) < p_hf_mrdx
      codes <- build_diagnoses(e$age, pts$sex[i], hf_mrdx, kvec, p_code,
                               e$n_dx_hint, pts$is_hf[i])
      fac <- if (e$is_transfer_in) sample(setdiff(facilities, home_fac), 1) else home_fac

      dad[row_i, 1:9] <- c(pts$patient_id[i], sprintf("E%07d", row_i),
                           as.character(e$age), pts$sex[i],
                           format(e$admit, "%Y-%m-%d %H:%M"),
                           format(e$disch, "%Y-%m-%d %H:%M"),
                           e$admit_cat, fac, e$disp)
      dad[row_i, 9 + seq_along(codes)] <- fmt(codes)
      stay_admit[row_i] <- as.numeric(e$admit)
      stay_disch[row_i] <- as.numeric(e$disch)
      stay_pt[row_i] <- i
    }
  }
  dad <- as.data.frame(dad, stringsAsFactors = FALSE)

  # -- phase 4 --------------------------------------------------------------
  birth_num <- as.numeric(as.POSIXct(paste(pts$birth, "00:00:00"), tz = "UTC"))
  win_lo <- pmax(as.numeric(start_dt), birth_num)
  win_hi <- ifelse(is.na(death_dt), as.numeric(end_dt),
                   pmin(as.numeric(end_dt), death_dt))
  years <- (win_hi - win_lo) / SECS_PER_DAY / 365.25
  cl <- gen_claims(pts, win_lo, win_hi, years, stay_admit, stay_disch, stay_pt)

  claims <- data.frame(
    CLAIM_ID = sprintf("C%08d", seq_len(nrow(cl))),
    PATIENT_ID = pts$patient_id[cl$pt],
    SERVICE_DATE = cl$date,
    SERVICE_LOCATION = cl$loc,
    PROVIDER_TYPE = cl$prov,
    DXCODE1 = cl$dx1, DXCODE2 = cl$dx2, DXCODE3 = cl$dx3,
    stringsAsFactors = FALSE)

  # -- write ----------------------------------------------------------------
  write.csv(dad, opts$dad_out, row.names = FALSE, quote = FALSE, na = "")
  write.csv(claims, opts$claims_out, row.names = FALSE, quote = FALSE, na = "")

  if (opts$truth) {
    # TRUE_* is the underlying condition, which is what drives readmission
    # risk here. DX_* is the subset that has actually been diagnosed, which
    # is the fairest reference standard for a coding algorithm -- a chart
    # review would not find what no clinician has recorded either.
    truth <- data.frame(PATIENT_ID = pts$patient_id,
                        AGE_AT_STUDY_START = pts$age0,
                        SEX = pts$sex,
                        TRUE_HF = as.integer(pts$is_hf),
                        stringsAsFactors = FALSE)
    for (nm in CONDITION_NAMES) {
      truth[[paste0("TRUE_", nm)]] <- as.integer(pts$conds[, nm])
      truth[[paste0("DX_", nm)]] <- as.integer(pts$known[, nm])
    }
    truth$DIED_IN_HOSPITAL <- as.integer(!is.na(death_dt))
    write.csv(truth, opts$truth_out, row.names = FALSE, quote = FALSE, na = "")
  }

  # -- validate and report --------------------------------------------------
  hf_prefixes <- fmt(c(HF_CODES, HF_HTN_CODE, HF_HTN_CKD_CODE))
  v <- validate_dad(dad, hf_prefixes)
  problems <- c(v$problems,
                validate_claims(claims, pts$patient_id, opts$start, opts$end))
  if (hf_capped)
    problems <- c(problems,
                  "could not hit --hf-rate exactly; adjust --hf-patient-rate")

  dadM <- as.matrix(dad[, paste0("DXCODE", seq_len(MAX_DX))])
  dx_counts <- rowSums(dadM != "")

  # Readmission sanity check, using the naive definition (any later
  # admission, transfers included) and a simple episode-of-care definition
  # (ignore an admission starting within 12 hours of the last discharge).
  ord <- order(dad$PATIENT_ID, dad$ADMIT_DTTM)
  d <- dad[ord, ]
  same_pt <- c(FALSE, d$PATIENT_ID[-1] == d$PATIENT_ID[-nrow(d)])
  gap_h <- c(NA, as.numeric(difftime(as.POSIXct(d$ADMIT_DTTM[-1], tz = "UTC"),
                                     as.POSIXct(d$DISCH_DTTM[-nrow(d)], tz = "UTC"),
                                     units = "hours")))
  naive_30 <- same_pt & !is.na(gap_h) & gap_h <= 24 * 30
  episode_30 <- naive_30 & gap_h > 12

  # How much of the truth each source recovers, condition by condition.
  pt_of_row <- match(dad$PATIENT_ID, pts$patient_id)
  cl_codes <- cbind(claims$DXCODE1, claims$DXCODE2, claims$DXCODE3)
  cl_pt <- match(claims$PATIENT_ID, pts$patient_id)
  capture <- function(nm) {
    cc <- CHRONIC[[nm]]
    dad_hit <- matrix(dadM %in% fmt(cc$icd10), nrow = nrow(dadM))
    pt_dad <- logical(n_pt)
    pt_dad[unique(pt_of_row[rowSums(dad_hit) > 0])] <- TRUE
    hit <- rowSums(matrix(cl_codes %in% cc$icd9, nrow = nrow(cl_codes))) > 0
    tab <- table(cl_pt[hit])
    pt_cl <- logical(n_pt)
    pt_cl[as.integer(names(tab))[tab >= 2]] <- TRUE
    ref <- pts$known[, nm]
    c(true = mean(pts$conds[, nm]), diagnosed = mean(ref),
      dad = mean(pt_dad), claims = mean(pt_cl),
      sens_dad = if (any(ref)) mean(pt_dad[ref]) else NA_real_,
      sens_claims = if (any(ref)) mean(pt_cl[ref]) else NA_real_,
      spec_dad = if (any(!ref)) mean(!pt_dad[!ref]) else NA_real_,
      spec_claims = if (any(!ref)) mean(!pt_cl[!ref]) else NA_real_)
  }
  cap <- t(vapply(CONDITION_NAMES, capture, numeric(8)))

  cat(sprintf("Wrote %d encounters for %d patients to %s\n",
              nrow(dad), n_pt, opts$dad_out))
  cat(sprintf("Wrote %d claims to %s\n", nrow(claims), opts$claims_out))
  if (opts$truth) cat(sprintf("Wrote patient-level truth to %s\n", opts$truth_out))
  cat(sprintf("  heart failure as DXCODE1 : %.1f%% of abstracts\n", 100 * v$hf_rate))
  cat(sprintf("  encounters per patient   : min %d, median %d, max %d\n",
              min(n_enc_by_pt), as.integer(median(n_enc_by_pt)), max(n_enc_by_pt)))
  cat(sprintf("  diagnoses per abstract   : min %d, median %d, max %d\n",
              min(dx_counts), as.integer(median(dx_counts)), max(dx_counts)))
  cat(sprintf("  claims per patient-year  : median %.1f\n",
              median((table(factor(claims$PATIENT_ID, levels = pts$patient_id)) /
                      pmax(years, 0.01))[years > 0.5])))
  cat(sprintf("  in-hospital deaths       : %.1f%% of abstracts\n",
              100 * mean(dad$DISCH_DISP == "07")))
  cat(sprintf("  acute-care transfers out : %.1f%% of abstracts\n",
              100 * mean(dad$DISCH_DISP == "01")))
  cat(sprintf("  30-day readmission       : %.1f%% naive, %.1f%% after a 12h\n",
              100 * mean(naive_30), 100 * mean(episode_30)))
  cat("                             episode-of-care rule\n")
  cat(sprintf("  conditions modelled      : %d (Charlson %d, Elixhauser %d, Tonelli %d)\n",
              length(CONDITION_NAMES),
              sum(vapply(CHRONIC, function(c) "charlson" %in% c$algos, logical(1))),
              sum(vapply(CHRONIC, function(c) "elixhauser" %in% c$algos, logical(1))),
              sum(vapply(CHRONIC, function(c) "tonelli" %in% c$algos, logical(1)))))
  cat(sprintf("  hypertension             : %.1f%% underlying, %.1f%% diagnosed\n",
              100 * cap["HTN", "true"], 100 * cap["HTN", "diagnosed"]))
  cat(sprintf("    from the DAD           : sens %.2f, spec %.2f\n",
              cap["HTN", "sens_dad"], cap["HTN", "spec_dad"]))
  cat(sprintf("    2+ claims in the window: sens %.2f, spec %.2f\n",
              cap["HTN", "sens_claims"], cap["HTN", "spec_claims"]))
  cat(sprintf("  ICD-10-CA code format    : %s\n",
              if (opts$decimal) "with decimal (I50.0)" else "DAD style (I500)"))

  if (opts$report) {
    cat("\nPer-condition capture (patient level; claims rule = 2+ claims,\n",
        "reference = diagnosed condition; algorithms: C Charlson, E Elixhauser,\n",
        " T Tonelli, all Quan 2005 / Tonelli 2015):\n\n", sep = "")
    cat(sprintf("%-34s %-5s %6s %6s %6s %6s %6s\n", "condition", "algo",
                "true", "dx", "DAD", "claims", "sensD"))
    for (nm in CONDITION_NAMES) {
      a <- CHRONIC[[nm]]$algos
      tag <- paste0(if ("charlson" %in% a) "C" else "-",
                    if ("elixhauser" %in% a) "E" else "-",
                    if ("tonelli" %in% a) "T" else "-")
      cat(sprintf("%-34s %-5s %5.1f%% %5.1f%% %5.1f%% %5.1f%% %6.2f\n",
                  substr(CHRONIC[[nm]]$label, 1, 34), tag,
                  100 * cap[nm, "true"], 100 * cap[nm, "diagnosed"],
                  100 * cap[nm, "dad"], 100 * cap[nm, "claims"],
                  cap[nm, "sens_dad"]))
    }
  }

  if (length(problems)) {
    cat("\nVALIDATION FAILURES:\n", file = stderr())
    for (m in head(problems, 20)) cat("  ", m, "\n", sep = "", file = stderr())
    return(1L)
  }
  cat("  validation               : all checks passed\n")
  0L
}

if (!interactive()) quit(status = main())
