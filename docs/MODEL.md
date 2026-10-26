# Model specification and parameter provenance

This document describes what `hbocrisk` implements, where each parameter comes
from, and where the implementation had to make a choice that the publications do not
pin down. References are listed at the end; tags such as `[LEE2022]` are used in the
source code as well.

## 1. Genetic model

**Major genes.** Pathogenic variants (PVs) are alleles of one "major gene" locus
(wild type plus one allele per gene), in dominance order [LEE2016, LEE2022]:

| model | genes (dominance order) |
|---|---|
| BOADICEA (breast) | BRCA1, BRCA2, PALB2, CHEK2, ATM, BARD1, RAD51C, RAD51D |
| EOC (ovarian) | BRCA1, BRCA2, RAD51D, RAD51C, BRIP1, PALB2 |

A genotype is an unordered allele pair (45 genotypes for 8 genes). Its *risk class* is
its highest-ranked PV. Founders are in Hardy–Weinberg equilibrium. With Ashkenazi
ancestry, the founder frequencies are the Ashkenazi ones, while the population
incidences are partitioned with the UK frequencies (CanRisk `PEDIGREE_ALLELE_FRQ`).

**Polygene.** Hypergeometric polygenic model with 2N = 6 alleles
(7 polygenotypes, x = (R−3)/√1.5) [CAN1978, LAN1997, LEE2019].

**Incidence.** For cancer c, risk class μ, polygenotype x and risk-factor category r:

    λ(t) = h_μ(t) · exp(β_PG(t)·x + g_r(t))

with β_PG(t) the residual polygenic standard deviation. The class baselines h_μ(t) are
solved age by age, so that among individuals still unaffected:

* the mean incidence of carriers of gene μ equals RR_μ(t) × population incidence,
  where the RRs (relative to the population, averaged over all other effects) are those in
  [LEE2022] Table 1, with BRCA2 breast RRs ≥ 60 years from [FIC2024];
* the mean incidence over all classes equals the population incidence
  [ANT2001, LEE2016, LEE2019].

Breast cancer (female and male) has a polygenic component in BOADICEA. Ovarian,
prostate and pancreatic cancers depend only on the major genotype. In the EOC model only
ovarian cancer is polygenic.

**Polygenic variances.**

* Female breast: σ²(t) = 4.83 − 0.05961·t [ANT2008 Table 6, LEE2022 supplement]. The
  variance explained by the explicitly modelled genes is subtracted:
  Var_μ(t) = log[((1−q)² + q(2−q)RR²) / ((1−q)² + q(2−q)RR)²] [LEE2022 suppl.]
  for PALB2, CHEK2, ATM, BARD1, RAD51C and RAD51D. For the CanRisk v2.4 BRCA2 update
  the BRCA2 term of the earlier RRs is added back and the new one subtracted, which keeps
  the total familial variance constant [FIC2024].
* Male breast: 1.4² = 1.96 [ANT2008, LEE2022] minus the PALB2 term.
* EOC model: 1.434 [LEE2022O] minus the RAD51D, RAD51C, BRIP1 and PALB2 terms
  (age-specific; the paper quotes 1.4156 before PALB2 was added).

**Other cancers.** Male breast, prostate and pancreatic RRs: [ANT2008] Table 7 and
[LEE2022] Table 1. Contralateral breast cancer incidence is half the genotype-specific
breast incidence (one breast at risk) [ANT2004, ANT2008]. Cancers are independent given
genotype (and polygene / risk factors).

**Tumour pathology.** For a first female breast cancer, the probability of the
observed subtype given the class is P(ER−), P(TN | ER−) and the CK5/6/CK14 distribution
among TN tumours ([LEE2022] tables s3–s4, [MAV2010] Table 2). The non-carrier proportions are
solved so that the population mixture matches the general-population distribution, as
in [MAV2010]. Pathology codes follow CanRisk: ER status first; PR and HER2 are used
jointly (TN or not); CK14 and CK5/6 are used only if both are known and the tumour is TN.

**Genetic tests.** 100% specificity. A negative mutation search for gene g has
likelihood 1 − sensitivity_g in genotypes carrying g; a negative direct test excludes
them. A positive test fixes the risk class to the highest-ranked positive gene (CanRisk
FAQ "Mutation carrier probabilities in families with … two or more mutations"). Default
sensitivities are those of CanRisk ≥ 4.1; the legacy values are in
`params.TEST_SENSITIVITY_LEGACY`.

## 2. Risk factors, mammographic density, PRS

Only the proband's risk factors are used. Factors that are not observed are left out
of the model rather than marginalised [LEE2019].

* BOADICEA: menarche, parity, age at first live birth, oral contraception, MHT, BMI,
  alcohol, age at menopause, BI-RADS density ([LEE2019] tables s2.1–s2.10), and height
  as a continuous factor ([LEE2022]: five bins with B(4, ½) masses plus a zero-mass
  category at the measured value; log-RR 0.1013 per SD; mean 162.81 cm, SD 6.452 cm).
* EOC model: parity, OC use (years), MHT, tubal ligation, endometriosis, BMI ([LEE2022O]
  table s2), and height (continuous).
* Age-dependent RRs (split ages, "from current age" for OC/MHT) are implemented as in
  the tables. For age at first live birth the published distribution sums to 0.72; the
  remaining 0.28 (nulliparous) is a reference category.
* The joint risk-factor distribution is the product of the marginals. When there are more than
  1,500 joint categories (default), adjacent categories, ordered by mean log-RR, are
  merged into mass-weighted bins. The proband's own category is never merged
  (cf. [LEE2019] suppl. §2). Tested error: under 0.05% relative.
* PRS: x_P = α·z + x_R. The proband's penetrance is multiplied by
  f(x_PRS | x_P ∈ bin k) ([LEE2019] eq. s1.14) with bins at the normal quantiles of the
  binomial masses (s1.13).
* PRS ancestry [FIC2025E]: the raw score is ancestry-independent; its mean μ_i, SD σ_i and
  α_i depend on the genetic ancestry i (European, African, East Asian, South Asian).
  `BCAC 309` (BCAC 313 without 4 SNPs not imputed in UK Biobank), `BCAC 307` (also
  without the 2 SNPs correlated with CHEK2 c.1100del) and `OCAC 36` carry these
  parameters and ancestry-specific allele frequencies q_ik (CanRisk files
  `*-UKB_all.prs`, SHARE-PRScalculation v4.2.12, GPL-3.0; downloaded at runtime, not
  bundled). If the main ancestry proportion P_i
  reaches the threshold (0.8; 0.7 for South Asian), z = (R − μ_i)/σ_i and α = α_i.
  Otherwise (mixed genetic ancestry)

      σ² = Σ_i P_i σ_i² + 4 Σ_k w_k² (Σ_i P_i q_ik² − (Σ_i P_i q_ik)²)
      z  = Σ_i P_i α_i (R − μ_i)/σ_i / (σ · Σ_i P_i α_i/σ_i),   α = Σ_i P_i α_i σ_i / σ

  (`PRSModel.normalise`; agrees with vcf2prs to 1e-7). Missing variants are imputed with
  the ancestry-specific (proportion-weighted) frequencies. Unknown ancestry uses the
  population parameters; the European-only reference PRS (e.g. BCAC 313) accept only
  European or unknown ancestry.

## 3. Pedigree likelihood and risks

The pedigree likelihood is computed by sum-product variable elimination over each
person's genotype and polygenotype variables. Transmission is factorised (Mendelian for the major
locus, hypergeometric for the polygene). A greedy minimum-factor-size elimination order
is used. This is equivalent to Elston–Stewart peeling for loop-free pedigrees and also
handles loops. MZ twins share one set of genetic variables. Individuals with no data and
no descendants are pruned.

Risks follow the reformulation of [LEE2022] (suppl. eq. s.2): the proband's posterior
over (genotype, polygenotype), given all pedigree information including her own current
status, is combined with her genotype-specific future survival:

    risk(a → T) = Σ_{g,x} P(g, x | data) · (1 − exp(−Σ_{u=a}^{T−1} λ_{g,x}(u)))

Time is discrete in years. Survival to age a is exp(−Σ_{u<a} λ(u)), and the probability of
diagnosis at age d is S(d)(1 − e^{−λ(d)}). Follow-up is censored at 80, and diagnoses at 80 or
later are ignored (CanRisk FAQ). An age of "AU" (affected, age unknown) is set to the age at
last follow-up. Individuals without both age and year of birth contribute only
genetic-test information.

CanRisk outputs reproduced: remaining lifetime risks (ages a+1…a+5, a+10, multiples of
5 up to 80), lifetime risk (proband age set to 20, to age 80), 10-year risk (age 40 to
50), the matching population baselines (proband only: sex, age, year of birth, first
breast cancer age), and carrier probabilities by risk class. Competing mortality is not
modelled, as in BOADICEA.

## 4. Population incidences

Birth-cohort-specific incidences are derived as in [LEE2014, LEE2022]. For each decadal
cohort (1910–19 … 1980–89; earlier and later births use the first and last cohort), each
single year of age is averaged over the birth years in the cohort, using the calendar-year
incidence for the 5-year age group. Years outside the data use the mean of the first or
last five years. The result is smoothed with LOWESS (bandwidth 0.2).

Source: IARC CI5plus (detailed, 2024 v1.1). For the UK (England), ONS registration counts for
2009–2017 replace the CI5plus counts, because CanRisk uses ONS data from 1992 onwards and
CI5plus has 2–4% fewer breast and 15–20% fewer ovarian registrations. Ovarian cancer is
C56 + C57 (ovary and fallopian tube). Germany (Hamburg, Bremen,
Schleswig-Holstein, Saarland), Austria and Switzerland are included as extras.
Finland and Sweden are not in CI5plus and are not included.

**UK ethnic groups (BOADICEA BC v7 / EOC v3)** [FIC2025E]. With `Settings.ethnicity`
(or the CanRisk header `##ethnicity=<ONS group>;<background>`), the UK incidences are
multiplied by the ethnic-group incidence rate ratios relative to the UK population
(Table 1; South and East Asian grouped; male breast cancer unchanged; 'Other' unchanged).
Non-White groups use the non-European CHEK2 frequency (0.00109) when the mutation
frequency set is the default `UK`. Baseline (population) risks use the same ethnic
group. As in CanRisk, an ethnicity can only be combined with the UK incidences.
`register_population()` lets you supply your own tables.

## 5. Known differences from CanRisk / not implemented

* The exact CanRisk incidence tables are not published. UK breast cancer risks are
  within about 0.5–1% (relative) of CanRisk 7.3.2 outputs. Ovarian cancer risks are about
  3–5% lower. Other populations have not been compared with CanRisk.
* BOADICEA v7 ethnicity: the incidence rate ratios, CHEK2 frequencies and the
  ancestry-specific PRS are implemented. The ethnicity-specific risk-factor and BI-RADS
  distributions (Tables S18, S21) and breast tumour subtype distributions (Tables
  S4–S6) are not; the default distributions are used and a warning is given when risk
  factors, BI-RADS or pathology are supplied. For unspecified ethnicity v7 equals the
  v6 parameterisation used here. Ovarian cancer risks of BRCA2 carriers by ethnicity
  agree with [FIC2025E] (to age 80: White 15.2/15, Black 10.3/10, Asian 11.2/11, Mixed
  8.1/8 %, hbocrisk/published). Breast cancer risks of non-White BRCA2 carriers are 3–4
  percentage points lower than published (White 58.8/58, Black 51.3/54, South Asian
  47.3/50, Mixed 42.9/47 %, born 2000); the cause is not known.
* Continuous mammographic density (Volpara/STRATUS, v7.2 [FIC2025M]) is not implemented:
  its regression/normalisation parameters are only in a data supplement that was not
  available. The v7.2 age-specific BI-RADS population distributions are not implemented
  either; the v5 distribution is used.
* The prostate cancer model (HOXB13) is not implemented.
* MENDEL's exact discrete-time conventions are not documented; the yearly conventions above
  are used.

## References

* [ANT2001] Antoniou AC et al. Genet Epidemiol 2001;21:1–18.
* [ANT2004] Antoniou AC et al. Br J Cancer 2004;91:1580–90.
* [ANT2008] Antoniou AC et al. Br J Cancer 2008;98:1457–66.
* [CAN1978] Cannings C, Thompson EA, Skolnick MH. Adv Appl Prob 1978;10:26–61.
* [LAN1997] Lange K. Genetics 1997;147:1423–30.
* [MAV2010] Mavaddat N et al. Breast Cancer Res 2010;12:R28.
* [LEE2014] Lee AJ et al. Br J Cancer 2014;110:535–45.
* [LEE2016] Lee AJ et al. Genet Med 2016;18:1190–8.
* [LEE2019] Lee A et al. Genet Med 2019;21:1708–18 (and supplementary material).
* [LEE2022] Lee A et al. J Med Genet 2022;59:1206–18 (BOADICEA v6, suppl.).
* [LEE2022O] Lee A et al. J Med Genet 2022;59 (EOC model; doi:10.1136/jmedgenet-2021-107904, suppl.).
* [MAV2023] Mavaddat N et al. Cancer Epidemiol Biomarkers Prev 2023 (alternative PRS in BOADICEA).
* [FIC2024] Ficorella L et al. BJC Reports 2024;2:53.
* [FIC2025E] Ficorella L et al. Br J Cancer 2025 (ethnically diverse UK population).
* [FIC2025M] Ficorella L et al. JCO Precis Oncol 2025;9:e2500203.
* CanRisk documentation (canrisk.atlassian.net) and the CanRisk web-service code
  (github.com/CCGE-BOADICEA/bws) for file formats, allele frequencies, default
  sensitivities and output conventions.
