# hbocrisk

**H**ereditary **B**reast and **O**varian **C**ancer **risk**: an open Python
re-implementation of the **BOADICEA** breast cancer risk model (v6 with the CanRisk v2.4 BRCA2 update; v7-compatible for unspecified ethnicity) and the
**CanRisk epithelial tubo-ovarian cancer (EOC) model** (v2). It is built from the
published model descriptions and supplementary material. The only dependencies are
`numpy` and `scipy`.

It computes, for the target individual of a pedigree:

* probabilities of carrying a pathogenic variant (BRCA1, BRCA2, PALB2, CHEK2, ATM, BARD1,
  RAD51C, RAD51D; EOC model: BRCA1, BRCA2, RAD51D, RAD51C, BRIP1, PALB2)
* remaining-lifetime, lifetime (20–80) and 10-year (40–50) breast or ovarian cancer risks,
  contralateral breast cancer risk for affected women, and the matching population
  baselines

It uses the following inputs:

* explicit family history of breast (first and contralateral), ovarian, prostate and
  pancreatic cancer, for pedigrees of any structure (including loops and MZ twins)
* genetic test results (mutation search or direct test) and breast tumour pathology
  (ER, PR, HER2, CK14, CK5/6)
* questionnaire-based risk factors, BI-RADS mammographic density and height (continuous)
* polygenic risk scores, given as z-score and α or computed from a VCF with the CanRisk
  reference PRS (BCAC 313 etc.)
* birth cohort, country-specific cancer incidences, and Ashkenazi Jewish ancestry

The model structure, every parameter with its source, and the known deviations from
CanRisk are documented in [docs/MODEL.md](docs/MODEL.md).

> **Not a medical device.** This is an independent research implementation. It has not
> been validated for clinical use. The reference implementation is the CanRisk tool
> (www.canrisk.org) from the University of Cambridge. BOADICEA is licensed by Cambridge
> Enterprise. Check the licensing situation before any non-research use.

## Installation

```bash
pip install git+https://github.com/christopher-schroeder/hbocrisk.git

# development install from a clone
git clone git@github.com:christopher-schroeder/hbocrisk.git
cd hbocrisk
pip install -e .            # from this directory
pip install -e .[test]      # with pytest
```

## Usage

### Command line

```bash
hbocrisk risk family.canrisk                 # BOADICEA (breast)
hbocrisk risk family.canrisk --model both    # breast and ovarian models
hbocrisk risk family.canrisk --population Germany --json
hbocrisk risk family.canrisk --prs-vcf sample.vcf --prs-name "BCAC 313" --prs-sample S1
hbocrisk populations                         # available incidence tables
hbocrisk serve --port 8080                   # JSON HTTP service (see below)
```

Input files use the CanRisk pedigree format (v1–v4) or the BOADICEA v4 import format.
Risk factors and PRS go in the file header (`##menarche=13`, `##BMI=24`,
`##birads=c`, `##PRS_BC=alpha=0.501,zscore=1.2`, ...).

### Python

```python
from hbocrisk import Pedigree, Person, GeneticTest, calculate, Settings

ped = Pedigree("FAM1", [
    Person("m", "F", age=52, yob=1960, bc1=45, pathology={"er": "N", "pr": "N", "her2": "N"}),
    Person("f", "M", age=60, yob=1958),
    Person("p", "F", father="f", mother="m", target=True, age=35, yob=1990,
           tests={"BRCA1": GeneticTest("S", "N"), "BRCA2": GeneticTest("S", "N")}),
], risk_factors={"parity": 1, "bmi": 23.5, "birads": "c", "height": 168},
   prs={"BC": (0.501, 1.1)})                         # (alpha, z-score)

res = calculate(ped, Settings(model="BC", population="UK"))
res.carrier_probabilities           # {'no_pathogenic_variant': ..., 'BRCA1': ..., ...}
res.remaining_lifetime["breast"]    # [(36, risk), ..., (80, risk)]
res.lifetime, res.ten_year, res.baseline_remaining_lifetime
res.to_dict()                       # JSON-serialisable
```

`Pedigree.from_dict` / `to_dict`, `read_canrisk` / `write_canrisk` convert between
formats. `hbocrisk.prs` computes PRS z-scores from VCFs or dosage dictionaries.
`register_population()` registers your own incidence tables, for example national
registry data.

### JSON service (backend for a browser interface)

`hbocrisk serve` starts a standard-library HTTP server:

* `POST /calculate` with `{"canrisk": "<file text>"}` or `{"pedigree": {...}}`, plus an
  optional `"settings": {"model": "BC"|"OC", "population": "UK",
  "mutation_frequencies": "UK", "sensitivities": {"BRCA1": 0.98}}`
* `GET /populations`

`hbocrisk.server.handle_request(payload)` gives the same logic for any web
framework.

## Incidence data

Cohort-specific incidences are built with `scripts/build_incidences.py` from IARC
CI5plus. For the UK they also use ONS England registrations 2009–2017, parsed with
`scripts/parse_ons_england.py`. Populations: UK, Australia, Canada, USA (SEER whites),
Denmark, Estonia, France, Iceland, Netherlands, New Zealand, Norway, Slovenia, Spain and,
not offered in CanRisk, Germany (Hamburg, Bremen, Schleswig-Holstein and Saarland
registries), Austria and Switzerland. Data sources: IARC CI5plus (Ferlay et al.) and ONS
(Open Government Licence).

## Validation

Reproduce with `python scripts/validate.py`. Reference values come from CanRisk
web-service outputs (`boadicea model 7.3.2`, `ovarian model 3.2.2`) and from the
published BOADICEA papers. The Lee et al. 2022 values come from v6 before the 2024 BRCA2
update; that update increases family-history effects slightly, which explains the
small positive differences in the "mother affected" rows.

| scenario | hbocrisk (%) | published (%) | rel. diff | source |
|---|---|---|---|---|
| population, age 25: lifetime BC (20-80) | 11.94 | 12 | -0.5% | CanRisk 7.3.2 web-service test output |
| population, age 25: 10-year BC (40-50) | 1.708 | 1.718 | -0.6% | CanRisk 7.3.2 web-service test output |
| population, age 25: BC risk to 80 | 11.92 | 11.99 | -0.5% | CanRisk 7.3.2 web-service test output |
| population, age 25: P(BRCA1 carrier) | 0.1272 | 0.1271 | +0.0% | CanRisk 7.3.2 web-service test output |
| population, age 25: P(BRCA2 carrier) | 0.2031 | 0.2031 | +0.0% | CanRisk 7.3.2 web-service test output |
| population, age 25: P(PALB2 carrier) | 0.1276 | 0.1276 | +0.0% | CanRisk 7.3.2 web-service test output |
| population, age 25: P(CHEK2 carrier) | 0.7426 | 0.7426 | +0.0% | CanRisk 7.3.2 web-service test output |
| population, age 25: P(ATM carrier) | 0.3574 | 0.3574 | +0.0% | CanRisk 7.3.2 web-service test output |
| population, age 64 (yob 1962): lifetime BC | 11.82 | 11.9 | -0.6% | CanRisk 7.3.2 web-service test output |
| population, age 64 (yob 1962): BC risk to 80 | 6.039 | 6.07 | -0.5% | CanRisk 7.3.2 web-service test output |
| population, age 64: P(BRCA1 carrier) | 0.04279 | 0.0408 | +4.9% | CanRisk 7.3.2 web-service test output |
| population, age 64: P(BRCA2 carrier) | 0.1074 | 0.1056 | +1.7% | CanRisk 7.3.2 web-service test output |
| population, age 64: P(CHEK2 carrier) | 0.6831 | 0.6825 | +0.1% | CanRisk 7.3.2 web-service test output |
| population, age 64: lifetime OC (EOC model) | 1.671 | 1.762 | -5.2% | CanRisk ovarian model 3.2.2 web-service test output |
| population, age 64: OC risk to 80 (EOC model) | 0.9575 | 0.9893 | -3.2% | CanRisk ovarian model 3.2.2 web-service test output |
| untested, unknown FH: lifetime BC | 11.94 | 12 | -0.5% | Lee et al. JMG 2022, table s1 (born 1985) |
| untested, mother BC at 50: lifetime BC | 19.51 | 19.3 | +1.1% | Lee et al. JMG 2022, table s1 (born 1985) |
| BARD1, unknown FH: lifetime BC | 23.3 | 23.5 | -0.8% | Lee et al. JMG 2022, table s1 (born 1985) |
| BARD1, mother BC at 50: lifetime BC | 34.06 | 33.7 | +1.1% | Lee et al. JMG 2022, table s1 (born 1985) |
| RAD51C, unknown FH: lifetime BC | 22.13 | 22.3 | -0.8% | Lee et al. JMG 2022, table s1 (born 1985) |
| RAD51C, mother BC at 50: lifetime BC | 32.6 | 32.2 | +1.2% | Lee et al. JMG 2022, table s1 (born 1985) |
| RAD51D, unknown FH: lifetime BC | 20.64 | 20.8 | -0.8% | Lee et al. JMG 2022, table s1 (born 1985) |
| RAD51D, mother BC at 50: lifetime BC | 30.71 | 30.4 | +1.0% | Lee et al. JMG 2022, table s1 (born 1985) |
| negative panel, unknown FH: lifetime BC | 11.57 | 11.6 | -0.3% | Lee et al. JMG 2022, table s1 (born 1985) |
| negative panel, mother BC at 50: lifetime BC | 18.29 | 17.9 | +2.2% | Lee et al. JMG 2022, table s1 (born 1985) |
| BRCA2 carrier: BC risk 20-80 | 57.5 | 58 | -0.9% | Ficorella et al. BJC Rep 2024 (born 2000) |
| BRCA2 carrier, mother BC at 40: BC risk 20-80 | 72.05 | 72.5 | -0.6% | Ficorella et al. BJC Rep 2024 (born 2000) |
| PALB2 carrier: lifetime OC (EOC model) | 4.934 | 5 | -1.3% | Lee et al. JMG 2022 |

Other checks in `tests/` (`pytest`):

* The variable-elimination likelihood matches brute-force enumeration.
* Calibration reproduces the population and carrier incidences exactly.
* Averaging risks over the risk-factor or PRS distribution gives the population risk.
* PRS normalisation matches the CanRisk reference parameters.
* All example pedigrees of the CanRisk web services run.

## Layout

```
src/hbocrisk/
  params.py       gene RRs, allele frequencies, sensitivities, polygenic variances, pathology
  riskfactors.py  BC and EOC risk factors, continuous height, category amalgamation
  genetics.py     major-gene locus and hypergeometric polygenic model, PRS likelihood
  incidence.py    cohort/country incidence tables
  models.py       model specifications, calibration of genotype-specific incidences
  likelihood.py   penetrance, genetic tests, pathology, variable-elimination peeling
  api.py          carrier probabilities and risk calculations (CanRisk conventions)
  pedigree.py     data model, CanRisk/BOADICEA file I/O
  prs.py          PRS from VCF/dosages, reference PRS parameters
  cli.py, server.py
scripts/          incidence builders, validation
docs/MODEL.md     model specification and parameter provenance
```

## License

hbocrisk is distributed under the MIT License; see [LICENSE](LICENSE).

The reference PRS definitions in `src/hbocrisk/data/prs/` (variant positions, allele
frequencies and log odds ratios from the publications of the respective PRS) were taken
from the CanRisk web services (github.com/CCGE-BOADICEA/bws, © University of Cambridge).
The incidence tables are derived from IARC CI5plus and, for England, from ONS cancer
registration statistics (Open Government Licence).

BOADICEA and CanRisk are names of the University of Cambridge. hbocrisk is an
independent project and is not affiliated with or endorsed by the University of
Cambridge or Cambridge Enterprise.
