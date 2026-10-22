"""Pedigree data model and readers for CanRisk (v1-v4) and BOADICEA v4 pedigree files.

File format specification: CanRisk "Pedigree Data v4 File Format Specification"
(canrisk.atlassian.net). Example::

    ##CanRisk 4.0
    ##menarche=13
    ##PRS_BC=alpha=0.501,zscore=1.2
    ##FamID Name Target IndivID FathID MothID Sex MZtwin Dead Age Yob BC1 BC2 OC PRO PAN Ashkn BRCA1 ...
    FAM1 me 1 ch1 m21 f21 F 0 0 35 1990 0 0 0 0 0 0 S:N 0:0 ...
"""
from dataclasses import dataclass, field, asdict
import re

CANCERS = ("bc1", "bc2", "oc", "prc", "pac")
PATHOLOGY_FIELDS = ("er", "pr", "her2", "ck14", "ck56")
GENES = ("BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1",
         "HOXB13")

COLUMNS = {
    1: ["FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead",
        "Age", "Yob", "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn", "BRCA1", "BRCA2", "PALB2",
        "ATM", "CHEK2", "RAD51D", "RAD51C", "BRIP1", "ER:PR:HER2:CK14:CK56"],
    2: ["FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead",
        "Age", "Yob", "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn", "BRCA1", "BRCA2", "PALB2",
        "ATM", "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1", "ER:PR:HER2:CK14:CK56"],
    4: ["FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead",
        "Age", "Yob", "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn", "BRCA1", "BRCA2", "PALB2",
        "ATM", "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1", "HOXB13",
        "ER:PR:HER2:CK14:CK56"],
}
COLUMNS[3] = COLUMNS[2]

BWA_COLUMNS = ["FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin",
               "Dead", "Age", "Yob", "1stBrCa", "2ndBrCa", "OvCa", "ProCa", "PanCa", "Ashkn",
               "BRCA1t", "BRCA1r", "BRCA2t", "BRCA2r", "PALB2t", "PALB2r", "ATMt", "ATMr",
               "CHEK2t", "CHEK2r", "ER", "PR", "HER2", "CK14", "CK56"]


class PedigreeError(ValueError):
    pass


@dataclass
class GeneticTest:
    type: str = "0"       # '0' untested, 'S' mutation search / panel, 'T' direct test
    result: str = "0"     # '0' untested, 'N' negative, 'P' positive

    @property
    def tested(self):
        return self.type in ("S", "T") and self.result in ("N", "P")


@dataclass
class Person:
    pid: str
    sex: str                              # 'F' or 'M'
    father: str = "0"
    mother: str = "0"
    name: str = ""
    target: bool = False
    mztwin: str = "0"
    dead: bool = False
    age: int = 0                          # age at last follow-up / death; 0 unknown
    yob: int = 0                          # year of birth; 0 unknown
    ashkenazi: bool = False
    # ages at diagnosis (None unaffected); 'AU' affected at unknown age
    bc1: object = None
    bc2: object = None
    oc: object = None
    prc: object = None
    pac: object = None
    tests: dict = field(default_factory=dict)        # gene -> GeneticTest
    pathology: dict = field(default_factory=dict)    # er/pr/her2/ck14/ck56 -> '0'|'N'|'P'

    @property
    def is_founder(self):
        return self.father in ("0", "", None) and self.mother in ("0", "", None)

    def cancer_age(self, c):
        """Age at diagnosis as int, or None. 'AU' (age unknown) is set to the age at last
        follow-up, as in CanRisk."""
        v = getattr(self, c)
        if v is None:
            return None
        if v == "AU":
            return self.age if self.age else None
        return int(v)

    def affected(self):
        return any(getattr(self, c) is not None for c in CANCERS)

    def pathology_code(self):
        """Breast tumour pathology code as used by BOADICEA (see params.PATHOLOGY_CODES)."""
        t = {k: self.pathology.get(k, "0") for k in PATHOLOGY_FIELDS}
        code = -1
        if t["er"] == "P":
            code = 0
        if t["er"] == "N":
            if (t["pr"] == "0" or t["her2"] == "0") and t["pr"] != "P" and t["her2"] != "P":
                code = 1
            elif t["pr"] == "P" or t["her2"] == "P":
                code = 2
            elif t["pr"] == "N" and t["her2"] == "N":
                if t["ck14"] == "0" or t["ck56"] == "0":
                    code = 3
                elif t["ck14"] == "N" and t["ck56"] == "N":
                    code = 4
                elif t["ck14"] == "P" and t["ck56"] == "P":
                    code = 6
                else:
                    code = 5
        return code

    def to_dict(self):
        d = asdict(self)
        d["tests"] = {g: [t.type, t.result] for g, t in self.tests.items()}
        return d

    @classmethod
    def from_dict(cls, d):
        d = dict(d)
        tests = {}
        for g, v in (d.pop("tests", None) or {}).items():
            if isinstance(v, GeneticTest):
                tests[g.upper()] = v
            elif isinstance(v, dict):
                tests[g.upper()] = GeneticTest(v.get("type", "0"), v.get("result", "0"))
            else:
                tests[g.upper()] = GeneticTest(*v)
        known = set(cls.__dataclass_fields__)
        unknown = set(d) - known
        if unknown:
            raise PedigreeError(f"unknown person fields: {sorted(unknown)}")
        p = cls(**d)
        p.tests = tests
        return p


@dataclass
class Pedigree:
    famid: str
    people: list
    risk_factors: dict = field(default_factory=dict)  # proband's risk factors (header)
    prs: dict = field(default_factory=dict)           # 'BC'/'OC'/'PC' -> (alpha, zscore)
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.by_id = {p.pid: p for p in self.people}

    @property
    def proband(self):
        t = [p for p in self.people if p.target]
        if len(t) != 1:
            raise PedigreeError(f"pedigree {self.famid} must have exactly one target "
                                f"individual (found {len(t)})")
        return t[0]

    def is_ashkenazi(self):
        return any(p.ashkenazi for p in self.people)

    def children(self, pid):
        return [p for p in self.people if p.father == pid or p.mother == pid]

    def validate(self):
        """Structural checks (as in the CanRisk web services)."""
        ids = [p.pid for p in self.people]
        if len(ids) != len(set(ids)):
            raise PedigreeError("duplicate individual identifiers")
        self.by_id = {p.pid: p for p in self.people}
        proband = self.proband
        for p in self.people:
            if p.sex not in ("F", "M"):
                raise PedigreeError(f"{p.pid}: sex must be 'F' or 'M'")
            fa, mo = p.father not in ("0", ""), p.mother not in ("0", "")
            if fa != mo:
                raise PedigreeError(f"{p.pid}: both parents or no parents must be given")
            if fa:
                if p.father not in self.by_id or p.mother not in self.by_id:
                    raise PedigreeError(f"{p.pid}: parent missing from pedigree")
                if self.by_id[p.father].sex != "M":
                    raise PedigreeError(f"{p.pid}: father must be male")
                if self.by_id[p.mother].sex != "F":
                    raise PedigreeError(f"{p.pid}: mother must be female")
            for c in CANCERS:
                a = getattr(p, c)
                if a is None or a == "AU":
                    continue
                if p.age and int(a) > p.age:
                    raise PedigreeError(f"{p.pid}: age at diagnosis ({c}) exceeds age")
            if p.sex == "M" and p.oc is not None:
                raise PedigreeError(f"{p.pid}: male with ovarian cancer")
            if p.sex == "F" and p.prc is not None:
                raise PedigreeError(f"{p.pid}: female with prostate cancer")
            if p.bc2 is not None and p.bc1 is None:
                raise PedigreeError(f"{p.pid}: contralateral breast cancer without first")
        # MZ twins
        twins = {}
        for p in self.people:
            if p.mztwin not in ("0", "", None):
                twins.setdefault(p.mztwin, []).append(p)
        for k, tw in twins.items():
            if len(tw) != 2:
                raise PedigreeError(f"MZ twin identifier {k} must appear exactly twice")
            a, b = tw
            if (a.father, a.mother) != (b.father, b.mother) or a.sex != b.sex:
                raise PedigreeError(f"MZ twins {a.pid}/{b.pid} must share parents and sex")
            if a.yob != b.yob:
                raise PedigreeError(f"MZ twins {a.pid}/{b.pid} must have the same year of birth")
        # connectivity to proband
        connected = {proband.pid}
        changed = True
        while changed:
            changed = False
            for p in self.people:
                if p.pid in connected:
                    for par in (p.father, p.mother):
                        if par not in ("0", "") and par not in connected:
                            connected.add(par)
                            changed = True
                elif p.father in connected or p.mother in connected:
                    connected.add(p.pid)
                    changed = True
        missing = [p.pid for p in self.people if p.pid not in connected]
        if missing:
            raise PedigreeError(f"individuals not connected to the target: {missing}")
        if not proband.yob or not proband.age:
            raise PedigreeError("the target must have an age and a year of birth")
        return True

    def to_dict(self):
        return {"famid": self.famid, "people": [p.to_dict() for p in self.people],
                "risk_factors": dict(self.risk_factors),
                "prs": {k: list(v) for k, v in self.prs.items()}}

    @classmethod
    def from_dict(cls, d):
        prs = {k.upper(): tuple(v) if not isinstance(v, dict) else (v["alpha"], v["zscore"])
               for k, v in (d.get("prs") or {}).items()}
        return cls(d.get("famid", "FAM"), [Person.from_dict(p) for p in d["people"]],
                   dict(d.get("risk_factors") or {}), prs)


# --------------------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------------------

def _age(v):
    v = str(v).strip()
    if v in ("0", "-1", "", "NA"):
        return None
    if v.upper() == "AU":
        return "AU"
    return int(v)


def _parse_prs(val):
    alpha = z = None
    for part in val.replace(" ", "").split(","):
        if "=" in part:
            k, v = part.split("=", 1)
            k = k.lower()
            if "alpha" in k:
                alpha = float(v)
            elif "zscore" in k or "beta" in k:
                z = float(v)
    if alpha is None or z is None:
        return None
    return (alpha, z)


def read_canrisk(text):
    """Parse a CanRisk (v1-v4) or BOADICEA v4 pedigree file. Returns a list of Pedigree."""
    lines = [ln.rstrip("\n\r") for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise PedigreeError("empty pedigree file")
    version = None
    m = re.match(r"^##\s*CanRisk\s+(\d)", lines[0], re.I)
    if m:
        version = int(m.group(1))
    elif lines[0].startswith("BOADICEA import pedigree file format 4.0"):
        version = "bwa"
    else:
        raise PedigreeError("first line must be '##CanRisk <version>' or a BOADICEA v4 header")

    pedigrees = []
    header = {}
    prs = {}
    records = []

    def flush():
        if records:
            pedigrees.append(_build(records, header, prs, version))

    for ln in lines[1:]:
        if ln.startswith("##"):
            body = ln[2:].strip()
            if body.lower().startswith("famid"):
                continue
            if records:                       # new family block
                flush()
                records.clear()
                header, prs = {}, {}
            if "=" in body:
                k, v = body.split("=", 1)
                k = k.strip().lower()
                if k.startswith("prs_"):
                    p = _parse_prs(v)
                    if p is not None:
                        prs[k[4:].upper()] = p
                else:
                    header[k] = v.strip()
            continue
        if version == "bwa" and ln.startswith("FamID"):
            continue
        records.append(re.split(r"\s+", ln.strip()))
    flush()
    return pedigrees


def _build(records, header, prs, version):
    cols = BWA_COLUMNS if version == "bwa" else COLUMNS[version]
    people = []
    for r in records:
        if len(r) != len(cols):
            raise PedigreeError(f"expected {len(cols)} columns, found {len(r)}: {' '.join(r)}")
        d = dict(zip(cols, r))
        tests, path = {}, {}
        if version == "bwa":
            for g in ("BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2"):
                tests[g] = GeneticTest(d[g + "t"], d[g + "r"])
            for k, c in zip(PATHOLOGY_FIELDS, ("ER", "PR", "HER2", "CK14", "CK56")):
                path[k] = d[c]
            bc1, bc2, oc, prc, pac = (d["1stBrCa"], d["2ndBrCa"], d["OvCa"], d["ProCa"],
                                      d["PanCa"])
        else:
            for g in GENES:
                if g in d:
                    t, res = (d[g].split(":") + ["0"])[:2]
                    tests[g] = GeneticTest(t, res)
            for k, v in zip(PATHOLOGY_FIELDS, d["ER:PR:HER2:CK14:CK56"].split(":")):
                path[k] = v
            bc1, bc2, oc, prc, pac = d["BC1"], d["BC2"], d["OC"], d["PRO"], d["PAN"]
        people.append(Person(
            pid=d["IndivID"], sex=d["Sex"].upper(), father=d["FathID"], mother=d["MothID"],
            name=d["Name"], target=d["Target"] == "1", mztwin=d["MZtwin"],
            dead=d["Dead"] == "1", age=int(d["Age"]), yob=int(d["Yob"]),
            ashkenazi=d["Ashkn"] == "1", bc1=_age(bc1), bc2=_age(bc2), oc=_age(oc),
            prc=_age(prc), pac=_age(pac), tests=tests, pathology=path))
    famid = records[0][0]
    return Pedigree(famid, people, dict(header), dict(prs), {"format": version})


def read_canrisk_file(path):
    with open(path) as fh:
        return read_canrisk(fh.read())


def write_canrisk(pedigree):
    """Serialise a Pedigree in CanRisk v4 format."""
    out = ["##CanRisk 4.0"]
    for k, v in pedigree.risk_factors.items():
        out.append(f"##{k}={v}")
    for k, (a, z) in pedigree.prs.items():
        out.append(f"##PRS_{k}=alpha={a},zscore={z}")
    out.append("##" + "\t".join(COLUMNS[4]))
    for p in pedigree.people:
        def a(v):
            return "0" if v is None else str(v)
        row = [pedigree.famid, p.name or p.pid, "1" if p.target else "0", p.pid,
               p.father or "0", p.mother or "0", p.sex, p.mztwin or "0",
               "1" if p.dead else "0", str(p.age), str(p.yob), a(p.bc1), a(p.bc2),
               a(p.oc), a(p.prc), a(p.pac), "1" if p.ashkenazi else "0"]
        for g in COLUMNS[4][17:27]:
            t = p.tests.get(g, GeneticTest())
            row.append(f"{t.type}:{t.result}")
        row.append(":".join(p.pathology.get(k, "0") for k in PATHOLOGY_FIELDS))
        out.append("\t".join(row))
    return "\n".join(out) + "\n"
