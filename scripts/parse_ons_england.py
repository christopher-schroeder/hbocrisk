"""Parse ONS Cancer Registration Statistics England 'Table 1' (counts by site, sex, age)."""
import sys, re, json
import pandas as pd
AGECOLS = ["Under 1","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49",
           "50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90 and over"]
SITES = {"C50": None, "C56": None, "C57": None, "C61": None, "C25": None}

def norm(s):
    s = re.sub(r"\s+", " ", str(s).replace("\xad", "-")).strip()
    s = re.sub(r"^(\d+) to (\d+)$", r"\1-\2", s)
    s = {"90+": "90 and over", "90 & over": "90 and over", "1 - 4": "1-4",
         "85 and": "85 and over", "85+": "85 and over"}.get(s, s)
    return s

def parse(path, year):
    xl = pd.ExcelFile(path)
    sheets = [s for s in xl.sheet_names if s.strip() in ("Table 1 Data",)] + \
             [s for s in xl.sheet_names if s.strip() in ("Table 1",)] + list(xl.sheet_names)
    hdr = None
    for sheet in sheets:
        df = pd.read_excel(path, sheet_name=sheet, header=None)
        for i in range(min(40, len(df))):
            vals = [norm(v) for v in df.iloc[i].tolist()]
            if "45-49" in vals and "50-54" in vals and "1-4" in vals:
                hdr = i; break
        if hdr is not None:
            print(" using sheet", sheet)
            break
    vals = [norm(v) for v in df.iloc[hdr].tolist()]
    agecol = {a: vals.index(a) for a in AGECOLS + ["85 and over"] if a in vals}
    bands = [["Under 1", "1-4"]] + [[a] for a in AGECOLS[2:18]] + \
        [["85-89", "90 and over"] if "85-89" in agecol else ["85 and over"]]
    out = {}
    cur_code = None
    for i in range(hdr + 1, len(df)):
        row = [norm(v) for v in df.iloc[i].tolist()]
        if "Number of Deaths" in " ".join(row) or "Deaths" == row[0]:
            continue
        codes = [v for v in row[:4] if re.fullmatch(r"C\d\d(-C\d\d)?(\.\d)?", v)]
        if codes:
            cur_code = codes[0]
        sexes = [{"M": "Males", "F": "Females"}.get(v, v) for v in row[:6]
                 if v in ("Males", "Females", "Persons", "M", "F")]
        if cur_code in SITES and sexes:
            try:
                def num(a):
                    v = row[agecol[a]].replace(",", "")
                    return 0.0 if v in ("nan", "-", ":", "") else float(v)
                counts = [sum(num(a) for a in b) for b in bands]
            except (ValueError, KeyError):
                continue
            key = f"{cur_code}_{sexes[0]}"
            if key not in out:  # first occurrence = registrations (deaths follow later)
                out[key] = counts
    return out

res = {}
for y in sys.argv[1:]:
    y = int(y)
    f = {2017: "2017.xlsx", 2008: "z2008/No.39, 2008/Table 1.xls"}.get(y, f"{y}.xls")
    res[y] = parse(f, y)
    print(y, {k: int(sum(v)) for k, v in res[y].items()})
json.dump(res, open("ons_counts.json", "w"))
