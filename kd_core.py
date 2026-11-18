"""Core logic of the qPCR knockdown analysis: file discovery, plate-layout parsing,
normalisation (relative expression to reference genes and a control sample) and QC.

Nothing in this module depends on Streamlit, so it can be tested on its own.
"""
from __future__ import annotations

import io
import os
import re
import tempfile
import zipfile
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

ROWS = list("ABCDEFGH")

HOUSEKEEPING = {
    "ACTB", "ACTIN", "GAPDH", "B2M", "HPRT", "HPRT1", "TBP", "PPIA", "RPLP0", "18S", "RN18S",
    "HMBS", "GUSB", "TUBB", "UBC", "YWHAZ", "PGK1", "SDHA", "RPL13A", "TFRC", "GAPD",
}
CONTROL_RE = re.compile(r"ctrl|control|mock|scr|scramble|neg|buf|untreated|\bwt\b", re.I)
NTC_RE = re.compile(r"^\s*(ntc|no[\s_-]*template.*|water|h2o|nfw)\s*$", re.I)
NORT_RE = re.compile(r"^\s*(-\s*rt|no[\s_-]*rt|minus[\s_-]*rt)\s*$", re.I)
DOSE_RE = re.compile(r"^(?P<group>.+?)[,;]\s*(?P<dose>[\d.]+\s*(?:nM|pM|µM|uM|mM|ng|µg|ug|%))\s*$", re.I)

RESULT_COLUMNS = ["Target", "Sample", "n_used", "n_wells", "Cq", "Cq_sd", "is_control", "dCt", "ddCt", "rel", "rel_lo",
                  "rel_hi", "sd_cycles", "knockdown_pct"]
SD_LIMIT = 0.3          # technical replicates with a larger SD are flagged
REF_RANGE_LIMIT = 1.0   # a reference gene whose mean Cq moves more than this between samples is flagged
NTC_MARGIN = 5.0        # NTC should be this many cycles later than the latest sample, or at/after Cq 35


# ----------------------------------------------------------------------------- small helpers

def normalize_well(val) -> str:
    """'A1', 'a01', 'SYBR:A1' -> 'A01'."""
    m = re.search(r"([A-Ha-h])\s*0?(\d{1,2})", str(val))
    return f"{m.group(1).upper()}{int(m.group(2)):02d}" if m else str(val).strip()


def gene_base(target: str) -> str:
    """'ACTB_h' -> 'ACTB'."""
    return re.split(r"[_\-\s]", str(target).strip())[0].upper()


def looks_like_reference(target: str) -> bool:
    return gene_base(target) in HOUSEKEEPING


def kind_of(sample: str) -> str:
    if NTC_RE.match(sample):
        return "NTC"
    if NORT_RE.match(sample):
        return "-RT"
    return "sample"


# ----------------------------------------------------------------------------- data classes

@dataclass
class Run:
    label: str
    cq: pd.DataFrame                      # Well, Cq, Target, Sample (the last two may be empty)
    amp: pd.DataFrame | None = None       # Cycle + one column per well
    melt: pd.DataFrame | None = None      # Temperature + one column per well
    hint: str = ""                        # text used to match the run to a layout file
    notes: list = field(default_factory=list)


@dataclass
class Layout:
    source: str                           # file name
    sheet: str
    targets: dict                         # well -> target (gene) name
    samples: dict                         # well -> sample name
    legend: dict = field(default_factory=dict)   # sample id -> description
    notes: list = field(default_factory=list)    # long free-text cells (instructions, warnings)

    @property
    def name(self) -> str:
        return f"{self.source} / {self.sheet}"


@dataclass
class Ingest:
    runs: list
    layouts: list
    warnings: list = field(default_factory=list)


# ----------------------------------------------------------------------------- reading instrument files

def _clean_matrix(path: str, axis: str) -> pd.DataFrame | None:
    df = pd.read_csv(path)
    df = df.loc[:, [c for c in df.columns if "unnamed" not in str(c).lower()]]
    ren = {c: (axis if axis.lower() in str(c).lower() else normalize_well(c)) for c in df.columns}
    df = df.rename(columns=ren)
    df = df.loc[:, ~df.columns.duplicated(keep="last")]
    return df if axis in df.columns else None


def read_cq_file(path: str):
    """Return (DataFrame[Well, Cq, Target, Sample], warnings)."""
    warnings = []
    df = pd.read_csv(path)
    if "Well" not in df.columns or "Cq" not in df.columns:
        raise ValueError(f"{os.path.basename(path)} has no 'Well' and 'Cq' columns")
    if "Fluor" in df.columns and df["Fluor"].nunique() > 1:
        first = df["Fluor"].dropna().iloc[0]
        df = df[df["Fluor"] == first]
        warnings.append(f"{os.path.basename(path)}: several fluorophores found, using {first} only.")
    out = pd.DataFrame({"Well": df["Well"].map(normalize_well), "Cq": pd.to_numeric(df["Cq"], errors="coerce")})
    for col in ("Target", "Sample"):
        out[col] = df[col].astype("string").str.strip().replace("", pd.NA).to_numpy() if col in df.columns else pd.NA
    return out.drop_duplicates("Well").reset_index(drop=True), warnings


# ----------------------------------------------------------------------------- plate layout (Excel)

def _cell_text(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float) and v.is_integer():
        v = int(v)
    return str(v).strip()


def _find_grids(ws):
    """Find 8x12 grids: a column that reads A..H downwards marks the top-left corner of a grid."""
    anchors = []
    for row in ws.iter_rows():
        for c in row:
            if _cell_text(c.value).upper() == "A":
                if all(_cell_text(ws.cell(c.row + i, c.column).value).upper() == L for i, L in enumerate(ROWS)):
                    anchors.append((c.row, c.column))
    return anchors


def parse_layouts(path: str) -> list:
    """Parse every sheet of an Excel file that holds a target grid and a sample grid."""
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True)
    layouts = []
    for ws in wb.worksheets:
        anchors = _find_grids(ws)
        if len(anchors) < 2:
            continue
        labelled = []
        for r, c in anchors:
            corner = _cell_text(ws.cell(r - 1, c).value).lower()
            labelled.append((corner, r, c))
        tgt = next((a for a in labelled if re.search(r"primer|target|gene|mix", a[0])), labelled[0])
        smp = next((a for a in labelled if "sample" in a[0] and a is not tgt), None)
        if smp is None:
            smp = next(a for a in labelled if a is not tgt)

        def grid(anchor):
            _, r0, c0 = anchor
            out = {}
            for i, letter in enumerate(ROWS):
                for j in range(12):
                    txt = _cell_text(ws.cell(r0 + i, c0 + 1 + j).value)
                    if txt and txt.lower() not in ("empty", "unused", "none", "nan", "-"):
                        out[f"{letter}{j + 1:02d}"] = txt
            return out

        targets, samples = grid(tgt), grid(smp)
        ids = set(samples.values())
        legend = {}
        long_text = []
        for row in ws.iter_rows():
            for c in row:
                txt = _cell_text(c.value)
                if not txt:
                    continue
                if txt in ids:
                    right = _cell_text(ws.cell(c.row, c.column + 1).value)
                    if right and right not in ids and not right.startswith("="):
                        legend.setdefault(txt, right)
                elif len(txt) >= 25 and not txt.startswith("="):
                    long_text.append(txt)
        legend_values = set(legend.values())
        notes = [t for t in long_text if t not in legend_values]
        layouts.append(Layout(os.path.basename(path), ws.title, targets, samples, legend, notes))
    return layouts


def layout_from_instrument(cq: pd.DataFrame) -> Layout | None:
    """Use Target/Sample names typed into the instrument software, if there are any."""
    ok = cq.dropna(subset=["Target", "Sample"])
    if ok.empty:
        return None
    return Layout("names in the Cq file", "", dict(zip(ok["Well"], ok["Target"].astype(str))),
                  dict(zip(ok["Well"], ok["Sample"].astype(str))))


# ----------------------------------------------------------------------------- upload ingestion

def _safe_extract(data: bytes, dest: str) -> None:
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        if sum(i.file_size for i in zf.infolist()) > 1_500_000_000:
            raise ValueError("ZIP archive is unreasonably large")
        root = os.path.realpath(dest)
        for m in zf.infolist():
            base = os.path.basename(m.filename)
            if m.is_dir() or m.filename.startswith("__MACOSX") or base.startswith(("._", "~$")):
                continue
            out = os.path.realpath(os.path.join(dest, m.filename))
            if not out.startswith(root + os.sep):
                continue
            os.makedirs(os.path.dirname(out), exist_ok=True)
            with zf.open(m) as src, open(out, "wb") as dst:
                dst.write(src.read())


def ingest(files: list) -> Ingest:
    """files: list of (name, bytes). ZIP archives, Bio-Rad CSV files and Excel layouts can be mixed."""
    runs, layouts, warnings = [], [], []
    with tempfile.TemporaryDirectory() as tmp:
        for i, (name, data) in enumerate(files):
            sub = os.path.join(tmp, str(i))
            os.makedirs(sub)
            if name.lower().endswith(".zip"):
                try:
                    _safe_extract(data, sub)
                except (zipfile.BadZipFile, ValueError) as e:
                    warnings.append(f"{name}: {e}")
            else:
                with open(os.path.join(sub, os.path.basename(name)), "wb") as f:
                    f.write(data)

        found = []
        for root, _, fnames in os.walk(tmp):
            for f in sorted(fnames):
                if not f.startswith(("._", "~$")):
                    found.append((root, f))

        for root, f in found:
            low = f.lower()
            path = os.path.join(root, f)
            if low.endswith(".csv") and "quantification cq results" in low:
                prefix = re.split(r"\s+-\s+", f)[0].strip()
                try:
                    cq, w = read_cq_file(path)
                except Exception as e:  # noqa: BLE001 - report any unreadable file to the user
                    warnings.append(str(e))
                    continue
                warnings += w
                amp = melt = None
                for r2, f2 in found:
                    if r2 != root or not f2.startswith(prefix) or not f2.lower().endswith(".csv"):
                        continue
                    if "quantification amplification results" in f2.lower():
                        amp = _clean_matrix(os.path.join(r2, f2), "Cycle")
                    elif "melt curve derivative results" in f2.lower():
                        melt = _clean_matrix(os.path.join(r2, f2), "Temperature")
                runs.append(Run(prefix, cq, amp, melt, hint=f"{os.path.basename(root)} {f}"))
            elif low.endswith(".xlsx"):
                try:
                    layouts += parse_layouts(path)
                except Exception as e:  # noqa: BLE001
                    warnings.append(f"{f}: could not be read ({e})")

    seen = {}
    for r in runs:  # make run labels unique
        seen[r.label] = seen.get(r.label, 0) + 1
        if seen[r.label] > 1:
            r.label = f"{r.label} ({seen[r.label]})"
    if not runs:
        warnings.append("No '... Quantification Cq Results.csv' file found.")
    return Ingest(runs, layouts, warnings)


# ----------------------------------------------------------------------------- matching runs to layouts

_MONTHS = {m: i for i, m in enumerate("jan feb mar apr may jun jul aug sep oct nov dec".split(), 1)}


def _dates(text: str) -> set:
    t = text.lower()
    out = set()
    for d, m in re.findall(r"(?<!\d)(\d{1,2})[\s_.\-]*(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)", t):
        out.add((int(d), _MONTHS[m]))
    for _, m, d in re.findall(r"(20\d{2})-(\d{2})-(\d{2})", t):
        out.add((int(d), int(m)))
    for d, m, _ in re.findall(r"(?<!\d)(\d{1,2})\.(\d{1,2})\.(\d{2,4})", t):
        out.add((int(d), int(m)))
    return out


def _plate_no(text: str):
    m = re.search(r"plate[\s_]*(\d+)", text.lower())
    return int(m.group(1)) if m else None


def match_score(run: Run, lay: Layout) -> int:
    lay_text = f"{lay.source} {lay.sheet}"
    score = 2 * len(_dates(run.hint) & _dates(lay_text))
    p1, p2 = _plate_no(run.hint), _plate_no(lay_text)
    if p1 is not None and p1 == p2:
        score += 3
    return score


def best_layout(run: Run, layouts: list) -> int | None:
    if not layouts:
        return None
    scores = [match_score(run, lay) for lay in layouts]
    return int(np.argmax(scores))


# ----------------------------------------------------------------------------- wells table

def build_wells(run: Run, layout: Layout | None) -> pd.DataFrame:
    """One row per used well: Run, Well, Target, Sample, Kind, Cq, WellKey."""
    cq = dict(zip(run.cq["Well"], run.cq["Cq"]))
    rows = []
    if layout is not None:
        for well in sorted(set(layout.targets) & set(layout.samples)):
            sample = layout.samples[well]
            rows.append((run.label, well, layout.targets[well], sample, kind_of(sample), cq.get(well, np.nan)))
    df = pd.DataFrame(rows, columns=["Run", "Well", "Target", "Sample", "Kind", "Cq"])
    df["Cq"] = pd.to_numeric(df["Cq"], errors="coerce")
    df["WellKey"] = df["Run"] + "|" + df["Well"]
    return df


def order_samples(samples, layouts: list) -> list:
    """Order samples like the layout legend (plate order); unknown ones follow in appearance order."""
    order = []
    for lay in layouts:
        for sid in lay.legend:
            if sid in samples and sid not in order:
                order.append(sid)
        for well in sorted(lay.samples, key=lambda w: (int(w[1:]), w[0])):
            sid = lay.samples[well]
            if sid in samples and sid not in order:
                order.append(sid)
    order += [s for s in samples if s not in order]
    return order


def find_experiments(wells: pd.DataFrame) -> list:
    """Group genes and samples into independent experiments: genes that share samples belong together."""
    s = wells[wells["Kind"] == "sample"]
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for t, smp in zip(s["Target"], s["Sample"]):
        parent[find(f"T:{t}")] = find(f"S:{smp}")
    groups = {}
    for t in dict.fromkeys(s["Target"]):
        groups.setdefault(find(f"T:{t}"), []).append(t)
    exps = []
    for targets in groups.values():
        sub = wells[wells["Target"].isin(targets)]
        samples = list(dict.fromkeys(sub.loc[sub["Kind"] == "sample", "Sample"]))
        exps.append({"targets": targets, "samples": samples, "wells": sub})
    return exps


# ----------------------------------------------------------------------------- quantification

def replicate_stats(wells: pd.DataFrame, excluded: set) -> pd.DataFrame:
    """Per target and sample: mean Cq, SD and number of replicates used."""
    s = wells[wells["Kind"] == "sample"].copy()
    s["Used"] = ~s["WellKey"].isin(excluded) & s["Cq"].notna()
    total = s.groupby(["Target", "Sample"]).size().rename("n_wells")
    used = s[s["Used"]].groupby(["Target", "Sample"])["Cq"].agg(mean_cq="mean", sd_cq="std", n_used="count")
    out = pd.concat([total, used], axis=1).reset_index()
    out["n_used"] = out["n_used"].fillna(0).astype(int)
    out["sd_cq"] = out["sd_cq"].fillna(0.0).where(out["n_used"] > 0)
    return out


def suggest_references(wells: pd.DataFrame, excluded: set = frozenset()):
    """Housekeeping-like targets, split into (stable, unstable) by how much their mean Cq moves between samples."""
    st = replicate_stats(wells, excluded)
    stable, unstable = [], []
    for t in dict.fromkeys(wells["Target"]):
        if not looks_like_reference(t):
            continue
        m = st[(st["Target"] == t) & (st["n_used"] > 0)]["mean_cq"]
        (stable if len(m) < 2 or m.max() - m.min() <= REF_RANGE_LIMIT else unstable).append(t)
    return stable, unstable


def reference_table(wells: pd.DataFrame, candidates: list, excluded: set) -> pd.DataFrame:
    st = replicate_stats(wells, excluded)
    rows = []
    for t in candidates:
        m = st[(st["Target"] == t) & (st["n_used"] > 0)]["mean_cq"]
        if m.empty:
            continue
        rng = float(m.max() - m.min())
        rows.append({"Gene": t, "Samples": len(m), "Mean Cq": float(m.mean()), "SD between samples": float(m.std(ddof=1)) if len(m) > 1 else 0.0,
                     "Range (cycles)": rng, "Status": "stable" if rng <= REF_RANGE_LIMIT else f"varies > {REF_RANGE_LIMIT:g} cycle"})
    return pd.DataFrame(rows)


def refs_from_notes(notes: list, targets: list, samples: list, always: bool = False) -> list:
    """Genes named in a layout note such as 'normalise to Actb only'.

    A note is used for an experiment only if it mentions one of its sample names (e.g. 'A2'),
    unless `always` is set (a single experiment on the whole plate).
    """
    tokens = {t for smp in samples for t in str(smp).split() if len(t) >= 2 and re.search(r"[A-Za-z]", t)}
    found = []
    for note in notes:
        if not always and not any(re.search(rf"(?<![\w-]){re.escape(t)}(?![\w-])", note) for t in tokens):
            continue
        for m in re.finditer(r"normali[sz]e\s+(?:\w+\s+)?to\s+([A-Za-z0-9]+(?:\s*(?:,|and|&)\s*[A-Za-z0-9]+)*)", note, re.I):
            names = {n.upper() for n in re.split(r"\s*(?:,|and|&)\s*", m.group(1))}
            found += [t for t in targets if gene_base(t) in names and t not in found]
    return found


def analyze(wells: pd.DataFrame, refs: list, calibrators: list, efficiency: dict | None = None,
            excluded: set = frozenset(), sample_order: list | None = None) -> pd.DataFrame:
    """Relative expression of every non-reference gene in every sample.

    Quantity of gene g in sample s relative to the calibrator: Q = E_g ** (Cq_cal - Cq_s). The reference
    factor is the geometric mean of Q over the reference genes, relative expression = Q_target / reference factor.
    With E = 2 for all genes this is the Livak 2^-ddCt method.
    """
    efficiency = efficiency or {}
    st = replicate_stats(wells, excluded)
    M = {(r.Target, r.Sample): r.mean_cq for r in st.itertuples() if r.n_used > 0}
    SD = {(r.Target, r.Sample): (r.sd_cq if not pd.isna(r.sd_cq) else 0.0) for r in st.itertuples() if r.n_used > 0}
    NU = {(r.Target, r.Sample): (r.n_used, r.n_wells) for r in st.itertuples()}
    E = lambda g: float(efficiency.get(g, 2.0))  # noqa: E731

    def cal(g):
        vals = [(M[(g, s)], SD[(g, s)]) for s in calibrators if (g, s) in M]
        if not vals:
            return np.nan, np.nan
        return float(np.mean([v[0] for v in vals])), float(np.sqrt(np.mean([v[1] ** 2 for v in vals])))

    cal_stats = {g: cal(g) for g in dict.fromkeys(wells["Target"])}
    samples = list(dict.fromkeys(wells.loc[wells["Kind"] == "sample", "Sample"]))
    if sample_order:
        samples = [s for s in sample_order if s in samples] + [s for s in samples if s not in sample_order]
    genes = [g for g in dict.fromkeys(wells["Target"]) if g not in refs]

    rows = []
    for g in genes:
        cg_mean, cg_sd = cal_stats[g]
        for s in samples:
            if (g, s) not in M or np.isnan(cg_mean):
                n_used, n_wells = NU.get((g, s), (0, 0))
                rows.append({"Target": g, "Sample": s, "n_used": n_used, "n_wells": n_wells})
                continue
            ref_ok = all((r, s) in M and not np.isnan(cal_stats[r][0]) for r in refs) and bool(refs)
            n_used, n_wells = NU[(g, s)]
            row = {"Target": g, "Sample": s, "n_used": n_used, "n_wells": n_wells, "Cq": M[(g, s)], "Cq_sd": SD[(g, s)],
                   "is_control": s in calibrators}
            if ref_ok:
                ref_cq = np.mean([M[(r, s)] for r in refs])
                ref_cal = np.mean([cal_stats[r][0] for r in refs])
                row["dCt"] = M[(g, s)] - ref_cq
                row["ddCt"] = row["dCt"] - (cg_mean - ref_cal)
                q_t = E(g) ** (cg_mean - M[(g, s)])
                log_rf = np.mean([(cal_stats[r][0] - M[(r, s)]) * np.log(E(r)) for r in refs])
                rel = q_t / np.exp(log_rf)
                var = SD[(g, s)] ** 2 + cg_sd ** 2 + sum(SD[(r, s)] ** 2 + cal_stats[r][1] ** 2 for r in refs) / len(refs) ** 2
                sd_cyc = float(np.sqrt(var))
                row.update({"rel": rel, "rel_lo": rel / E(g) ** sd_cyc, "rel_hi": rel * E(g) ** sd_cyc,
                            "sd_cycles": sd_cyc, "knockdown_pct": (1.0 - rel) * 100.0})
            rows.append(row)
    out = pd.DataFrame(rows, columns=None if rows else RESULT_COLUMNS)
    if out.empty:
        return out
    out["Sample"] = pd.Categorical(out["Sample"], categories=samples, ordered=True)
    return out.sort_values(["Target", "Sample"]).reset_index(drop=True)


# ----------------------------------------------------------------------------- QC

def replicate_flags(wells: pd.DataFrame, excluded: set, sd_limit: float = SD_LIMIT) -> pd.DataFrame:
    """Replicate groups whose SD is above the limit, with the single best well to exclude."""
    s = wells[(wells["Kind"] == "sample") & wells["Cq"].notna() & ~wells["WellKey"].isin(excluded)]
    rows = []
    for (t, smp), g in s.groupby(["Target", "Sample"]):
        if len(g) < 2:
            continue
        sd = g["Cq"].std()
        if sd <= sd_limit:
            continue
        suggestion = ""
        if len(g) >= 3:
            med = g["Cq"].median()
            far = (g["Cq"] - med).abs().idxmax()
            rest = g.drop(index=far)["Cq"]
            if rest.std() <= sd_limit:
                suggestion = f"exclude {g.loc[far, 'Well']} (Cq {g.loc[far, 'Cq']:.2f}): SD {sd:.2f} -> {rest.std():.2f}"
        rows.append({"Target": t, "Sample": smp, "n": len(g), "Mean Cq": g["Cq"].mean(), "SD": sd,
                     "Spread": g["Cq"].max() - g["Cq"].min(), "Suggestion": suggestion})
    return pd.DataFrame(rows, columns=["Target", "Sample", "n", "Mean Cq", "SD", "Spread", "Suggestion"])


def missing_wells(wells: pd.DataFrame) -> pd.DataFrame:
    m = wells[(wells["Kind"] == "sample") & wells["Cq"].isna()]
    return m[["Run", "Well", "Target", "Sample"]].reset_index(drop=True)


def ntc_table(wells: pd.DataFrame, excluded: set = frozenset()) -> pd.DataFrame:
    rows = []
    for t in dict.fromkeys(wells["Target"]):
        sub = wells[wells["Target"] == t]
        samples = sub[(sub["Kind"] == "sample") & sub["Cq"].notna() & ~sub["WellKey"].isin(excluded)]
        latest = samples.groupby("Sample")["Cq"].mean().max() if not samples.empty else np.nan
        for kind in ("NTC", "-RT"):
            ctl = sub[sub["Kind"] == kind]
            if ctl.empty:
                continue
            lowest = ctl["Cq"].min()
            if pd.isna(lowest):
                status = "Clean (no Cq)"
            elif lowest >= 35 or (not pd.isna(latest) and lowest >= latest + NTC_MARGIN):
                status = "Pass"
            else:
                status = "Check"
            rows.append({"Target": t, "Control": kind, "Wells": len(ctl), "With Cq": int(ctl["Cq"].notna().sum()),
                         "Lowest Cq": lowest, "Latest sample Cq": latest, "Status": status})
    return pd.DataFrame(rows, columns=["Target", "Control", "Wells", "With Cq", "Lowest Cq", "Latest sample Cq", "Status"])


# ----------------------------------------------------------------------------- plot helpers

def split_sample(sample: str, legend: dict):
    """(group, dose, description) used to draw grouped bars, e.g. 'siGENE-1', '10 nM'."""
    desc = legend.get(sample, "")
    m = DOSE_RE.match(desc)
    if m:
        return m.group("group").strip(), m.group("dose").replace(" ", "").replace("uM", "µM"), desc
    return (desc or sample), "", desc


def plot_frame(results: pd.DataFrame, target: str, legend: dict) -> pd.DataFrame:
    """Rows of one gene ready for a grouped bar chart."""
    d = results[(results["Target"] == target) & results["rel"].notna()].copy() if "rel" in results else pd.DataFrame()
    if d.empty:
        return d
    parts = [split_sample(str(s), legend) for s in d["Sample"]]
    d["Group"] = [p[0] for p in parts]
    d["Dose"] = [p[1] for p in parts]
    d["Description"] = [p[2] for p in parts]
    return d


# ----------------------------------------------------------------------------- synthetic demo

def demo_ingest(seed: int = 1) -> Ingest:
    """Two synthetic plates that mirror a typical siRNA experiment (target gene on one plate, reference on the other)."""
    rng = np.random.default_rng(seed)
    samples = ["1 Ctrl", "2 Mock", "3 Scr1", "4 Scr2", "5 si1-10", "6 si1-30", "7 si2-10", "8 si2-30"]
    legend = {"1 Ctrl": "Control", "2 Mock": "Transfection reagent only", "3 Scr1": "Scramble 1", "4 Scr2": "Scramble 2",
              "5 si1-10": "siGENE-1, 10 nM", "6 si1-30": "siGENE-1, 30 nM", "7 si2-10": "siGENE-2, 10 nM", "8 si2-30": "siGENE-2, 30 nM"}
    true_rel = dict(zip(samples, [1.0, 1.05, 0.98, 1.02, 0.45, 0.20, 0.70, 0.35]))
    cols = [f"{r}{c:02d}" for c in range(1, 13) for r in ROWS]  # column-major, like pipetting

    def plate(genes_by_cols, base, effect):
        targets, smp, cq = {}, {}, {}
        for gi, (gene, cstart) in enumerate(genes_by_cols):
            wells = [f"{r}{c:02d}" for c in range(cstart, cstart + 2) for r in ROWS]
            for i, s in enumerate(samples):
                for rep in range(3):
                    w = wells[i * 3 + rep] if i * 3 + rep < 16 else None
                    if w is None:
                        continue
                    targets[w], smp[w] = gene, s
                    shift = -np.log2(true_rel[s]) if effect.get(gene) else 0.0
                    cq[w] = base[gene] + shift + rng.normal(0, 0.08)
        return targets, smp, cq

    # only 5 samples fit per 2-column block at 3 replicates, so use 3 blocks of 2 columns per gene
    def plate2(gene_blocks, base, effect):
        targets, smp, cq = {}, {}, {}
        for gene, blocks in gene_blocks.items():
            k = 0
            for cstart in blocks:
                wells = [f"{r}{c:02d}" for c in range(cstart, cstart + 2) for r in ROWS]
                for i in range(4):
                    s = samples[k * 4 + i]
                    for rep in range(3):
                        w = wells[i * 3 + rep]
                        targets[w], smp[w] = gene, s
                        shift = -np.log2(true_rel[s]) + rng.normal(0, 0.05) if effect.get(gene) else rng.normal(0, 0.05)
                        cq[w] = base[gene] + shift + rng.normal(0, 0.08)
                k += 1
        return targets, smp, cq

    out_runs, out_layouts = [], []
    specs = [
        ("demo plate 2", {"ACTB_h": [1, 3]}, {"ACTB_h": 17.0}, {}),
        ("demo plate 1", {"GENE_h": [1, 3], "GAPDH_h": [7, 9]}, {"GENE_h": 26.0, "GAPDH_h": 19.0}, {"GENE_h": True}),
    ]
    for label, gene_blocks, base, effect in specs:
        targets, smp, cq = plate2(gene_blocks, base, effect)
        for gene, cstarts in gene_blocks.items():  # one NTC per gene
            w = f"E{cstarts[0] + 1:02d}"
            if w not in targets:
                targets[w], smp[w], cq[w] = gene, "NTC", np.nan
        table = pd.DataFrame({"Well": cols})
        table["Cq"] = [cq.get(w, np.nan) for w in cols]
        table["Target"] = pd.NA
        table["Sample"] = pd.NA
        out_runs.append(Run(label, table, None, None, hint=label))
        out_layouts.append(Layout(f"{label}.xlsx", "Sheet1", targets, smp, dict(legend), []))
    return Ingest(out_runs, out_layouts, [])


# ----------------------------------------------------------------------------- curves for QC plots

def curves_for(wells: pd.DataFrame, runs: dict, target: str, kind: str = "melt"):
    """(x values, x label, [(well key, sample, well kind, y values)]) for one gene; kind is 'melt' or 'amp'."""
    out, x, xcol = [], None, None
    for run_label, g in wells[wells["Target"] == target].groupby("Run"):
        run = runs.get(run_label)
        df = None if run is None else (run.melt if kind == "melt" else run.amp)
        if df is None:
            continue
        xcol = df.columns[0] if x is None else xcol
        if x is None:
            x = df[df.columns[0]].to_numpy()
        for r in g.itertuples():
            if r.Well in df.columns:
                out.append((r.WellKey, r.Sample, r.Kind, df[r.Well].to_numpy()))
    return x, xcol, out
