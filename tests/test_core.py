import io
import zipfile

import numpy as np
import pandas as pd
import pytest

import kd_core as kd
import kd_report


@pytest.fixture(scope="module")
def demo():
    ing = kd.demo_ingest()
    wells = pd.concat([kd.build_wells(r, lay) for r, lay in zip(ing.runs, ing.layouts)], ignore_index=True)
    return ing, wells


def test_normalize_well():
    assert kd.normalize_well("A1") == "A01"
    assert kd.normalize_well("SYBR:h12") == "H12"


def test_demo_recovers_known_fold_changes(demo):
    _, wells = demo
    exp = kd.find_experiments(wells)
    assert len(exp) == 1 and set(exp[0]["targets"]) == {"ACTB_h", "GENE_h", "GAPDH_h"}
    res = kd.analyze(exp[0]["wells"], ["ACTB_h", "GAPDH_h"], ["1 Ctrl"])
    got = res[res["Target"] == "GENE_h"].set_index("Sample")["rel"]
    true = {"5 si1-10": 0.45, "6 si1-30": 0.20, "7 si2-10": 0.70, "8 si2-30": 0.35}
    for s, v in true.items():
        assert got[s] == pytest.approx(v, rel=0.12)
    assert got["1 Ctrl"] == pytest.approx(1.0)


def test_livak_formula_with_one_reference():
    rows = []
    for s, t_cq, r_cq in [("ctrl", 25.0, 18.0), ("kd", 27.0, 18.0)]:
        for rep in range(3):
            rows += [("r", f"A0{rep + 1}", "GOI", s, "sample", t_cq), ("r", f"B0{rep + 1}", "REF", s, "sample", r_cq)]
    wells = pd.DataFrame(rows, columns=["Run", "Well", "Target", "Sample", "Kind", "Cq"])
    wells["WellKey"] = wells["Run"] + "|" + wells["Well"]
    res = kd.analyze(wells, ["REF"], ["ctrl"]).set_index("Sample")
    assert res.loc["kd", "ddCt"] == pytest.approx(2.0)
    assert res.loc["kd", "rel"] == pytest.approx(0.25)          # 2^-2
    assert res.loc["kd", "knockdown_pct"] == pytest.approx(75.0)
    res90 = kd.analyze(wells, ["REF"], ["ctrl"], {"GOI": 1.9, "REF": 1.9}).set_index("Sample")
    assert res90.loc["kd", "rel"] == pytest.approx(1.9 ** -2)   # efficiency-corrected


def test_excluding_a_well_changes_the_mean(demo):
    _, wells = demo
    sub = wells[wells["Target"] == "GENE_h"]
    key = sub[(sub["Sample"] == "1 Ctrl")]["WellKey"].iloc[0]
    full = kd.replicate_stats(sub, set())
    cut = kd.replicate_stats(sub, {key})
    a = full[full["Sample"] == "1 Ctrl"].iloc[0]
    b = cut[cut["Sample"] == "1 Ctrl"].iloc[0]
    assert b["n_used"] == a["n_used"] - 1


def test_reference_stability_flags_a_moving_gene(demo):
    _, wells = demo
    w = wells.copy()
    w.loc[(w["Target"] == "GAPDH_h") & (w["Sample"] == "8 si2-30"), "Cq"] += 3.0
    stable, unstable = kd.suggest_references(w)
    assert "GAPDH_h" in unstable and "ACTB_h" in stable


def test_refs_from_notes_only_for_matching_experiment():
    notes = ["A20: GAPDH is knocked down in A2, so normalise to Actb only."]
    targets = ["ACTB_m", "GAPDH_m", "CAMK1D_m"]
    assert kd.refs_from_notes(notes, targets, ["A1 FAM-Scr", "A2 siGAPDH"]) == ["ACTB_m"]
    assert kd.refs_from_notes(notes, ["ACTB_h"], ["1 Ctrl", "3 Scr1"]) == []
    assert kd.refs_from_notes(notes, ["ACTB_h"], ["1 Ctrl"], always=True) == ["ACTB_h"]


def test_template_layout_roundtrip(tmp_path):
    path = tmp_path / "t.xlsx"
    path.write_bytes(kd_report.template_xlsx())
    lay = kd.parse_layouts(str(path))
    assert len(lay) == 1
    assert lay[0].targets["A01"] == "GENE_h" and lay[0].samples["A01"] == "1 Ctrl"
    assert lay[0].legend == {"1 Ctrl": "Control"}  # only samples that are on the plate


def test_ingest_skips_macos_junk_and_matches_layout(tmp_path):
    ing = kd.demo_ingest()
    cq = ing.runs[0].cq.copy()
    cq["Target"] = cq["Sample"] = ""
    csv = cq.assign(Fluor="SYBR").to_csv(index=False).encode()
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("run 17nov26/x - Quantification Cq Results.csv", csv)
        z.writestr("__MACOSX/run 17nov26/._x - Quantification Cq Results.csv", b"\x00\x05binary")
    layout = tmp_path / "plate_17Nov26.xlsx"
    layout.write_bytes(kd_report.template_xlsx())
    out = kd.ingest([("a.zip", buf.getvalue()), ("plate_17Nov26.xlsx", layout.read_bytes())])
    assert len(out.runs) == 1 and len(out.layouts) == 1 and not out.warnings


def test_layout_matching_uses_date_and_plate_number():
    run = kd.Run("x plate 2", pd.DataFrame(columns=["Well", "Cq"]), hint="qpcr 17nov26 x plate 2 - Quantification Cq Results.csv")
    a = kd.Layout("layout_18Nov26.xlsx", "Plate1", {}, {})
    b = kd.Layout("layout_17Nov26.xlsx", "Plate2", {}, {})
    assert kd.best_layout(run, [a, b]) == 1


def test_reports_build(demo):
    ing, wells = demo
    exp = kd.find_experiments(wells)[0]
    refs, cals = ["ACTB_h", "GAPDH_h"], ["3 Scr1", "4 Scr2"]
    legend = {k: v for lay in ing.layouts for k, v in lay.legend.items()}
    samples = kd.order_samples(exp["samples"], ing.layouts)
    e = {"name": "demo", "targets": exp["targets"], "samples": samples, "refs": refs, "calibrators": cals, "efficiency": {},
         "excluded": set(), "results": kd.analyze(exp["wells"], refs, cals, sample_order=samples),
         "ref_table": kd.reference_table(exp["wells"], refs, set()), "flags": kd.replicate_flags(exp["wells"], set()),
         "ntc": kd.ntc_table(exp["wells"]), "legend": legend, "notes": [], "warnings": [], "wells": exp["wells"], "sd_limit": kd.SD_LIMIT}
    runs = {r.label: r for r in ing.runs}
    assert kd_report.build_pdf("demo", [e], wells, runs)[:4] == b"%PDF"
    assert kd_report.build_excel([e], wells, list(runs), set())[:2] == b"PK"
