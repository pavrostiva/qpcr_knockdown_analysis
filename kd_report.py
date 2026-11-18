"""Excel and PDF reports for the knockdown analysis."""
from __future__ import annotations

import io
from datetime import date
from xml.sax.saxutils import escape

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from kd_core import ROWS, curves_for, plot_frame

NAVY = colors.HexColor("#1F4E79")
GREY = colors.HexColor("#6B7280")
LIGHT = colors.HexColor("#F3F4F6")
DOSE_COLORS = ["#1F4E79", "#4F8FCB", "#9DC3E6", "#C9A227", "#7A5195"]
CONTROL_COLOR = "#9CA3AF"

_ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=_ss["Title"], fontSize=21, leading=25, textColor=NAVY, alignment=0, spaceAfter=2)
H2 = ParagraphStyle("H2", parent=_ss["Heading2"], fontSize=13, leading=16, textColor=NAVY, spaceBefore=6, spaceAfter=4)
BODY = ParagraphStyle("Body", parent=_ss["BodyText"], fontSize=9, leading=12)
SMALL = ParagraphStyle("Small", parent=BODY, fontSize=8, leading=10, textColor=GREY)
CELL = ParagraphStyle("Cell", parent=BODY, fontSize=8, leading=9.5)
CELL_C = ParagraphStyle("CellC", parent=CELL, alignment=TA_CENTER)
HEAD = ParagraphStyle("Head", parent=CELL_C, textColor=colors.white, fontName="Helvetica-Bold")
PAGE = landscape(A4)
MARGIN = 1.4 * cm
CONTENT_W = PAGE[0] - 2 * MARGIN


def _f(v, spec):
    return "n.d." if v is None or (isinstance(v, float) and np.isnan(v)) else format(v, spec)


# ----------------------------------------------------------------------------- shared tables

def results_table(exp: dict) -> pd.DataFrame:
    """Readable results of one experiment (one row per gene and sample)."""
    r = exp["results"].copy()
    if r.empty:
        return r
    r["Description"] = [exp["legend"].get(str(s), "") for s in r["Sample"]]
    r["Sample"] = r["Sample"].astype(str)
    r["Role"] = np.where(r["is_control"].fillna(False), "control (=1)", "")
    cols = ["Target", "Sample", "Description", "Role", "n_used", "n_wells", "Cq", "Cq_sd", "dCt", "ddCt",
            "rel", "rel_lo", "rel_hi", "knockdown_pct"]
    r = r.reindex(columns=cols)
    return r.rename(columns={"n_used": "Replicates used", "n_wells": "Replicates on plate", "Cq": "Mean Cq", "Cq_sd": "SD Cq",
                             "dCt": "dCt", "ddCt": "ddCt", "rel": "Relative expression", "rel_lo": "Rel. expr. low (-SD)",
                             "rel_hi": "Rel. expr. high (+SD)", "knockdown_pct": "Knockdown (%)"})


def well_map(wells: pd.DataFrame, legend: dict, excluded: set) -> pd.DataFrame:
    """Every used well with the sample and gene assigned to it."""
    w = wells.copy()
    w["Row"] = w["Well"].str[0]
    w["Col"] = w["Well"].str[1:].astype(int)
    w["Description"] = [legend.get(s, "") for s in w["Sample"]]
    w["Excluded"] = np.where(w["WellKey"].isin(excluded), "yes", "")
    return w[["Run", "Well", "Row", "Col", "Target", "Sample", "Description", "Kind", "Cq", "Excluded"]].sort_values(
        ["Run", "Col", "Row"]).reset_index(drop=True)


# ----------------------------------------------------------------------------- Excel

def build_excel(experiments: list, wells_all: pd.DataFrame, run_labels: list, excluded: set) -> bytes:
    legend = {}
    for e in experiments:
        legend.update(e["legend"])
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        settings = pd.DataFrame([{
            "Experiment": e["name"], "Reference genes": ", ".join(e["refs"]),
            "Control (calibrator) samples": ", ".join(e["calibrators"]),
            "Efficiencies used": ", ".join(f"{g}: {(v - 1) * 100:.0f}%" for g, v in e["efficiency"].items()) or "100% for all genes",
            "Wells excluded": ", ".join(sorted(k for k in e["excluded"])) or "none",
        } for e in experiments])
        settings.to_excel(xw, sheet_name="Settings", index=False)
        pd.concat([results_table(e).assign(Experiment=e["name"]) for e in experiments], ignore_index=True)[
            ["Experiment"] + list(results_table(experiments[0]).columns)].to_excel(xw, sheet_name="Results", index=False)
        well_map(wells_all, legend, excluded).to_excel(xw, sheet_name="Well map", index=False)

        ws = xw.book.create_sheet("Plate view")
        r0 = 1
        for label in run_labels:
            sub = wells_all[wells_all["Run"] == label].set_index("Well")
            ws.cell(r0, 1, label).font = Font(bold=True, size=12, color="1F4E79")
            for c in range(1, 13):
                h = ws.cell(r0 + 1, c + 1, c)
                h.font = Font(bold=True, color="FFFFFF")
                h.fill = PatternFill("solid", start_color="1F497D")
                h.alignment = Alignment(horizontal="center")
            for i, letter in enumerate(ROWS):
                ws.cell(r0 + 2 + i, 1, letter).font = Font(bold=True)
                for c in range(1, 13):
                    well = f"{letter}{c:02d}"
                    cell = ws.cell(r0 + 2 + i, c + 1)
                    cell.alignment = Alignment(wrap_text=True, horizontal="center", vertical="center")
                    if well in sub.index:
                        row = sub.loc[well]
                        cq = "no Cq" if pd.isna(row["Cq"]) else f"{row['Cq']:.2f}"
                        cell.value = f"{row['Target']}\n{row['Sample']}\n{cq}"
                        if f"{label}|{well}" in excluded:
                            cell.fill = PatternFill("solid", start_color="E5E7EB")
                    ws.row_dimensions[r0 + 2 + i].height = 45
            r0 += 12
        for c in range(2, 14):
            ws.column_dimensions[get_column_letter(c)].width = 16

        qc_row = 0
        for e in experiments:
            for title, df in (("Reference genes", e["ref_table"]), ("Replicates with SD above the limit", e["flags"]),
                              ("Negative controls", e["ntc"])):
                pd.DataFrame({e["name"]: [title]}).to_excel(xw, sheet_name="QC", startrow=qc_row, index=False)
                df.to_excel(xw, sheet_name="QC", startrow=qc_row + 1, index=False)
                qc_row += len(df) + 4

        for sheet in xw.book.worksheets:
            if sheet.title == "Plate view":
                continue
            for col in sheet.columns:
                width = max(len(str(c.value or "")) for c in col[:60])
                sheet.column_dimensions[get_column_letter(col[0].column)].width = min(max(width + 2, 10), 60)
            sheet.freeze_panes = "A2" if sheet.title != "QC" else None
    return buf.getvalue()


# ----------------------------------------------------------------------------- figures

def _png(fig, width):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=190, bbox_inches="tight")
    plt.close(fig)
    w_px, h_px = plt.imread(io.BytesIO(buf.getvalue())).shape[1::-1]
    buf.seek(0)
    return Image(buf, width=width, height=width * h_px / w_px)


def bar_chart(pf: pd.DataFrame, target: str, height: float = 3.1):
    """Grouped bars: groups (siRNA / control) on x, one bar per dose, centred inside each group."""
    fig, ax = plt.subplots(figsize=(11, height))
    groups = list(dict.fromkeys(pf["Group"]))
    treated = [d for d in dict.fromkeys(pf["Dose"]) if d]
    width = 0.8 / max(pf.groupby("Group")["Dose"].nunique().max(), 1)
    top = float(pf["rel_hi"].max())
    seen = set()
    for g in groups:
        sub = pf[pf["Group"] == g]
        doses = list(dict.fromkeys(sub["Dose"]))
        for j, dose in enumerate(doses):
            row = sub[sub["Dose"] == dose].iloc[0]
            x = groups.index(g) + (j - (len(doses) - 1) / 2) * width
            color = CONTROL_COLOR if not dose else DOSE_COLORS[treated.index(dose) % len(DOSE_COLORS)]
            label = (dose or "control / no dose") if dose not in seen else None
            seen.add(dose)
            ax.bar(x, row["rel"], width * 0.92, color=color, label=label,
                   yerr=[[max(row["rel"] - row["rel_lo"], 0)], [max(row["rel_hi"] - row["rel"], 0)]], capsize=2.5, error_kw={"lw": 0.9})
            ax.text(x, row["rel_hi"] + top * 0.015, f"{row['rel']:.2f}", ha="center", va="bottom", fontsize=6.5, color="#374151")
    ax.set_ylim(0, top * 1.12)
    ax.axhline(1.0, color="#C00000", ls="--", lw=0.9)
    ax.set_xticks(range(len(groups)), [g if len(g) < 22 else g[:20] + "…" for g in groups], rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Relative expression\n(control = 1)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25)
    if len(seen) > 1:
        ax.legend(title="dose", fontsize=7, title_fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    return fig


def melt_figure(exp: dict, wells: pd.DataFrame, runs: dict):
    targets = exp["targets"]
    fig, axes = plt.subplots(1, len(targets), figsize=(3.9 * min(len(targets), 3), 2.6), squeeze=False)
    for ax, t in zip(axes[0], targets):
        x, xcol, curves = curves_for(wells, runs, t, "melt")
        for _, _, kind, y in curves:
            ax.plot(x, y, color="#d62728" if kind != "sample" else "#2ca02c", lw=0.8, alpha=0.8,
                    ls="--" if kind != "sample" else "-")
        if not curves:
            ax.text(0.5, 0.5, "no melt data", ha="center", va="center", transform=ax.transAxes, color="#888")
        ax.set_title(t, fontsize=9)
        ax.set_xlabel("°C", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------- PDF

def _table(rows, widths, header=True):
    t = Table(rows, colWidths=widths, repeatRows=1 if header else 0)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY), ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT]),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D1D5DB")), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]))
    return t


def _df_table(df: pd.DataFrame, widths_cm: list, fmt: dict | None = None):
    fmt = fmt or {}
    rows = [[Paragraph(escape(str(c)), HEAD) for c in df.columns]]
    for _, r in df.iterrows():
        rows.append([Paragraph(escape(_f(v, fmt.get(c, ".2f")) if isinstance(v, (float, np.floating)) else str(v)),
                               CELL if i == 0 else CELL_C) for i, (c, v) in enumerate(r.items())])
    if len(widths_cm) != len(df.columns):
        widths_cm = [1] * len(df.columns)
    scale = CONTENT_W / (sum(widths_cm) * cm)
    return _table(rows, [w * cm * scale for w in widths_cm])


def _footer(title):
    def draw(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(GREY)
        canvas.drawString(MARGIN, 0.8 * cm, f"qPCR knockdown analysis  |  {title}")
        canvas.drawRightString(PAGE[0] - MARGIN, 0.8 * cm, f"Page {doc.page}")
        canvas.restoreState()
    return draw


def build_pdf(title: str, experiments: list, wells_all: pd.DataFrame, runs: dict) -> bytes:
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=PAGE, leftMargin=MARGIN, rightMargin=MARGIN, topMargin=MARGIN,
                            bottomMargin=1.4 * cm, title=f"qPCR knockdown analysis - {title}")
    story = [Paragraph("qPCR knockdown analysis", H1),
             Paragraph(f"Run(s): <b>{escape(title)}</b> &nbsp;|&nbsp; Generated: {date.today().isoformat()}", SMALL),
             Spacer(1, 6),
             Paragraph("Method: relative expression = efficiency-corrected 2<super>-ddCt</super> (Livak / Pfaffl). Cq values are averaged over "
                       "technical replicates, normalised to the reference gene(s) and to the control sample(s). Error bars span "
                       "+/- 1 SD of the technical replicates propagated through the calculation; they do not include biological variation.",
                       BODY)]
    for e in experiments:
        story += [Paragraph(f"Experiment: {escape(' | '.join(e['targets']))}", H2),
                  Paragraph(f"<b>Reference gene(s):</b> {escape(', '.join(e['refs']) or 'none')} &nbsp;&nbsp; "
                            f"<b>Control (calibrator):</b> {escape(', '.join(e['calibrators']) or 'none')} &nbsp;&nbsp; "
                            f"<b>Wells excluded:</b> {escape(', '.join(sorted(e['excluded'])) or 'none')}", BODY)]
        for note in e["notes"]:
            story.append(Paragraph("Layout note: " + escape(note), SMALL))
        for warn in e["warnings"]:
            story.append(Paragraph("<font color='#C62828'>Warning: " + escape(warn) + "</font>", BODY))
        if len(e["ref_table"]):
            story += [Spacer(1, 4), _df_table(e["ref_table"], [3, 2, 2.4, 3.2, 2.6, 4], {"Mean Cq": ".2f", "SD between samples": ".2f", "Range (cycles)": ".2f"})]

    for e in experiments:
        res = results_table(e)
        for target in dict.fromkeys(res["Target"]) if not res.empty else []:
            pf = plot_frame(e["results"], target, e["legend"])
            story.append(PageBreak())
            story.append(Paragraph(f"{escape(target)}", H1))
            if pf.empty:
                story.append(Paragraph("No relative expression could be calculated (check reference genes and control samples).", BODY))
                continue
            story += [_png(bar_chart(pf, target), CONTENT_W), Spacer(1, 4)]
            t = res[res["Target"] == target].copy()
            t["Relative expression (range)"] = [f"{_f(a, '.2f')}  ({_f(b, '.2f')} - {_f(c, '.2f')})" for a, b, c in
                                                zip(t["Relative expression"], t["Rel. expr. low (-SD)"], t["Rel. expr. high (+SD)"])]
            t["Mean Cq +/- SD"] = [f"{_f(a, '.2f')} +/- {_f(b, '.2f')}" for a, b in zip(t["Mean Cq"], t["SD Cq"])]
            t["Replicates"] = [f"{a}/{b}" for a, b in zip(t["Replicates used"], t["Replicates on plate"])]
            t["Knockdown"] = [_f(v, ".0f") + ("%" if v == v else "") for v in t["Knockdown (%)"]]
            out = t[["Sample", "Description", "Role", "Replicates", "Mean Cq +/- SD", "ddCt", "Relative expression (range)", "Knockdown"]]
            out = out.assign(ddCt=[_f(v, ".2f") for v in t["ddCt"]])
            story.append(_df_table(out, [2.6, 5.2, 2.4, 2, 3.4, 1.6, 4.6, 2.2]))

    for e in experiments:
        story += [PageBreak(), Paragraph(f"Quality control: {escape(' | '.join(e['targets']))}", H1)]
        story.append(Paragraph("Negative controls", H2))
        story.append(_df_table(e["ntc"], [3, 2, 2, 2, 2.6, 3.4, 3], {"Lowest Cq": ".2f", "Latest sample Cq": ".2f"}) if len(e["ntc"])
                     else Paragraph("No NTC wells in this experiment.", BODY))
        story.append(Paragraph(f"Replicates with SD above the limit ({e['sd_limit']:g} cycles)", H2))
        if len(e["flags"]):
            story.append(_df_table(e["flags"].head(14), [2.6, 2.8, 1, 2, 1.8, 1.8, 7], {"Mean Cq": ".2f", "SD": ".2f", "Spread": ".2f"}))
        else:
            story.append(Paragraph("All replicate groups are within the limit.", BODY))
        story += [Paragraph("Melt curves (green: samples, red dashed: negative controls)", H2),
                  _png(melt_figure(e, e["wells"], runs), CONTENT_W * 0.9)]

    footer = _footer(title)
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    return buf.getvalue()


# ----------------------------------------------------------------------------- layout template

def template_xlsx() -> bytes:
    """An empty plate layout in the format the app reads: a gene grid, a sample grid and a legend."""
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Plate1"
    ws["A1"] = "qPCR plate 1: describe the plate here (free notes are shown in the app)"
    for top, title, example in ((4, "Primer mix", "GENE_h"), (15, "Sample", "1 Ctrl")):
        ws.cell(top + 1, 1, title).font = Font(bold=True)
        for c in range(1, 13):
            ws.cell(top + 1, c + 1, c).font = Font(bold=True)
        for i, letter in enumerate(ROWS):
            ws.cell(top + 2 + i, 1, letter).font = Font(bold=True)
        ws.cell(top + 2, 2, example)
    ws["O5"] = "Samples (legend)"
    ws["O5"].font = Font(bold=True)
    for i, (sid, desc) in enumerate([("1 Ctrl", "Control"), ("2 Scr", "Scramble"), ("3 si1-10", "siGENE-1, 10 nM"),
                                     ("4 si1-30", "siGENE-1, 30 nM"), ("NTC", "No-template control")]):
        ws.cell(6 + i, 15, sid)
        ws.cell(6 + i, 16, desc)
    ws.column_dimensions["A"].width = 12
    ws.column_dimensions["O"].width = 18
    ws.column_dimensions["P"].width = 28
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()
