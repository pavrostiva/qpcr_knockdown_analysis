import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import kd_core as kd
import kd_report

st.set_page_config(page_title="qPCR knockdown analysis", page_icon="🧬", layout="wide")

DOSE_COLORS = ["#1F4E79", "#4F8FCB", "#9DC3E6", "#C9A227", "#7A5195"]
CONTROL_COLOR = "#9CA3AF"


@st.cache_data(show_spinner="Reading files...")
def load_files(files: tuple) -> kd.Ingest:
    return kd.ingest(list(files))


def bar_figure(pf: pd.DataFrame, target: str, log_scale: bool) -> go.Figure:
    groups = list(dict.fromkeys(pf["Group"]))
    doses = list(dict.fromkeys(pf["Dose"]))
    treated = [d for d in doses if d]
    fig = go.Figure()
    for dose in doses:
        sub = pf[pf["Dose"] == dose]
        color = CONTROL_COLOR if not dose else DOSE_COLORS[treated.index(dose) % len(DOSE_COLORS)]
        fig.add_trace(go.Bar(
            x=sub["Group"], y=sub["rel"], name=dose or "no dose", marker_color=color,
            error_y=dict(type="data", symmetric=False, array=(sub["rel_hi"] - sub["rel"]).clip(lower=0),
                         arrayminus=(sub["rel"] - sub["rel_lo"]).clip(lower=0)),
            customdata=np.stack([sub["Sample"].astype(str), sub["Description"], sub["Cq"], sub["knockdown_pct"], sub["n_used"]], axis=-1),
            hovertemplate="<b>%{customdata[0]}</b><br>%{customdata[1]}<br>Relative expression: %{y:.2f}<br>"
                          "Knockdown: %{customdata[3]:.0f}%<br>Mean Cq: %{customdata[2]:.2f} (n=%{customdata[4]})<extra></extra>",
        ))
    fig.add_hline(y=1.0, line_dash="dash", line_color="#C00000")
    fig.update_layout(barmode="group", height=380, margin=dict(t=30, b=10), yaxis_title="Relative expression (control = 1)",
                      xaxis=dict(categoryorder="array", categoryarray=groups), legend_title="dose" if treated else None,
                      showlegend=len(doses) > 1 or bool(doses[0]))
    if log_scale:
        fig.update_yaxes(type="log")
    return fig


def curve_figure(x, xlabel, curves, title) -> go.Figure:
    fig = go.Figure()
    for _, sample, kind, y in curves:
        ntc = kind != "sample"
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=f"{sample}", line=dict(
            width=1, color="#d62728" if ntc else "#2ca02c", dash="dash" if ntc else "solid")))
    fig.update_layout(title=title, xaxis_title=xlabel, height=300, showlegend=False, margin=dict(t=40, b=10))
    return fig


# ----------------------------------------------------------------------------- sidebar

with st.sidebar:
    st.markdown("## 🧬 qPCR knockdown analysis")
    st.caption("Relative expression (ΔΔCt) after siRNA / CRISPR knockdown.")
    st.markdown("---")
    use_demo = st.checkbox("🧪 Load demo data", value=False, help="Synthetic siRNA experiment, no files needed.")
    uploaded = []
    if not use_demo:
        uploaded = st.file_uploader(
            "Upload run files", type=["zip", "csv", "xlsx"], accept_multiple_files=True,
            help="Drop everything at once: the ZIP export(s) of each run and the Excel plate layout(s).")
        st.info("📋 **What to upload**\n- one ZIP per run (or its `... Quantification Cq Results.csv`; amplification and melt files "
                "are optional)\n- the plate layout `.xlsx` with a gene grid and a sample grid\n\nSeveral plates can be analysed together.")
        st.download_button("📄 Plate layout template", kd_report.template_xlsx(), "plate_layout_template.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch")
    st.markdown("---")
    st.caption("Files are processed in memory and are not stored. For confidential data, run the app locally (see README).")
    st.caption(f"Version **{kd.__version__}** · {kd.__date__}")

# ----------------------------------------------------------------------------- loading

if use_demo:
    ing = kd.demo_ingest()
elif uploaded:
    ing = load_files(tuple((f.name, f.getvalue()) for f in uploaded))
else:
    st.title("🧬 qPCR knockdown analysis")
    st.markdown(
        "Upload your Bio-Rad CFX export(s) and the plate layout in the sidebar. You get the **well → gene → sample table**, "
        "**relative expression** of every gene in every sample (ΔΔCt), **knockdown %**, quality control, and a **PDF / Excel report**.\n\n"
        "Not sure what to upload? Tick *Load demo data* first.")
    st.stop()

for w in ing.warnings:
    st.warning(w)
if not ing.runs:
    st.error("No `... Quantification Cq Results.csv` found. Upload the run export (ZIP or CSV).")
    st.stop()

st.title("🧬 qPCR knockdown analysis")

# ----------------------------------------------------------------------------- step 1: match layouts

st.markdown("### 1 · Runs and plate layouts")
INSTRUMENT = "Use names typed into the instrument software"
runs = {r.label: r for r in ing.runs}
chosen = {}
cols = st.columns(min(len(ing.runs), 2))
for i, run in enumerate(ing.runs):
    options = [lay.name for lay in ing.layouts]
    inst = kd.layout_from_instrument(run.cq)
    if inst is not None:
        options.append(INSTRUMENT)
    with cols[i % len(cols)]:
        if not options:
            st.error(f"**{run.label}**: no plate layout found. Upload the layout `.xlsx` too.")
            st.stop()
        default = kd.best_layout(run, ing.layouts) if ing.layouts else len(options) - 1
        pick = st.selectbox(f"Layout for run: {run.label[:70]}", options, index=default, key=f"lay_{run.label}")
        chosen[run.label] = inst if pick == INSTRUMENT else ing.layouts[options.index(pick)]

wells_all = pd.concat([kd.build_wells(runs[label], chosen[label]) for label in runs], ignore_index=True)
if wells_all.empty:
    st.error("The layout does not match any well of the run.")
    st.stop()
legend_all = {}
for lay in chosen.values():
    legend_all.update(lay.legend)

n_used = wells_all.groupby("Run").size()
st.caption("Wells used per run: " + " · ".join(f"{k[:40]}: **{v}**" for k, v in n_used.items()))

missing = kd.missing_wells(wells_all)
if len(missing):
    st.warning(f"{len(missing)} sample well(s) have no Cq (no amplification or not pipetted). They are listed in the QC section.")

with st.expander("🗺️ Well map: which gene and sample is in each well", expanded=False):
    wm = kd_report.well_map(wells_all, legend_all, set())
    st.dataframe(wm, width="stretch", hide_index=True)
    st.download_button("📥 Download well map (CSV)", wm.to_csv(index=False), "well_map.csv", mime="text/csv")

# ----------------------------------------------------------------------------- step 2: experiments

experiments = kd.find_experiments(wells_all)
st.markdown("### 2 · Analysis")
if len(experiments) > 1:
    st.caption(f"{len(experiments)} independent experiments were found (genes that share samples are analysed together).")

exp_out = []


def render_experiment(idx: int, exp: dict) -> dict:
    eid = f"e{idx}_" + "|".join(sorted(exp["targets"]))
    wells = exp["wells"]
    layouts_here = [chosen[r] for r in dict.fromkeys(wells["Run"])]
    samples = kd.order_samples(exp["samples"], layouts_here)
    legend = {s: legend_all.get(s, "") for s in samples}
    notes = list(dict.fromkeys(n for lay in layouts_here for n in lay.notes))
    warnings = []

    if notes:
        with st.expander("📝 Notes found in the plate layout file", expanded=any("normali" in n.lower() for n in notes)):
            for n in notes:
                st.write("• " + n)

    stable, unstable = kd.suggest_references(wells)
    note_refs = kd.refs_from_notes(notes, exp["targets"], samples, always=len(experiments) == 1)
    default_refs = note_refs or stable or unstable
    if note_refs:
        st.info(f"Reference gene(s) preselected from the layout note: **{', '.join(note_refs)}**.")

    c1, c2 = st.columns(2)
    refs = c1.multiselect("Reference gene(s)", exp["targets"], default=default_refs, key=f"{eid}_refs",
                          help="Housekeeping genes that do not change with the treatment. Several genes are averaged (geometric mean).")
    ctrl_default = [s for s in samples if "scr" in s.lower()] or [s for s in samples if kd.CONTROL_RE.search(s)][:1] or samples[:1]
    cals = c2.multiselect("Control sample(s) = 1.0", samples, default=ctrl_default, key=f"{eid}_cal",
                          format_func=lambda s: f"{s}  ({legend[s]})" if legend.get(s) else s,
                          help="Relative expression is shown against these samples (several are averaged), e.g. the scrambled siRNA.")

    with st.expander("⚙️ Primer efficiencies (optional) and excluded wells"):
        st.caption("Default is 100 % (ideal doubling per cycle). Enter the values from the primer validation tool if you have them.")
        ecols = st.columns(min(len(exp["targets"]), 4))
        eff = {}
        for i, g in enumerate(exp["targets"]):
            pct = ecols[i % len(ecols)].number_input(f"{g} efficiency (%)", 50.0, 110.0, 100.0, 1.0, key=f"{eid}_eff_{g}")
            if abs(pct - 100.0) > 1e-9:
                eff[g] = 1.0 + pct / 100.0
        sw = wells[wells["Kind"] == "sample"].dropna(subset=["Cq"])
        labels = {r.WellKey: f"{r.Well} · {r.Target} · {r.Sample} · Cq {r.Cq:.2f}" for r in sw.itertuples()}
        excluded = set(st.multiselect("Exclude outlier wells from the calculation", list(labels), key=f"{eid}_excl",
                                      format_func=lambda k: labels[k]))

    flags = kd.replicate_flags(wells, excluded)
    ntc = kd.ntc_table(wells, excluded)
    candidates = list(dict.fromkeys([t for t in exp["targets"] if kd.looks_like_reference(t)] + refs))
    ref_table = kd.reference_table(wells, candidates, excluded)
    if len(ref_table) and len(ntc):
        status = ntc[ntc["Control"] == "NTC"].set_index("Target")["Status"]
        ref_table["NTC"] = ref_table["Gene"].map(status).fillna("")

    if not refs:
        st.warning("Choose at least one reference gene to calculate relative expression.")
    if not cals:
        st.warning("Choose at least one control sample.")
    for g in refs:
        row = ref_table[ref_table["Gene"] == g]
        if len(row) and row.iloc[0]["Status"] != "stable":
            warnings.append(f"Reference gene {g} varies by {row.iloc[0]['Range (cycles)']:.1f} cycles between samples; "
                            "the treatment may affect it.")
        if len(row) and "NTC" in row and row.iloc[0]["NTC"] == "Check":
            warnings.append(f"Reference gene {g}: the negative control amplifies late but within the sample range; check the primers.")
    if len(refs) > 1:
        means = ref_table[ref_table["Gene"].isin(refs)]["Mean Cq"]
        if len(means) > 1 and means.max() - means.min() > 5:
            warnings.append(f"The reference genes differ by {means.max() - means.min():.0f} cycles in mean Cq. "
                            "Such different abundances make averaging questionable; check that every primer pair works.")
    ref_runs = set(wells.loc[wells["Target"].isin(refs), "Run"])
    goi_runs = set(wells.loc[~wells["Target"].isin(refs), "Run"])
    if refs and goi_runs and not goi_runs <= ref_runs:
        warnings.append("Reference and target genes were measured on different plates. Plate-to-plate variation then adds to "
                        "the error; the same cDNA must have been used.")
    for w in warnings:
        st.warning(w)

    results = kd.analyze(wells, refs, cals, eff, excluded, samples) if refs and cals else pd.DataFrame()
    if refs and cals and all(t in refs for t in exp["targets"]):
        st.info("This experiment contains only reference gene(s) (" + ", ".join(exp["targets"]) + "). Upload the run with the target "
                "gene(s) too, or deselect a reference gene, to calculate relative expression.")

    if len(results) and "rel" in results:
        log_scale = st.checkbox("Logarithmic y axis", key=f"{eid}_log")
        for g in [t for t in exp["targets"] if t not in refs]:
            pf = kd.plot_frame(results, g, legend)
            st.markdown(f"#### {g}")
            if pf.empty:
                st.info("No relative expression could be calculated for this gene (missing Cq in the gene, reference or control).")
                continue
            tested = pf[~pf["is_control"].fillna(False)]
            if len(tested):
                best = tested.loc[tested["knockdown_pct"].idxmax()]
                st.caption(f"Strongest knockdown: **{best['knockdown_pct']:.0f} %** in {best['Sample']}"
                           + (f" ({best['Description']})" if best["Description"] else "") + ". Error bars: ±1 SD of technical replicates.")
            st.plotly_chart(bar_figure(pf, g, log_scale), width="stretch", key=f"{eid}_bar_{g}")
            show = kd_report.results_table({"results": results, "legend": legend})
            show = show[show["Target"] == g].drop(columns=["Target"])
            st.dataframe(show, width="stretch", hide_index=True, column_config={
                "Mean Cq": st.column_config.NumberColumn(format="%.2f"), "SD Cq": st.column_config.NumberColumn(format="%.2f"),
                "dCt": st.column_config.NumberColumn(format="%.2f"), "ddCt": st.column_config.NumberColumn(format="%.2f"),
                "Relative expression": st.column_config.NumberColumn(format="%.3f"),
                "Rel. expr. low (-SD)": st.column_config.NumberColumn(format="%.3f"),
                "Rel. expr. high (+SD)": st.column_config.NumberColumn(format="%.3f"),
                "Knockdown (%)": st.column_config.NumberColumn(format="%.0f")})

    with st.expander("🔬 Quality control", expanded=bool(len(flags) or (len(ntc) and (ntc["Status"] == "Check").any()))):
        st.markdown("**Reference genes** (mean Cq of each sample; a stable gene moves less than 1 cycle)")
        st.dataframe(ref_table, width="stretch", hide_index=True) if len(ref_table) else st.caption("No housekeeping-like genes detected.")
        st.markdown(f"**Negative controls** (pass: no Cq, Cq ≥ 35, or {kd.NTC_MARGIN:g}+ cycles after the latest sample)")
        st.dataframe(ntc, width="stretch", hide_index=True) if len(ntc) else st.caption("No NTC wells in the layout.")
        st.markdown(f"**Replicates with SD > {kd.SD_LIMIT:g} cycles**")
        st.dataframe(flags, width="stretch", hide_index=True) if len(flags) else st.success("All replicate groups are tight.")
        miss = missing[missing["Target"].isin(exp["targets"])]
        if len(miss):
            st.markdown("**Wells without Cq**")
            st.dataframe(miss, width="stretch", hide_index=True)
        tcol, kcol = st.columns([2, 1])
        tsel = tcol.selectbox("Curves of gene", exp["targets"], key=f"{eid}_curve")
        ksel = kcol.radio("Curve type", ["Melt peaks", "Amplification"], horizontal=True, key=f"{eid}_ckind")
        x, xlabel, curves = kd.curves_for(wells, runs, tsel, "melt" if ksel == "Melt peaks" else "amp")
        if curves:
            st.plotly_chart(curve_figure(x, xlabel, curves, f"{tsel}: {ksel.lower()}"), width="stretch", key=f"{eid}_curves")
            st.caption("Green: samples, red dashed: negative controls. A single sharp melt peak means one specific product.")
        else:
            st.caption("The export has no curve file for this gene.")

    return {"name": " | ".join(exp["targets"]), "targets": exp["targets"], "samples": samples, "refs": refs, "calibrators": cals,
            "efficiency": eff, "excluded": excluded, "results": results, "ref_table": ref_table, "flags": flags, "ntc": ntc,
            "legend": legend, "notes": notes, "warnings": warnings, "wells": wells, "sd_limit": kd.SD_LIMIT}


if len(experiments) > 1:
    tabs = st.tabs([" · ".join(e["targets"]) for e in experiments])
    for i, (tab, e) in enumerate(zip(tabs, experiments)):
        with tab:
            exp_out.append(render_experiment(i, e))
else:
    exp_out.append(render_experiment(0, experiments[0]))

# ----------------------------------------------------------------------------- step 3: downloads

st.markdown("### 3 · Reports")
ready = [e for e in exp_out if len(e["results"]) and "rel" in e["results"]]
if not ready:
    st.info("Choose reference genes and control samples above to enable the reports.")
    st.stop()

title = " + ".join(dict.fromkeys(r.label for r in ing.runs))[:120]
signature = tuple((e["name"], tuple(e["refs"]), tuple(e["calibrators"]), tuple(sorted(e["efficiency"].items())),
                   tuple(sorted(e["excluded"]))) for e in exp_out)
all_excluded = set().union(*[e["excluded"] for e in exp_out])
col_pdf, col_dl, col_xls = st.columns(3)
with col_pdf:
    if st.button("📄 Build PDF report", width="stretch"):
        with st.spinner("Rendering the PDF..."):
            st.session_state["pdf"] = (signature, kd_report.build_pdf(title, ready, wells_all, runs))
state = st.session_state.get("pdf")
with col_dl:
    if state and state[0] == signature:
        st.download_button("⬇️ Download PDF", state[1], "qPCR_knockdown_report.pdf", mime="application/pdf", width="stretch")
    elif state:
        st.caption("Settings changed since the last build. Build the PDF again.")
with col_xls:
    st.download_button("📊 Download Excel report", kd_report.build_excel(ready, wells_all, list(runs), all_excluded),
                       "qPCR_knockdown_report.xlsx", width="stretch",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
st.caption("The Excel file also contains the well map (well → gene → sample → Cq) and a plate view of each run.")
