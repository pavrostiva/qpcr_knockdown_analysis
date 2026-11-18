<div align="center">

# qPCR Knockdown Analysis

**Did the knockdown work? Upload the run, get the answer.**
Relative expression (ΔΔCt), knockdown %, quality control and a ready PDF / Excel report for Bio-Rad CFX runs.

<br>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit](https://img.shields.io/badge/UI-Streamlit-FF4B4B.svg)](https://streamlit.io)

</div>

---

## What is this for?

You knocked a gene down with siRNA (or CRISPR, shRNA, ASO...) and measured it by qPCR. Now you want to know:

- How much is the gene reduced compared with the control? (**knockdown %**)
- Which siRNA and which dose works best?
- Can I trust the numbers? (reference genes, replicates, negative controls, melt curves)

The Bio-Rad software gives you raw Cq values. This app turns them into **relative expression** with error bars, shows the plate **well by well** (gene and sample in every well), and flags anything suspicious.

It is the sister tool of [qpcr_primer_validation](https://github.com/pavrostiva/qpcr_primer_validation), which checks the primers *before* the experiment.

---

## How to use it

1. **Export the run** from the CFX software as CSV (the `... Quantification Cq Results.csv` file is required; amplification and melt files are optional but recommended).
2. **Describe the plate** in an Excel file (see [Plate layout](#plate-layout)). The app offers an empty template.
3. **Upload everything at once**: the ZIP(s) or CSV(s) of the run(s) and the layout `.xlsx`. Several plates can be analysed together.
4. The app pairs every run with its layout (by date and plate number; you can change the choice) and shows:

| Section | What you get |
|---|---|
| **Well map** | Which gene and sample sits in each well, with the Cq. Downloadable as a table |
| **Analysis** | Relative expression of every gene in every sample, grouped bars by siRNA and dose, knockdown % |
| **Quality control** | Reference-gene stability, negative controls, noisy replicates (with a suggestion which well to drop), melt and amplification curves |
| **Reports** | PDF (one page per gene) and Excel (results, well map, plate view, QC) |

> **Try it first:** tick *Load demo data* in the sidebar.

---

## Plate layout

One Excel sheet per plate with two 8×12 grids: **genes** (primer mixes) on top and **samples** below. Rows are labelled **A–H** in the first column, columns **1–12**:

| Primer mix | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| **A** | `TARGET_h` | `TARGET_h` | `REF_h` | `REF_h` |
| **B** | `TARGET_h` | `TARGET_h` | `REF_h` | `REF_h` |

| Sample | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| **A** | `1 Ctrl` | `3 Scr1` | `1 Ctrl` | `3 Scr1` |
| **B** | `1 Ctrl` | `4 Scr2` | `1 Ctrl` | `4 Scr2` |

- Leave unused wells empty. Replicates can be anywhere.
- **`NTC`** marks a no-template control, **`-RT`** a no-reverse-transcriptase control.
- **Legend (optional, recommended).** Put the sample names in a column and a description next to them, e.g. `5 si1-10` → `siGENE-1, 10 nM`. The app uses it to group bars by siRNA and dose (a description of the form `name, 10 nM` is split into group and dose) and to order the samples.
- **Free-text notes** in the sheet are shown in the app. A note like *"GAPDH is knocked down in sample B2, so normalise to Actb only"* preselects ACTB as the reference for the experiment whose samples it mentions.
- If you typed `Target` and `Sample` names into the instrument software, you can use those instead of a layout file.

Genes and samples that share each other are analysed together even if they are on different plates, for example the reference gene on one plate and the target on another. Independent experiments (for instance two cell lines) appear as separate tabs.

---

## How it calculates

For each gene and sample the technical replicates are averaged. Then, relative to the **control sample(s)** you choose (for example the scrambled siRNA; several controls are averaged):

```
quantity of gene in sample  = E ^ (Cq in control - Cq in sample)        E = 2 for 100 % efficiency
reference factor            = geometric mean of that quantity over the reference gene(s)
relative expression         = quantity of target gene / reference factor
knockdown %                 = (1 - relative expression) x 100
```

With 100 % efficiency this is the classic **2^-ΔΔCt** (Livak) method. You can enter a measured efficiency per gene (for example from the primer validation tool), which gives the Pfaffl-type correction.

**Error bars** are ±1 SD of the technical replicates, propagated through the calculation (target, reference genes and control). They do **not** include biological variability. With one cDNA per condition no statistical test is possible, so the app does not show p-values: treat the result as a technical readout, and repeat with independent biological replicates before drawing conclusions.

### Checks the app does for you

| Check | Rule |
|---|---|
| **Reference gene stability** | Mean Cq of a reference gene should move less than 1 cycle between samples; otherwise the treatment may affect it |
| **Reference genes very different in abundance** | Warning when two reference genes are more than 5 cycles apart |
| **Replicates** | SD above 0.3 cycles is flagged; if dropping one well fixes it, the app tells you which |
| **Negative controls** | Pass: no Cq, Cq ≥ 35, or 5+ cycles later than the latest sample |
| **Different plates** | Warning when reference and target gene were not on the same plate |
| **Melt curves** | Shown per gene; one sharp peak means one specific product |

---

## Run it on your own computer

If your data is confidential, run the app locally: then nothing leaves your machine.

```bash
git clone https://github.com/pavrostiva/qpcr_knockdown_analysis.git
cd qpcr_knockdown_analysis
pip install -r requirements.txt
streamlit run app.py
```

The app opens at `http://localhost:8501`.

### Put it online (Streamlit Community Cloud)

1. Push the repository to GitHub.
2. Go to [share.streamlit.io](https://share.streamlit.io), choose **Create app**, pick the repository, branch `main` and the file `app.py`.
3. Deploy. You get a public link you can share with the lab.

On the hosted version uploaded files are processed in memory on Streamlit's servers. Do not put patient-identifying data there.

---

## Limitations

- Bio-Rad CFX exports (`Quantification Cq Results`). Other instruments need a different file parser.
- One fluorophore per run (SYBR Green style).
- Technical replicates only; no statistics between biological replicates.
- The effect on a reference gene can only be detected through its Cq shift; always choose reference genes that your treatment does not touch.

## Development

```bash
pip install -r requirements.txt pytest
pytest
```

The calculations live in `kd_core.py` (no Streamlit dependency), the reports in `kd_report.py`, the interface in `app.py`.

## License

MIT
