# Mumbai monsoon Markov chain — stochastic modeling project

## What is here

| File | What it is |
|------|-----------|
| `PROJECT_PLAN.md` | Every decision we made, the team split and the timeline |
| `analysis.ipynb` | **The project.** Runs top to bottom: data → model → theory → simulation → validation → upgrades → Poisson. Use it for the live demo |
| `analysis.html` | Same notebook as a web page (opens without Python) |
| `markov.py` | All the maths, our own code (numpy/scipy only) |
| `report_overleaf.zip` | **Final report as a LaTeX project** — upload to Overleaf (New Project → Upload Project). Compiles with the default pdfLaTeX |
| `report_draft.md` / `.docx` / `.pdf` | The same report as Markdown (the source), Word and PDF |
| `build_latex.py` | Regenerates `latex/` and the zip from `report_draft.md` |
| `SLIDES.md` | 15-minute talk plan, slide by slide, with likely questions |
| `figures/` | All plots, 200 dpi |
| `data/mumbai_daily_rain.csv` | Clean daily rainfall 1985–2024 (IMD + gauge) |
| `download_data.py`, `prepare_data.py` | Rebuild the data from scratch |
| `build_notebook.py` | Regenerates `analysis.ipynb` (edit this, not the notebook, if you change the analysis) |

## Data sources

- **India Meteorological Department**, 0.25° gridded daily rainfall (Pai et al., 2014), from imdpune.gov.in via `imdlib`. Grid cell 19.00 N, 72.75 E.
- **NOAA GHCN-Daily**, station IN012070800 Bombay/Santacruz (cross-check only).

`data/raw/` holds the 40 yearly IMD files (~1 GB). Delete them to save space once
`mumbai_daily_rain.csv` exists; `download_data.py` re-fetches them.

## Run it

```bash
python -m venv .venv
.venv/Scripts/python -m pip install numpy pandas scipy matplotlib jupyter nbconvert nbformat imdlib networkx
.venv/Scripts/python download_data.py
.venv/Scripts/python prepare_data.py
.venv/Scripts/python build_notebook.py
.venv/Scripts/jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
```

The `.venv` folder here is already set up, so for a demo just open the notebook:

```bash
.venv/Scripts/jupyter notebook analysis.ipynb
```

The random seed is fixed (2026), so every number is reproducible.
