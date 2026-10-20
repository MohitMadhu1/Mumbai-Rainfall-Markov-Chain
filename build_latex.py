"""Build the Overleaf-ready LaTeX project from report_draft.md.

    .venv/Scripts/python build_latex.py   ->  latex/  and  report_overleaf.zip

The body is converted with pandoc, then wrapped in a plain pdfLaTeX preamble
(standard packages only). Unicode symbols are replaced with LaTeX commands so the
default Overleaf compiler (pdfLaTeX) needs no extra configuration.
"""

import os
import re
import shutil
import zipfile

import pypandoc

OUT = "latex"

PREAMBLE = r"""\documentclass[11pt,a4paper]{article}

% ---------------------------------------------------------------- packages
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage[a4paper,margin=2.3cm]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{longtable,booktabs,array,calc}
\usepackage{caption}
\usepackage{float}
\usepackage{enumitem}
\usepackage{xcolor}
\usepackage{textcomp}
\usepackage[colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black,citecolor=blue!50!black]{hyperref}

\setlength{\parskip}{4pt}
\setlength{\parindent}{0pt}
\setlist{itemsep=2pt,topsep=3pt}
\captionsetup{font=small,labelfont=bf}
\graphicspath{{figures/}}

% pandoc helpers
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\makeatletter
\newsavebox\pandoc@box
\newcommand*\pandocbounded[1]{%
  \sbox\pandoc@box{#1}%
  \Gscale@div\@tempa{\textheight}{\dimexpr\ht\pandoc@box+\dp\pandoc@box\relax}%
  \Gscale@div\@tempb{\linewidth}{\wd\pandoc@box}%
  \ifdim\@tempb\p@<\@tempa\p@\let\@tempa\@tempb\fi
  \ifdim\@tempa\p@<\p@\scalebox{\@tempa}{\usebox\pandoc@box}%
  \else\usebox{\pandoc@box}\fi}
\makeatother

% ---------------------------------------------------------------- title
\title{\textbf{A Markov Chain Model of Monsoon Rainfall in Mumbai}\\[4pt]
  \large Data-Driven Stochastic Modeling and Simulation of Real-World Systems}
\author{Member A \and Member B \and Member C}  % <-- put your names here
\date{2026}

\begin{document}
\maketitle
"""

ENDING = r"""
\appendix
\section{Code}

All code is our own, written in Python with \texttt{numpy} and \texttt{scipy} only; no Markov-chain
or statistics package does the modelling for us. The complete code is available at:

\begin{center}
% >>> PASTE YOUR GOOGLE DRIVE LINK BETWEEN THE BRACES BELOW <<<
\url{https://drive.google.com/PASTE-YOUR-LINK-HERE}
\end{center}

The folder contains:
\begin{itemize}
  \item \texttt{analysis.ipynb} -- the Jupyter notebook that runs the whole study from top to bottom
        (data, estimation, theory, simulation, validation, upgrades, Poisson analysis and robustness
        checks) and produces every table and figure in this report;
  \item \texttt{markov.py} -- the core functions: maximum-likelihood estimation, stationary
        distribution, first-passage times, the season simulator, all hypothesis tests and the exact
        flood-run recursion;
  \item \texttt{viz.py} -- the shared plot style used by every figure;
  \item \texttt{data/mumbai\_daily\_rain.csv} -- the cleaned daily rainfall series (IMD grid cell and
        Santacruz gauge, 1985--2024);
  \item \texttt{download\_data.py} and \texttt{prepare\_data.py} -- scripts that rebuild the data
        from the original IMD and NOAA sources.
\end{itemize}
The random seed is fixed, so every number in this report is reproduced exactly by running the notebook.

\end{document}
"""

# Unicode -> LaTeX (text mode). Order matters: longer sequences first.
SUPERS = {"⁻": "-", "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
          "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}
TEXT_MAP = [
    ("χ²", r"$\chi^2$"), ("χ", r"$\chi$"), ("π", r"$\pi$"), ("λ", r"$\lambda$"),
    ("≥", r"$\geq$"), ("≤", r"$\leq$"), ("≈", r"$\approx$"), ("≫", r"$\gg$"),
    ("×", r"$\times$"), ("±", r"$\pm$"), ("→", r"$\rightarrow$"), ("∝", r"$\propto$"),
    ("°", r"\textdegree{}"), ("·", r"$\cdot$"), ("−", r"$-$"),
]
MATH_MAP = [("≥", r"\geq "), ("≤", r"\leq "), ("≈", r"\approx "), ("×", r"\times "),
            ("π", r"\pi "), ("λ", r"\lambda "), ("χ", r"\chi "), ("−", "-")]


def superscripts(t):
    # runs like 10⁻¹⁰ -> 10$^{-10}$
    return re.sub("[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+",
                  lambda m: "$^{" + "".join(SUPERS[c] for c in m.group(0)) + "}$", t)


def fix_unicode(tex):
    # split on inline/display math so math and text get different treatment
    parts = re.split(r"(\\\(.*?\\\)|\\\[.*?\\\])", tex, flags=re.S)
    out = []
    for i, part in enumerate(parts):
        if i % 2:  # math
            for a, b in MATH_MAP:
                part = part.replace(a, b)
        else:
            part = superscripts(part)
            for a, b in TEXT_MAP:
                part = part.replace(a, b)
        out.append(part)
    tex = "".join(out)
    tex = tex.replace("$$", "")           # merge adjacent inline math like $\approx$$^{-10}$
    return tex


def figures_as_floats(tex):
    # pandoc emits figure environments already; pin them near their text
    tex = tex.replace(r"\begin{figure}", r"\begin{figure}[H]")
    # the alt= key needs a 2023+ graphicx; drop it so older TeX installs compile too
    return re.sub(r",alt=\{[^{}]*\}", "", tex)


def main():
    s = open("report_draft.md", encoding="utf-8").read()
    body = s[: s.index("# Appendix: Code")]
    body = re.sub(r"^---.*?---\n", "", body, flags=re.S)       # title block -> \maketitle
    body = re.sub(r"!\[Figure \d+\. ", "![", body)              # LaTeX numbers figures itself
    body = re.sub(r"^# Abstract\n", "", body, flags=re.M)       # handled as abstract env
    abstract, rest = body.split("\n# 1. Introduction", 1)
    rest = "# 1. Introduction" + rest
    rest = re.sub(r"^(#+) \d+(\.\d+)*\.? ", r"\1 ", rest, flags=re.M)   # LaTeX numbers sections

    conv = lambda md: pypandoc.convert_text(md, "latex", format="markdown",
                                            extra_args=["--wrap=preserve"])
    tex_abs = fix_unicode(conv(abstract.strip()))
    tex_rest = figures_as_floats(fix_unicode(conv(rest)))
    tex_rest = tex_rest.replace(r"\section{References}\label{references}",
                                r"\section*{References}\label{references}"
                                r"\addcontentsline{toc}{section}{References}")

    left = sorted({c for c in tex_abs + tex_rest if ord(c) > 127 and c not in "–—‘’“”éñ"})
    if left:
        raise SystemExit(f"unhandled unicode: {left}")

    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    os.makedirs(f"{OUT}/figures")
    with open(f"{OUT}/main.tex", "w", encoding="utf-8") as f:
        f.write(PREAMBLE)
        f.write("\\begin{abstract}\n" + tex_abs + "\n\\end{abstract}\n\n")
        f.write(tex_rest)
        f.write(ENDING)
    used = set(re.findall(r"figures/([\w\-]+\.png)", s))
    for fn in used:
        shutil.copy(f"figures/{fn}", f"{OUT}/figures/{fn}")

    with zipfile.ZipFile("report_overleaf.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(OUT):
            for fn in files:
                full = os.path.join(root, fn)
                z.write(full, os.path.relpath(full, OUT))
    print(f"main.tex + {len(used)} figures -> report_overleaf.zip")


if __name__ == "__main__":
    main()
