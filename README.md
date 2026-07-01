# Sequence Aligner

Build and tune your own DNA and protein sequence aligner. This repo has two parts:

1. A **native Windows app** (WinUI 3, C#) where you tune the scoring and watch the alignment change live.
2. A **Python engine** with the same algorithms, plus a one-click runner that writes all results and charts to a folder.

Both run the same methods: global (Needleman-Wunsch) and local (Smith-Waterman), with a linear or affine gap model, for DNA (match / mismatch) or protein (BLOSUM62 / PAM250).

## The Windows app

Folder: `SequenceAlignerApp/`

Features:
- Two sequence boxes, or load a FASTA file, or pick a built-in sample (toy DNA, hemoglobin, SARS spike).
- Sliders for match reward, mismatch penalty, and gap penalty. The alignment redraws as you move them.
- Global or local, DNA or protein, BLOSUM62 or PAM250, linear or affine gaps.
- Colored alignment view (green match, red mismatch, gray gap), with score, percent identity, and gap counts.
- Translate DNA to protein, a dot plot, and copy or save the alignment.

Build and run (needs the .NET 9 SDK; the Windows App SDK restores automatically):

```
cd SequenceAlignerApp
dotnet build
dotnet run
```

## The Python engine

- `app/engine.py` is the alignment engine (no third-party dependencies).
- `run.py` (or double-click `run.bat` on Windows) runs every analysis and writes all results and charts into an **Output** folder.
- `app/make_figures.py` and `app/make_3d.py` draw the charts.

Setup:

```
pip install -r requirements.txt
python run.py
```

## What is in here

| Path | What it is |
|------|-----------|
| `SequenceAlignerApp/` | The WinUI 3 native app (C#) |
| `app/engine.py` | The Python alignment engine |
| `run.py`, `run.bat` | One-click runner that fills the Output folder |
| `Alignment/` | The original algorithm files and the sequence data |
| `figures/` | The charts as images |

## The science, in one line

Changing the scoring changes the alignment, so the score alone cannot tell you which alignment is right. Percent identity, sensible gaps, and biology have to.
