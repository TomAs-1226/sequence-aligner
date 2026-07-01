# Sequence Aligner

Build and tune your own DNA and protein sequence aligner. This repo has two parts:

1. A **native Windows app** (WinUI 3, C#) where you tune the scoring and watch the alignment change live.
2. A **Python engine** with the same algorithms, plus a one-click runner that writes all results and charts to a folder.

Both run the same methods: global (Needleman-Wunsch) and local (Smith-Waterman), with a linear or affine gap model, for DNA (match / mismatch) or protein (BLOSUM62 / PAM250).

## Download and run (no build needed)

Get the latest Windows build from the [Releases page](https://github.com/TomAs-1226/sequence-aligner/releases): download `SequenceAlignerApp-win-x64.zip`, unzip it, and run `SequenceAlignerApp.exe`. It is self-contained, so you do not need to install .NET or anything else. (The 3D viewer and the AI structure feature need an internet connection.)

## The Windows app

Folder: `SequenceAlignerApp/`

Features:
- Two sequence boxes, load a FASTA file, or pick a built-in sample (toy DNA, hemoglobin, cytochrome c, SARS spike, or whole SARS genomes).
- Sliders for match reward, mismatch penalty, and gap penalty. The alignment redraws as you move them.
- Global, local, or semi-global (free end gaps); DNA or protein; BLOSUM62 or PAM250; linear or affine gaps.
- Handles two whole 30,000-letter genomes with a banded aligner.
- Colored alignment view with a match line, plus score, percent identity, gap counts, and sequence stats.
- Translate DNA to protein, reverse complement, a dot plot, and copy or save the alignment.
- A 3D viewer that predicts your protein's structure (ESMFold) or loads a PDB, and colors it by where the two sequences differ.
- An AI secondary-structure predictor (helix / sheet / coil) that runs on the GPU or NPU with Windows ML.
- An AI protein-analysis panel driven by a transformer language model trained from scratch: it scores how related two proteins are, predicts a protein's family, and finds similar proteins in a built-in database, also on the GPU or NPU.

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
| `TECHNICAL_DOCUMENTATION.md` | Full write-up of every algorithm and model (also as `.docx`) |
| `SUBMISSION_AND_DEMO_GUIDE.md` | What to submit and how to run a live demo (also as `.docx`) |

## Documentation

- **[TECHNICAL_DOCUMENTATION.md](TECHNICAL_DOCUMENTATION.md)** explains every algorithm and model in detail: the five aligners (global, local, semi-global, affine, banded), the supporting algorithms and scoring matrices, the GPU batched aligner, the two on-device AI models, and the app. A Word version (`TECHNICAL_DOCUMENTATION.docx`) is included for reading offline.
- **[SUBMISSION_AND_DEMO_GUIDE.md](SUBMISSION_AND_DEMO_GUIDE.md)** lists exactly what to hand in and gives a rehearsed, timed live-demo script for both the notebook and the app.
- A shorter plain-English overview is in [`Alignment/DOCUMENTATION.md`](Alignment/DOCUMENTATION.md).

## The science, in one line

Changing the scoring changes the alignment, so the score alone cannot tell you which alignment is right. Percent identity, sensible gaps, and biology have to.
