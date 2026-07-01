# What to Submit and How to Run a Live Demo

This is the practical guide. It tells you exactly which files to hand in, and gives you two rehearsed ways to demo the project live (the Colab notebook and the native Windows app).

## Part A: What to Submit

Submit these four things. They are all you need for a complete grade.

1. **The Colab notebook link (`SequenceAligner.ipynb`, sharing turned ON).**
   - This is the main deliverable. It contains all 16 parts, from the alignment engine (Part 1) through the interactive app (Part 3), the good-vs-bad parameter arguments (Parts 5, 6, 8), the DNA-vs-protein comparison (Part 7), and the stretch and extra parts (affine gaps, the 3D viewer, the family trees, the whole-genome alignment).
   - **Turn sharing on before you submit.** In Colab: click **Share** (top right), set access to **"Anyone with the link"** with the **Viewer** role, then copy the link. If sharing is off, the grader opens the link and sees nothing but a permission wall.
   - Paste that link into your submission text.

2. **The slides: `SequenceAligner_Slides.pptx`.**
   - This is the presentation deck. It carries the "good alignment vs bad alignment" story visually, which is the central argument of the project.
   - Submit the `.pptx` file itself (or a Google Slides / PDF export if the assignment asks for a link instead of a file).

3. **The write-up, embedded in the notebook.**
   - You are not submitting a separate essay document. The write-up lives **inside** `SequenceAligner.ipynb` as markdown text cells: the title and organization cell at the very top, the explanation cells before every part, and the **Wrap-up** cell at the very bottom.
   - The one-line thesis in the Wrap-up is the thing to be able to say out loud: changing the scoring changes the alignment, so the raw score alone cannot tell you which alignment is right. Percent identity, sensible gaps, conserved regions, and biology have to.
   - Nothing extra to attach here. Just make sure those markdown cells are present when you share the notebook.

4. **The DATA: the `Alignment/Data` folder.**
   - This is the real biological sequence data the notebook aligns. It has three datasets:
     - **`Data/Hemoglobin`** (hemoglobin sequences, including `Homo_sapiens_hemoglobin.fasta`, `Gorilla_gorilla_hemoglobin.fasta`, `Bos_taurus_hemoglobin.fasta`, `Danio_rerio_hemoglobin.fasta`, and the combined `hemoglobin_protein.fasta`).
     - **`Data/CytochromeC`** (cytochrome c across seven species, e.g. `Homo_sapiens_cytc.fasta`, `Mus_musculus_cytc.fasta`, `Saccharomyces_cerevisiae_cytc.fasta`).
     - **`Data/Coronaviruses`** (the SARS sequences: the whole `SARS-CoV-2_genome.fasta` and `SARS-CoV_genome.fasta`, the spike gene, and the spike protein `.fasta` files).
   - **Important: the notebook does not need you to upload any of this.** The sequences used by the notebook are also baked directly into the code cells (see Part 2, "The sequences"), so the notebook runs top to bottom with **no file uploads**. You submit the `Data` folder so the grader can see the original source data, not because the notebook depends on it at run time.

### Optional extras (mention these, do not stress about them)

Two GitHub repositories exist and are nice to link, but they are extras on top of the four required items above.

- **Public app repo:** `https://github.com/TomAs-1226/sequence-aligner`. This hosts the native Windows app and, on its **Releases** page, a downloadable build (`SequenceAlignerApp-win-x64.zip`). Link this if you want the grader to be able to download and run the desktop app themselves.
- **Private group repo:** `https://github.com/TomAs-1226/group-a1-alignment-algorithm-application`. This is the group's working repo (notebook, slides, data, and Python engine together). It is private, so only share it if the grader has, or can be given, access.

## Part B: Live Demo of the Notebook

This is the demo you will actually give. It needs only a browser and internet.

1. **Open the notebook in Google Colab.** Use the shared link, or go to `colab.research.google.com`, choose **File then Open notebook then GitHub / Upload**, and open `SequenceAligner.ipynb`.
2. **Run everything.** In the top menu click **Runtime**, then **Run all**. Approve any "run anyway" prompt that appears for a notebook from outside Google.
3. **Wait for it to finish.** Most parts run in a few seconds. Two parts are slow on purpose: Part 15 aligns two whole 30,000-letter SARS genomes and takes roughly 20 to 30 seconds, and Part 11 (the 3D structure) and Part 16 (the GPU batch) may take a little longer. Let the run complete before you start clicking.
4. **Scroll to Part 3, "The interactive app," and move the sliders.** This is the centerpiece. Drag the **match**, **mismatch**, and **gap** sliders and watch the colored alignment redraw live (green = match, orange = mismatch, gray = gap). Flip the **mode** dropdown between **Global** and **Local**, and the **type** dropdown between **DNA** and **Protein** (Protein lets you pick **BLOSUM62** or **PAM250**). The point to narrate: the letters never changed, only the scoring did, yet the alignment changed.
5. **Show the good-vs-bad story on the slides and in the notebook.** Point to Part 5 (the toy DNA pair) and Part 6 (the real hemoglobin beta pair). Same two sequences, two parameter sets, two very different answers, and the higher score is the worse biology. This is where the slides and the notebook say the same thing.
6. **Show the DNA-vs-protein part (Part 7, the SARS spike).** Explain the payoff: synonymous (silent) mutations show up when you align the DNA but vanish when you align the protein, because they do not change the amino acid.
7. **Show the 3D part (Part 11).** Scroll to the spike structure colored by what changed between the two sequences. In Colab you can rotate the interactive version with the mouse. This is the visual "wow" moment to end on.

**Suggested 7-minute running order:**

| Time | What you do |
|------|-------------|
| 0:00 to 1:00 | Open in Colab, hit **Runtime then Run all**, say the one-line thesis while it runs. |
| 1:00 to 3:00 | Part 3: drag the sliders live; switch Global/Local and DNA/Protein. |
| 3:00 to 4:30 | Parts 5 and 6 (plus slides): the good-vs-bad parameter story, higher score = worse biology. |
| 4:30 to 6:00 | Part 7: DNA vs protein on the SARS spike (silent mutations). |
| 6:00 to 7:00 | Part 11: the 3D spike, rotate it, close with the Wrap-up thesis. |

## Part C: Live Demo of the App

Use this if you want to show the native Windows desktop app instead of (or after) the notebook. Windows only.

1. **Download the release zip.** Go to the public repo's Releases page (`https://github.com/TomAs-1226/sequence-aligner/releases`) and download **`SequenceAlignerApp-win-x64.zip`**.
2. **Unzip it** to a normal folder (for example your Desktop). It is self-contained, so you do **not** need to install .NET or anything else first.
3. **Run `SequenceAlignerApp.exe`** from the unzipped folder. If Windows SmartScreen warns about an unknown publisher, choose **More info then Run anyway**.
4. **Pick a built-in sample.** Use the sample selector to load a ready-made pair (toy DNA, hemoglobin, cytochrome c, SARS spike, or whole SARS genomes) so you are not typing sequences by hand.
5. **Move the sliders.** Drag **match reward**, **mismatch penalty**, and **gap penalty**, and watch the colored alignment, the score, percent identity, and gap counts update live. Switch between global, local, and semi-global, and between DNA and protein (BLOSUM62 or PAM250, linear or affine gaps).
6. **Open the 3D and AI panels.** Use the **3D viewer** to predict the protein's structure (ESMFold) or load a PDB, colored by where the two sequences differ. Use the **AI** panel for the secondary-structure prediction (helix / sheet / coil) and the protein-language-model similarity readout. These panels need an internet connection.

## Part D: Troubleshooting (quick list)

- **The sliders do nothing / the app section is blank in Colab.** You did not run the notebook first. Do **Runtime then Run all** and wait for it to finish, then the widgets in Part 3 become live.
- **The 3D viewer or the AI features do not work.** They need an **internet connection**. The 3D structure prediction uses the online ESMFold service, and the AI panels download or query models. On a locked-down school network these may be blocked, so test them on open internet beforehand.
- **The app will not open the 3D viewer ("WebView2 runtime may be missing").** The desktop app renders 3D through Microsoft's **WebView2 runtime**. It ships with recent Windows, but if it is missing, install the free "Microsoft Edge WebView2 Runtime" and restart the app.

**One-line backup plan:** if the app or the internet features fail during your live demo, fall back to the Colab notebook. It contains the same algorithms and shows the same story, and its core parts run without any of the online AI services.