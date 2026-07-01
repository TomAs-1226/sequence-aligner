// Slide deck for the sequence-aligner project. Plain language, no em dashes.
const pptxgen = require("pptxgenjs");
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.3 x 7.5 in
pres.author = "Thomas";
pres.title = "Build and tune your own sequence aligner";

const W = 13.3, H = 7.5;
const DARK = "12332E", INK = "1A2B28", GREEN = "2E7D46", RED = "C43D2E",
      AMBER = "D98A26", MUTED = "5B6B66", CARD = "F2F5F3", LINEC = "D9E2DD";

const shadow = () => ({ type: "outer", color: "000000", blur: 7, offset: 3, angle: 90, opacity: 0.10 });

function header(slide, kicker, title, dotColor) {
  slide.addShape(pres.shapes.OVAL, { x: 0.6, y: 0.62, w: 0.16, h: 0.16, fill: { color: dotColor || GREEN } });
  slide.addText(kicker.toUpperCase(), { x: 0.85, y: 0.52, w: 11, h: 0.35, margin: 0,
    fontFace: "Arial", fontSize: 12, color: MUTED, charSpacing: 2 });
  slide.addText(title, { x: 0.58, y: 0.9, w: 12.1, h: 0.9, margin: 0,
    fontFace: "Cambria", fontSize: 30, bold: true, color: INK });
}

function contentSlide() {
  const s = pres.addSlide();
  s.background = { color: "FFFFFF" };
  return s;
}

function card(slide, x, y, w, h, fill) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.08,
    fill: { color: fill || CARD }, line: { color: LINEC, width: 1 }, shadow: shadow() });
}

// ---------- Slide 1: title ----------
(() => {
  const s = pres.addSlide();
  s.background = { color: DARK };
  s.addText("Build and tune your own sequence aligner", {
    x: 0.9, y: 2.0, w: 11.5, h: 1.6, fontFace: "Cambria", fontSize: 44, bold: true, color: "FFFFFF" });
  s.addText("Turning global and local alignment into a tool, then using it to ask which alignment is actually right.", {
    x: 0.95, y: 3.5, w: 10.6, h: 0.9, fontFace: "Arial", fontSize: 18, color: "CADCD6" });
  // alignment strip motif
  const cols = ["G","G"],pairs = [
    ["A","A"],["C","C"],["G","G"],["T","-"],["A","A"],["C","G"],["G","G"],["T","T"],
    ["A","A"],["-","C"],["G","G"],["C","C"],["A","A"],["T","T"],["G","G"],["A","A"]];
  const x0 = 0.95, y0 = 5.1, cw = 0.34;
  pairs.forEach((p, i) => {
    const [a, b] = p;
    const col = (a === "-" || b === "-") ? "6E8079" : (a === b ? GREEN : AMBER);
    [0,1].forEach(r => {
      s.addShape(pres.shapes.RECTANGLE, { x: x0 + i*cw, y: y0 + r*0.34, w: cw-0.04, h: 0.30, fill: { color: col } });
      s.addText(p[r], { x: x0 + i*cw, y: y0 + r*0.34, w: cw-0.04, h: 0.30, margin: 0,
        align: "center", valign: "middle", fontFace: "Courier New", fontSize: 11, color: "FFFFFF" });
    });
  });
  s.addText("Thomas    |    CMU Pre-College Computational Biology    |    Module 2", {
    x: 0.95, y: 6.7, w: 11, h: 0.4, fontFace: "Arial", fontSize: 13, color: "9DB3AC" });
  s.addNotes("One line opener: the algorithm is only half the story. The scoring you choose decides the answer, and this project is about how to tell a good answer from a bad one.");
})();

// ---------- Slide 2: the question ----------
(() => {
  const s = contentSlide();
  header(s, "The problem", "Change the scoring, change the alignment", RED);
  s.addText([
    { text: "The same two sequences can line up sensibly or turn into nonsense, depending only on three numbers you pick:", options: { breakLine: true, fontSize: 18, color: INK, paraSpaceAfter: 10 } },
  ], { x: 0.6, y: 2.0, w: 7.0, h: 1.4, fontFace: "Arial", valign: "top" });
  s.addText([
    { text: "So which alignment is right, and how would you know?", options: { italic: true, fontSize: 20, color: GREEN } },
  ], { x: 0.6, y: 3.6, w: 7.0, h: 1.0, fontFace: "Cambria" });

  card(s, 8.1, 1.9, 4.6, 4.2);
  s.addText("The three knobs", { x: 8.4, y: 2.1, w: 4.0, h: 0.5, fontFace: "Cambria", fontSize: 18, bold: true, color: INK });
  const knobs = [
    ["Match reward", "points for two letters that agree", GREEN],
    ["Mismatch penalty", "points lost when two letters differ", AMBER],
    ["Gap penalty", "points lost for every dash inserted", RED],
  ];
  knobs.forEach((k, i) => {
    const y = 2.75 + i*1.05;
    s.addShape(pres.shapes.OVAL, { x: 8.4, y: y+0.05, w: 0.22, h: 0.22, fill: { color: k[2] } });
    s.addText(k[0], { x: 8.75, y: y-0.05, w: 3.7, h: 0.35, margin: 0, fontFace: "Arial", fontSize: 15, bold: true, color: INK });
    s.addText(k[1], { x: 8.75, y: y+0.3, w: 3.7, h: 0.5, margin: 0, fontFace: "Arial", fontSize: 12.5, color: MUTED });
  });
  s.addNotes("Set up the whole talk: there is no single correct alignment, only the best one for the scoring you chose.");
})();

// ---------- Slide 3: the app ----------
(() => {
  const s = contentSlide();
  header(s, "What I built", "A tunable aligner that runs my own code", GREEN);
  const feats = [
    "Sliders for match, mismatch, and gap",
    "Switch between global and local alignment",
    "Works on DNA and on protein",
    "Protein scoring matrices: BLOSUM62 and PAM250",
    "Reads out the alignment, score, percent identity, and gaps",
    "Redraws live as you move the sliders",
  ];
  s.addText(feats.map((f, i) => ({ text: f, options: { bullet: { code: "2022" }, breakLine: true, paraSpaceAfter: 9 } })),
    { x: 0.6, y: 2.0, w: 6.7, h: 4.4, fontFace: "Arial", fontSize: 16, color: INK, valign: "top" });
  card(s, 7.7, 1.95, 5.0, 4.5, "FFFFFF");
  s.addImage({ path: "blosum62_heatmap.png", x: 7.85, y: 2.1, sizing: { type: "contain", w: 4.7, h: 3.9 } });
  s.addText("BLOSUM62, the protein scoring matrix the app uses", { x: 7.8, y: 6.0, w: 4.8, h: 0.4, margin: 0, align: "center", fontFace: "Arial", fontSize: 11.5, color: MUTED });
  s.addNotes("Demo moment. Show the sliders moving and the alignment recoloring. Mention it runs the global and local code I wrote, generalized to proteins.");
})();

// ---------- Slide 3b: how it works (DP grid) ----------
(() => {
  const s = contentSlide();
  header(s, "How it works", "The aligner fills a grid, then walks back", GREEN);
  s.addImage({ path: "figures/fig7_dp_grid.png", x: 0.5, y: 2.0, w: 3.9, h: 3.9 / 0.957 });
  s.addImage({ path: "figures/fig9_dotplot.png", x: 4.55, y: 2.0, w: 3.9, h: 3.9 });
  card(s, 8.7, 2.0, 4.1, 4.07, "F2F5F3");
  s.addText([
    { text: "Left: the scoring grid. ", options: { bold: true, color: INK } },
    { text: "Each box holds the best score for lining up everything up to that point. The red path back through it is the alignment.", options: { breakLine: true, paraSpaceAfter: 12, color: INK } },
    { text: "Right: a dot plot. ", options: { bold: true, color: INK } },
    { text: "A dot marks every spot where the two sequences share a short piece. Related sequences make a diagonal line, and breaks in it are insertions or deletions.", options: { color: INK } },
  ], { x: 8.95, y: 2.25, w: 3.6, h: 3.6, fontFace: "Arial", fontSize: 13.5, valign: "top" });
  s.addNotes("Two charts: how the grid is filled and walked back, and how a dot plot shows shared regions at a glance.");
})();

// ---------- Slide 4: good vs bad toy ----------
(() => {
  const s = contentSlide();
  header(s, "Good vs bad parameters", "A toy example makes it obvious", RED);
  s.addImage({ path: "figures/fig1_toy_good_bad.png", x: 0.6, y: 1.85, sizing: { type: "contain", w: 12.1, h: 3.7 } });
  card(s, 0.6, 5.9, 12.1, 1.15, "FBEDEB");
  s.addText([
    { text: "The bad settings score almost 6 times higher and reach 100 percent identity, but they are wrong. ", options: { color: INK } },
    { text: "Only biology tells them apart: one clean insertion is believable, five scattered gaps are not.", options: { color: RED, bold: true } },
  ], { x: 0.85, y: 6.05, w: 11.6, h: 0.85, fontFace: "Arial", fontSize: 15, valign: "middle" });
  s.addNotes("This is the core good-versus-bad story. Point at the scattered gaps in the bad panel.");
})();

// ---------- Slide 5: good vs bad HBB ----------
(() => {
  const s = contentSlide();
  header(s, "Good vs bad parameters", "The same trap on a real gene (hemoglobin beta)", RED);
  const iw = 8.6, ih = iw / 2.684; // fig2 true aspect, so no empty letterbox frame
  s.addImage({ path: "figures/fig2_hbb_good_bad.png", x: 0.6, y: 1.95, w: iw, h: ih });
  card(s, 9.4, 1.95, 3.3, ih, "F2F5F3");
  s.addText("Here both the score and the identity go up for the bad settings. The giveaway is the gap structure. 29 separate insert or delete events in one small protein is not believable.",
    { x: 9.62, y: 2.12, w: 2.86, h: ih - 0.35, fontFace: "Arial", fontSize: 14, color: INK, valign: "top" });
  const stat = (x, label, a, b) => {
    card(s, x, 5.55, 3.9, 1.35, "FFFFFF");
    s.addText(label, { x: x+0.22, y: 5.68, w: 3.5, h: 0.35, margin: 0, fontFace: "Arial", fontSize: 12.5, color: MUTED });
    s.addText([{ text: a + "  ", options: { color: MUTED } }, { text: "→ " + b, options: { color: RED, bold: true } }],
      { x: x+0.22, y: 6.02, w: 3.5, h: 0.7, margin: 0, fontFace: "Cambria", fontSize: 26, valign: "middle" });
  };
  stat(0.6, "raw score", "419", "497");
  stat(4.75, "separate gap runs", "1", "29");
  s.addNotes("Stronger than the toy: even identity favors the wrong answer here. Gap structure is what exposes it.");
})();

// ---------- Slide 6: DNA vs protein ----------
(() => {
  const s = contentSlide();
  header(s, "DNA vs protein", "Why some mutations hide", GREEN);
  const ih = 4.6, iw = ih * 1.775; // fig3 codon bar
  s.addImage({ path: "figures/fig3_sars_codons.png", x: 0.6, y: 1.95, w: iw, h: ih });
  card(s, 9.1, 1.95, 3.6, ih, "F2F5F3");
  s.addText("SARS-CoV vs SARS-CoV-2 spike", { x: 9.35, y: 2.15, w: 3.1, h: 0.7, margin: 0, fontFace: "Cambria", fontSize: 16, bold: true, color: INK });
  s.addText([
    { text: "I aligned the spike gene as DNA, then as the protein it codes for.", options: { breakLine: true, paraSpaceAfter: 10 } },
    { text: "438 of the 713 DNA changes are synonymous. They change the DNA but not the protein.", options: { breakLine: true, paraSpaceAfter: 10, color: GREEN, bold: true } },
    { text: "The mutations that matter are the 275 that change the amino acid.", options: {} },
  ], { x: 9.35, y: 3.0, w: 3.1, h: 3.4, fontFace: "Arial", fontSize: 14, color: INK, valign: "top" });
  s.addNotes("The centerpiece. Synonymous changes show up in DNA and vanish in protein because the genetic code is redundant.");
})();

// ---------- Slide 7: 3D ----------
(() => {
  const s = contentSlide();
  header(s, "Stretch: seeing it in 3D", "Where the changes land on the real spike", RED);
  const iw = 8.2, ih = iw / 1.812; // fig6 true aspect (real spike trimer, two views)
  s.addImage({ path: "figures/fig6_spike_3d.png", x: (13.3 - iw) / 2, y: 1.7, w: iw, h: ih });
  card(s, 0.6, 6.15, 12.1, 1.05, "F2F5F3");
  s.addText([
    { text: "The same 297 changed amino acids (red) drawn on the real spike structure. Most of it is gray and conserved. ", options: { color: INK } },
    { text: "The red clusters at the head, the part that grabs human cells and that antibodies target. ", options: { color: GREEN, bold: true } },
    { text: "The notebook has this as a picture and as a viewer you can rotate.", options: { color: MUTED, italic: true } },
  ], { x: 0.85, y: 6.28, w: 11.6, h: 0.8, fontFace: "Arial", fontSize: 13.5, valign: "middle" });
  s.addNotes("Switch to the notebook here and rotate the interactive 3D structure. Red is changed, gray is conserved.");
})();

// ---------- Slide: more views of the SARS changes (domain map + codon bar) ----------
(() => {
  const s = contentSlide();
  header(s, "The SARS changes", "Two more views along the sequence", GREEN);
  const dw = 12.0, dh = dw / 4.455; // fig4 domain map, wide
  s.addImage({ path: "figures/fig4_spike_domains.png", x: (13.3 - dw) / 2, y: 1.85, w: dw, h: dh });
  s.addImage({ path: "figures/fig10_conservation_line.png", x: 0.7, y: 4.75, w: 5.0, h: 5.0 / 2.195 });
  card(s, 6.0, 4.85, 6.7, 2.15, "F2F5F3");
  s.addText([
    { text: "Top: ", options: { bold: true, color: INK } },
    { text: "the 297 changed amino acids as red ticks along the 1273-amino-acid chain, with the functional regions marked. ", options: { breakLine: true, paraSpaceAfter: 10, color: INK } },
    { text: "Bottom: ", options: { bold: true, color: INK } },
    { text: "percent identity in a sliding window. The S2 fusion stalk stays high, while the host-facing NTD and receptor-binding domain dip the most.", options: { color: GREEN, bold: true } },
  ], { x: 6.25, y: 5.0, w: 6.2, h: 1.85, fontFace: "Arial", fontSize: 13.5, valign: "top" });
  s.addNotes("Two linear views that complement the 3D: where along the chain the changes fall, and the sliding-window identity.");
})();

// ---------- Slide 8: quality answer ----------
(() => {
  const s = contentSlide();
  header(s, "The central question", "How do you judge a good alignment?", GREEN);
  s.addText("The score cannot be the judge. It is defined by the knobs you picked, so raising the match reward makes any alignment score higher without making it more correct.",
    { x: 0.6, y: 1.95, w: 12.1, h: 0.9, fontFace: "Arial", fontSize: 16, color: INK, valign: "top" });
  const items = [
    ["Percent identity", "does not depend on your settings"],
    ["Gap structure", "a few long gaps beat many scattered ones"],
    ["Conserved regions", "do the known important parts line up"],
    ["Better than chance", "is the score real or just luck"],
    ["Independent evidence", "does a pro tool or the 3D structure agree"],
    ["Believable history", "are the implied mutations plausible"],
  ];
  items.forEach((it, i) => {
    const col = i % 2, row = Math.floor(i / 2);
    const x = 0.6 + col * 6.15, y = 3.05 + row * 1.15;
    card(s, x, y, 5.9, 1.0, "FFFFFF");
    s.addShape(pres.shapes.OVAL, { x: x+0.22, y: y+0.32, w: 0.34, h: 0.34, fill: { color: GREEN } });
    s.addText(String(i+1), { x: x+0.22, y: y+0.32, w: 0.34, h: 0.34, margin: 0, align: "center", valign: "middle", fontFace: "Arial", fontSize: 14, bold: true, color: "FFFFFF" });
    s.addText([
      { text: it[0] + ".  ", options: { bold: true, color: INK } },
      { text: it[1], options: { color: MUTED } },
    ], { x: x+0.72, y: y+0.15, w: 5.0, h: 0.7, margin: 0, fontFace: "Arial", fontSize: 14, valign: "middle" });
  });
  s.addNotes("This is the graded core. End on the bottom line: a good alignment is one whose story is biologically believable and holds up across reasonable settings, not the one with the biggest number.");
})();

// ---------- Slide 9: gap problem ----------
(() => {
  const s = contentSlide();
  header(s, "Stretch: the gap problem", "One gap penalty cannot do the job", RED);
  s.addImage({ path: "figures/fig5_gap_problem.png", x: 0.6, y: 1.85, sizing: { type: "contain", w: 12.1, h: 3.6 } });
  card(s, 0.6, 5.75, 12.1, 1.3, "FBEDEB");
  s.addText([
    { text: "With a flat penalty, a length-9 gap costs exactly the same as nine length-1 gaps, so no single value can both allow a real insertion and block scattered gaps. ", options: { color: INK } },
    { text: "Affine gaps fix it: opening a gap is expensive, extending one is cheap.", options: { color: RED, bold: true } },
  ], { x: 0.85, y: 5.9, w: 11.6, h: 1.0, fontFace: "Arial", fontSize: 14.5, valign: "middle" });
  s.addNotes("Explain that a real insertion is one event, so extending it should be cheap. Affine scoring matches the biology.");
})();

// ---------- Slide 10: sanity check ----------
(() => {
  const s = contentSlide();
  header(s, "Stretch: sanity check", "Does my aligner agree with the pros?", GREEN);
  card(s, 0.6, 2.2, 12.1, 2.2, "F2F5F3");
  s.addText([
    { text: "777.0", options: { fontSize: 54, bold: true, color: GREEN } },
    { text: "  optimal score", options: { fontSize: 20, color: MUTED } },
  ], { x: 1.1, y: 2.6, w: 5.6, h: 1.4, fontFace: "Cambria", valign: "middle", margin: 0 });
  s.addText([
    { text: "99.32%", options: { fontSize: 54, bold: true, color: GREEN } },
    { text: "  percent identity", options: { fontSize: 20, color: MUTED } },
  ], { x: 6.9, y: 2.6, w: 5.6, h: 1.4, fontFace: "Cambria", valign: "middle", margin: 0 });
  s.addText("My engine and Biopython (a tool professionals use) gave the exact same optimal score and identity on the same sequences with the same settings. Same optimum means the engine is correct. When there is a tie, two tools can show different-looking alignments that still score the same.",
    { x: 0.6, y: 4.7, w: 12.1, h: 1.8, fontFace: "Arial", fontSize: 16, color: INK, valign: "top" });
  s.addNotes("Independent confirmation the engine is right, not just plausible.");
})();

// ---------- Slide 11: closing ----------
(() => {
  const s = pres.addSlide();
  s.background = { color: DARK };
  s.addText("The score rewards what you tell it to.", { x: 1.0, y: 2.5, w: 11.3, h: 0.9, fontFace: "Cambria", fontSize: 36, bold: true, color: "FFFFFF" });
  s.addText("Biology decides what is true.", { x: 1.0, y: 3.4, w: 11.3, h: 0.9, fontFace: "Cambria", fontSize: 36, bold: true, color: "8FD1A6" });
  s.addText("A tunable aligner, a good-versus-bad story, the DNA-versus-protein comparison, and a 3D view of the mutations that matter. Everything runs on my own alignment code.",
    { x: 1.0, y: 4.6, w: 10.8, h: 1.0, fontFace: "Arial", fontSize: 17, color: "CADCD6" });
  s.addText("Thank you. Questions welcome.", { x: 1.0, y: 5.9, w: 11, h: 0.5, fontFace: "Arial", fontSize: 15, color: "9DB3AC" });
  s.addNotes("Close on the one-liner, then invite the pointed questions.");
})();

pres.writeFile({ fileName: process.env.DECK_OUT || "SequenceAligner_Slides.pptx" }).then(f => console.log("wrote", f));
