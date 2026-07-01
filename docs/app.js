// Front end for the sequence aligner. Talks to worker.js, which runs the real
// Python engine. Keeps the UI responsive and renders the colored alignment.

const $ = (id) => document.getElementById(id);

const state = {
  seqtype: "dna", matrix: "BLOSUM62", mode: "global", gapmodel: "linear",
  match: 1, mismatch: 1, gap: 2, gap_open: 8, gap_extend: 1,
  seq1: "", seq2: "",
};

let worker = null;
let reqId = 0;      // latest request wins
let ready = false;

// ---------- worker ----------
function startWorker() {
  worker = new Worker("worker.js");
  worker.onmessage = (ev) => {
    const msg = ev.data;
    if (msg.type === "ready") {
      ready = true;
      $("overlay").classList.add("hidden");
      align();
    } else if (msg.type === "result") {
      if (msg.id === reqId) showResult(msg.data);
    } else if (msg.type === "translate") {
      applyTranslate(msg.data);
    } else if (msg.type === "revcomp") {
      $("seq2").value = msg.data; state.seq2 = msg.data; align();
    } else if (msg.type === "error") {
      $("stUsed").textContent = "Could not align: " + msg.message;
    } else if (msg.type === "fatal") {
      $("overlay").innerHTML =
        '<div class="big">The Python engine did not load</div>' +
        '<p>Please check your internet connection and refresh. Details: ' + msg.message + "</p>";
    }
  };
}

// ---------- align ----------
let timer = null;
function align(now) {
  if (!ready) return;
  state.seq1 = clean($("seq1").value, state.seqtype);
  state.seq2 = clean($("seq2").value, state.seqtype);
  clearTimeout(timer);
  const go = () => {
    $("stUsed").textContent = "aligning...";
    const id = ++reqId;
    worker.postMessage({ type: "align", id, params: { ...state } });
    drawDotplot(state.seq1, state.seq2, state.seqtype);
  };
  if (now) go(); else timer = setTimeout(go, 130);
}

function clean(s, type) {
  s = (s || "").toUpperCase();
  return type === "dna" ? s.replace(/[^ACGTUN]/g, "") : s.replace(/[^A-Z*]/g, "");
}

// ---------- render ----------
function showResult(d) {
  if (d.error) { $("aln").textContent = d.error; $("stScore").textContent = "-"; return; }
  $("stScore").textContent = d.score;
  $("stIdent").textContent = d.identity + "%";
  $("stGaps").textContent = d.gaps;
  $("stUsed").textContent =
    `${d.used}. Aligned ${d.len1} vs ${d.len2} letters into ${d.cols} columns.`;
  $("aln").innerHTML = renderAln(d.row1, d.match, d.row2);
  window._lastAln = d.row1 + "\n" + d.match + "\n" + d.row2;
}

function colorSeq(s, other) {
  let out = "";
  for (let k = 0; k < s.length; k++) {
    const c = s[k], o = other[k];
    const cls = c === "-" ? "g" : (c === o ? "m" : "x");
    out += `<span class="${cls}">${c}</span>`;
  }
  return out;
}

function renderAln(row1, match, row2) {
  const W = 60; let html = "";
  for (let i = 0; i < row1.length; i += W) {
    const a = row1.slice(i, i + W), m = match.slice(i, i + W), b = row2.slice(i, i + W);
    html += `<span class="lab">position ${i + 1}</span>\n`;
    html += colorSeq(a, b) + "\n";
    html += `<span class="bar">${m}</span>\n`;
    html += colorSeq(b, a) + "\n\n";
  }
  return html;
}

// ---------- dot plot ----------
function drawDotplot(s1, s2, type) {
  const cv = $("dotplot"), ctx = cv.getContext("2d"), S = cv.width;
  ctx.fillStyle = "#0c1a17"; ctx.fillRect(0, 0, S, S);
  const n = s1.length, m = s2.length;
  const k = type === "dna" ? 4 : 3;
  if (n < k || m < k) return;
  const idx = {};
  for (let j = 0; j <= m - k; j++) {
    const key = s2.slice(j, j + k);
    (idx[key] || (idx[key] = [])).push(j);
  }
  ctx.fillStyle = "#4cc27a";
  const sx = (S - 2) / n, sy = (S - 2) / m;
  const dot = Math.max(1, Math.min(sx, sy));
  let pts = 0; const MAX = 80000;
  for (let i = 0; i <= n - k && pts < MAX; i++) {
    const arr = idx[s1.slice(i, i + k)];
    if (!arr) continue;
    for (const j of arr) {
      ctx.fillRect(1 + i * sx, 1 + j * sy, dot, dot);
      if (++pts >= MAX) break;
    }
  }
}

// ---------- tools ----------
function applyTranslate(protPairJson) {
  const p = JSON.parse(protPairJson);
  $("seq1").value = p[0]; $("seq2").value = p[1];
  setType("protein"); state.gap = 10; syncSliders(); align(true);
}
$("btnTranslate").onclick = () => {
  if (state.seqtype !== "dna") { flash("btnTranslate", "already protein"); return; }
  worker.postMessage({ type: "translate", seq: $("seq1").value + "||" + $("seq2").value });
};
$("btnRevcomp").onclick = () => {
  if (state.seqtype !== "dna") { flash("btnRevcomp", "DNA only"); return; }
  worker.postMessage({ type: "revcomp", seq: $("seq2").value });
};
$("btnCopy").onclick = () => {
  if (!window._lastAln) return;
  navigator.clipboard.writeText(window._lastAln).then(() => flash("btnCopy", "copied"));
};
function flash(id, txt) { const b = $(id), o = b.textContent; b.textContent = txt; setTimeout(() => (b.textContent = o), 1100); }

// ---------- controls wiring ----------
function seg(segId, key, after) {
  const el = $(segId);
  el.querySelectorAll("button").forEach((b) => {
    b.onclick = () => {
      el.querySelectorAll("button").forEach((x) => x.classList.remove("on"));
      b.classList.add("on");
      state[key] = b.dataset.v;
      if (after) after();
      align(true);
    };
  });
}

function setType(v) {
  state.seqtype = v;
  $("typeSeg").querySelectorAll("button").forEach((b) =>
    b.classList.toggle("on", b.dataset.v === v));
  $("dnaSliders").style.display = v === "dna" ? "" : "none";
  $("matrixWrap").style.display = v === "protein" ? "" : "none";
}

function updateGapModel() {
  const affine = state.gapmodel === "affine";
  $("gapRow").style.display = affine ? "none" : "";
  $("affineSliders").style.display = affine ? "" : "none";
}

function updateAffineAvailability() {
  // affine gaps are wired for global mode; for other modes we use linear
  const affineBtn = $("gapSeg").querySelector('button[data-v="affine"]');
  const allow = state.mode === "global";
  affineBtn.style.opacity = allow ? "1" : "0.4";
  affineBtn.disabled = !allow;
  if (!allow && state.gapmodel === "affine") {
    state.gapmodel = "linear";
    $("gapSeg").querySelectorAll("button").forEach((x) =>
      x.classList.toggle("on", x.dataset.v === "linear"));
    updateGapModel();
  }
}

function slider(id) {
  const el = $(id);
  el.oninput = () => {
    state[id] = parseFloat(el.value);
    $(id + "Val").textContent = el.value;
    align();
  };
}

function syncSliders() {
  ["match", "mismatch", "gap", "gap_open", "gap_extend"].forEach((id) => {
    $(id).value = state[id]; $(id + "Val").textContent = state[id];
  });
}

// ---------- samples ----------
function loadSamples() {
  const sel = $("sample");
  window.SAMPLES.forEach((s, i) => {
    const o = document.createElement("option");
    o.value = i; o.textContent = s.label; sel.appendChild(o);
  });
  sel.onchange = () => applySample(window.SAMPLES[+sel.value]);
  applySample(window.SAMPLES[0]);
}

function applySample(s) {
  $("seq1").value = s.seq1; $("seq2").value = s.seq2;
  $("sampleNote").textContent = s.note || "";
  setType(s.type);
  state.gap = s.type === "dna" ? 2 : 10;
  syncSliders();
  align(true);
}

// ---------- presets ----------
$("presetGood").onclick = () => {
  state.gapmodel = "linear";
  $("gapSeg").querySelectorAll("button").forEach((x) => x.classList.toggle("on", x.dataset.v === "linear"));
  updateGapModel();
  state.match = 1; state.mismatch = 1; state.gap = state.seqtype === "dna" ? 2 : 10;
  syncSliders(); align(true);
};
$("presetBad").onclick = () => {
  state.gapmodel = "linear";
  $("gapSeg").querySelectorAll("button").forEach((x) => x.classList.toggle("on", x.dataset.v === "linear"));
  updateGapModel();
  state.match = 1; state.mismatch = 1; state.gap = 0;
  syncSliders(); align(true);
};

// ---------- init ----------
seg("typeSeg", "seqtype", () => { state.gap = state.seqtype === "dna" ? 2 : 10; setType(state.seqtype); syncSliders(); });
seg("matrixSeg", "matrix");
seg("modeSeg", "mode", updateAffineAvailability);
seg("gapSeg", "gapmodel", () => { updateGapModel(); updateAffineAvailability(); });
["match", "mismatch", "gap", "gap_open", "gap_extend"].forEach(slider);
$("seq1").oninput = () => align();
$("seq2").oninput = () => align();

updateGapModel();
loadSamples();
startWorker();
