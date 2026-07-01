// Runs the real Python engine in the browser, off the main thread so the
// page never freezes. Loads Pyodide (CPython compiled to WebAssembly), then
// engine.py and web.py, exactly the same code as the desktop program.
importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.4/full/pyodide.js");

let pyodide = null;
let web = null;

async function init() {
  pyodide = await loadPyodide();
  const [engineSrc, webSrc] = await Promise.all([
    fetch("engine.py").then((r) => r.text()),
    fetch("web.py").then((r) => r.text()),
  ]);
  pyodide.FS.writeFile("engine.py", engineSrc);
  pyodide.FS.writeFile("web.py", webSrc);
  pyodide.runPython("import web");
  web = pyodide.pyimport("web");
  postMessage({ type: "ready" });
}

const ready = init().catch((err) =>
  postMessage({ type: "fatal", message: String(err) })
);

onmessage = async (ev) => {
  await ready;
  const msg = ev.data;
  try {
    if (msg.type === "align") {
      const out = web.web_align_json(JSON.stringify(msg.params));
      postMessage({ type: "result", id: msg.id, data: JSON.parse(out) });
    } else if (msg.type === "sweep") {
      const out = web.sweep_json(JSON.stringify(msg.params));
      postMessage({ type: "sweep", id: msg.id, data: JSON.parse(out) });
    } else if (msg.type === "significance") {
      const out = web.significance_json(JSON.stringify(msg.params));
      postMessage({ type: "significance", id: msg.id, data: JSON.parse(out) });
    } else if (msg.type === "trace") {
      const out = web.trace_align_json(JSON.stringify(msg.params));
      postMessage({ type: "trace", id: msg.id, data: JSON.parse(out) });
    } else if (msg.type === "translate") {
      postMessage({ type: "translate", id: msg.id, data: web.do_translate(msg.seq) });
    } else if (msg.type === "revcomp") {
      postMessage({ type: "revcomp", id: msg.id, data: web.do_revcomp(msg.seq) });
    }
  } catch (err) {
    postMessage({ type: "error", id: msg.id, message: String(err) });
  }
};
