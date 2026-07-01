"""
Train a small neural network to predict per-residue secondary structure
(helix / sheet / coil) from a protein sequence, and export it to ONNX so the
WinUI app can run it on the GPU or NPU.

Training data is REAL: we download experimental protein structures from the PDB
and read the secondary structure off the coordinates with biotite. So the model
learns structure-from-sequence, a genuine (small) AI task.
"""
import io
import json
import os
import time
import urllib.parse
import urllib.request

import numpy as np
import torch
import torch.nn as nn
import biotite.structure as struc
import biotite.structure.io.pdb as pdb

AA = "ACDEFGHIKLMNPQRSTVWY"
AAI = {c: i for i, c in enumerate(AA)}
X = len(AA)
THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}
CACHE = "ss_cache"
os.makedirs(CACHE, exist_ok=True)


def enc(seq):
    return [AAI.get(c, X) for c in seq]


def get_pdb_ids(rows=400):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.selected_polymer_entity_types", "operator": "exact_match", "value": "Protein (only)"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "entity_poly.rcsb_sample_sequence_length", "operator": "range", "value": {"from": 60, "to": 250}}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less", "value": 2.0}}]},
        "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": rows}, "results_content_type": ["experimental"]}}
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + urllib.parse.quote(json.dumps(q))
    res = json.loads(urllib.request.urlopen(url, timeout=45).read().decode())
    return [r["identifier"] for r in res["result_set"]]


def get_seq_ss(pid):
    key = os.path.join(CACHE, pid + ".pdb")
    if os.path.exists(key):
        txt = open(key).read()
    else:
        try:
            txt = urllib.request.urlopen("https://files.rcsb.org/download/%s.pdb" % pid, timeout=25).read().decode("utf-8", "ignore")
            open(key, "w").write(txt)
        except Exception:
            return None
    try:
        st = pdb.PDBFile.read(io.StringIO(txt)).get_structure(model=1)
        st = st[struc.filter_amino_acids(st)]
        if len(st) == 0:
            return None
        st = st[st.chain_id == st.chain_id[0]]
        sse = struc.annotate_sse(st)
        _, names = struc.get_residues(st)
        if len(sse) != len(names) or len(names) < 40:
            return None
        seq = "".join(THREE2ONE.get(n, "X") for n in names)
        m = {"a": 0, "b": 1, "c": 2}
        lab = np.array([m.get(x, 2) for x in sse], dtype=np.int64)
        return seq, lab
    except Exception:
        return None


print("Getting PDB IDs...")
ids = get_pdb_ids(400)
print(f"  {len(ids)} candidate structures")

data = []
t0 = time.time()
for pid in ids:
    if len(data) >= 240:
        break
    r = get_seq_ss(pid)
    if r:
        data.append(r)
    if len(data) % 40 == 0 and len(data) > 0:
        print(f"  {len(data)} structures ({time.time()-t0:.0f}s)")
print(f"dataset: {len(data)} proteins, {sum(len(s) for s, _ in data)} residues")

rs = np.random.RandomState(0)
rs.shuffle(data)
split = int(0.85 * len(data))
train, test = data[:split], data[split:]


class SSNet(nn.Module):
    def __init__(self, vocab=X + 1, emb=32, hid=64, k=11):
        super().__init__()
        self.emb = nn.Embedding(vocab, emb)
        self.c1 = nn.Conv1d(emb, hid, k, padding=k // 2)
        self.c2 = nn.Conv1d(hid, hid, k, padding=k // 2)
        self.c3 = nn.Conv1d(hid, hid, 3, padding=1)
        self.out = nn.Conv1d(hid, 3, 1)

    def forward(self, x):
        e = self.emb(x).transpose(1, 2)
        h = torch.relu(self.c1(e))
        h = torch.relu(self.c2(h))
        h = torch.relu(self.c3(h))
        return self.out(h).transpose(1, 2)


model = SSNet()
opt = torch.optim.Adam(model.parameters(), lr=1e-3)
lossf = nn.CrossEntropyLoss()
for epoch in range(40):
    model.train()
    rs.shuffle(train)
    tot = 0.0
    for s, lab in train:
        x = torch.tensor([enc(s)], dtype=torch.long)
        logits = model(x)[0]
        loss = lossf(logits, torch.tensor(lab, dtype=torch.long))
        opt.zero_grad(); loss.backward(); opt.step()
        tot += loss.item()
    if epoch % 8 == 0 or epoch == 39:
        model.eval()
        cor = tot_res = 0
        with torch.no_grad():
            for s, lab in test:
                pred = model(torch.tensor([enc(s)], dtype=torch.long))[0].argmax(1).numpy()
                cor += int((pred == lab).sum()); tot_res += len(lab)
        print(f"epoch {epoch}: loss {tot/len(train):.3f}, test Q3 accuracy {100*cor/tot_res:.1f}%")

model.eval()
torch.onnx.export(model, torch.zeros(1, 40, dtype=torch.long), "ss_model.onnx",
                  input_names=["seq"], output_names=["logits"],
                  dynamic_axes={"seq": {1: "L"}, "logits": {1: "L"}}, opset_version=13,
                  dynamo=False)
print("exported ss_model.onnx (%d KB)" % (os.path.getsize("ss_model.onnx") // 1024))

import onnxruntime as ort
sess = ort.InferenceSession("ss_model.onnx", providers=["CPUExecutionProvider"])
s0 = test[0][0]
x0 = np.array([enc(s0)], dtype=np.int64)
op = sess.run(None, {"seq": x0})[0][0].argmax(1)
tp = model(torch.tensor(x0))[0].argmax(1).numpy()
print("ONNX == torch:", bool((op == tp).all()))
print("example prediction:", "".join("HEC"[i] for i in op[:60]))
