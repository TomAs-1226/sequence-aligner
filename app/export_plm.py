"""Export the trained protein language model's encoder to ONNX as an EMBEDDING
model (sequence -> one vector), for the app to run on the NPU/GPU via Windows ML.
Also precompute a small database of family embeddings the app can search against."""
import numpy as np
import torch
import torch.nn as nn
import urllib.parse
import urllib.request

AA = "ACDEFGHIKLMNPQRSTVWY"
TOK = {c: i for i, c in enumerate(AA)}
PAD, MASK, UNK, VOCAB, MAXLEN = 20, 21, 22, 23, 200


def encode(s):
    return [TOK.get(c, UNK) for c in s]


FIXED = 200  # the app pads/truncates every sequence to this length


class PLM(nn.Module):
    def __init__(self, d=192, nhead=6, layers=4, ff=512):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, d)
        self.pos = nn.Embedding(MAXLEN + 1, d)
        layer = nn.TransformerEncoderLayer(d, nhead, ff, batch_first=True, dropout=0.1)
        self.enc = nn.TransformerEncoder(layer, layers, enable_nested_tensor=False)
        self.head = nn.Linear(d, VOCAB)


class Embedder(nn.Module):
    """Fixed length 200 with a padding mask (1 = pad). Outputs an L2-normalized
    embedding, mean-pooled over the real residues."""
    def __init__(self, plm):
        super().__init__()
        self.plm = plm

    def forward(self, x, mask):                # x: [1,200] int64, mask: [1,200] float
        padb = mask > 0.5
        pos = torch.arange(FIXED, device=x.device)
        h = self.plm.emb(x) + self.plm.pos(pos)[None]
        h = self.plm.enc(h, src_key_padding_mask=padb)
        w = (1.0 - mask).unsqueeze(-1)
        v = (h * w).sum(dim=1) / w.sum(dim=1).clamp(min=1.0)
        return v / (v.norm(dim=1, keepdim=True) + 1e-6)


def pad_seq(seq):
    ids = encode(seq)[:FIXED]
    mask = [0.0] * len(ids) + [1.0] * (FIXED - len(ids))
    ids = ids + [PAD] * (FIXED - len(ids))
    return np.array(ids, dtype=np.int64), np.array(mask, dtype=np.float32)


model = PLM()
model.load_state_dict(torch.load("plm.pt", map_location="cpu"))
model.eval()
emb = Embedder(model).eval()


def embed(seq):
    ids, mask = pad_seq(seq)
    with torch.no_grad():
        return emb(torch.tensor([ids]), torch.tensor([mask]))[0].numpy()


# ---- export to ONNX (fixed shapes) ----
dx = torch.zeros(1, FIXED, dtype=torch.long)
dm = torch.zeros(1, FIXED, dtype=torch.float32)
torch.onnx.export(emb, (dx, dm), "plm_embed.onnx", input_names=["seq", "mask"],
                  output_names=["embedding"], opset_version=17, dynamo=False)
print("exported plm_embed.onnx")

import onnx
mg = onnx.load("plm_embed.onnx")
print("ONNX ops:", sorted({n.op_type for n in mg.graph.node}))

import onnxruntime as ort
sess = ort.InferenceSession("plm_embed.onnx", providers=["CPUExecutionProvider"])
ids, mask = pad_seq("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFF")
oe = sess.run(None, {"seq": np.array([ids]), "mask": np.array([mask])})[0]
te = embed("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFF")
print("ONNX == torch embedding:", bool(np.allclose(oe[0], te, atol=1e-4)))
print("embedding dim:", oe.shape[-1])


# ---- build a small family database the app can search against ----
def fetch(query, cap):
    url = "https://rest.uniprot.org/uniprotkb/stream?query=%s&format=fasta" % urllib.parse.quote(query)
    out, cur, name = [], None, None
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=120) as r:
            for raw in r:
                line = raw.decode("utf-8", "ignore").rstrip("\n")
                if line.startswith(">"):
                    if cur and 40 <= len(cur) <= MAXLEN and set(cur) <= set(AA):
                        out.append((name, cur))
                    name = line[1:].split(" ", 1)[-1][:40]
                    cur = ""
                    if len(out) >= cap:
                        break
                elif cur is not None:
                    cur += line.strip()
    except Exception as ex:
        print("fetch warn:", ex)
    return out


FAMS = {"globin": 'protein_name:globin', "cytochrome c": 'protein_name:"cytochrome c"',
        "protein kinase": 'protein_name:"protein kinase"', "histone": 'protein_name:histone',
        "ribonuclease": 'protein_name:ribonuclease', "insulin": 'protein_name:insulin',
        "lysozyme": 'protein_name:lysozyme', "myoglobin": 'protein_name:myoglobin'}
db = []
for fam, q in FAMS.items():
    for name, seq in fetch("reviewed:true AND (%s) AND length:[40 TO 200]" % q, 40):
        db.append((fam, name, embed(seq)))
print("database proteins:", len(db))

# write a compact text database: family<TAB>name<TAB>comma-separated-embedding
with open("../SequenceAlignerApp/plm_db.tsv", "w", encoding="utf-8") as f:
    for fam, name, v in db:
        f.write(fam + "\t" + name + "\t" + ",".join("%.4f" % x for x in v) + "\n")
print("wrote SequenceAlignerApp/plm_db.tsv (%d entries)" % len(db))
