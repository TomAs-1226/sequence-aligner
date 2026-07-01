"""
Train a small protein language model (a transformer) from scratch on the GPU.

This is the heavy-compute showcase. The model reads protein sequences and learns
to predict masked-out amino acids (the same idea as BERT, but for proteins). Once
trained, the encoder turns any protein into a vector (an "embedding"). Proteins
in the same family end up near each other in this space, even when their
sequences are too different for a plain alignment to notice.

Runs on the Intel Arc GPU (xpu), or CUDA / CPU if that is what is present.
"""
import time
import urllib.parse
import urllib.request

import numpy as np
import torch
import torch.nn as nn

AA = "ACDEFGHIKLMNPQRSTVWY"
TOK = {c: i for i, c in enumerate(AA)}
PAD, MASK, UNK = 20, 21, 22
VOCAB = 23
MAXLEN = 200


def device():
    if hasattr(torch, "xpu") and torch.xpu.is_available():
        return "xpu"
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"


def sync(d):
    if d == "xpu":
        torch.xpu.synchronize()
    elif d == "cuda":
        torch.cuda.synchronize()


def encode(s):
    return [TOK.get(c, UNK) for c in s]


def fetch(query, cap):
    url = "https://rest.uniprot.org/uniprotkb/stream?query=%s&format=fasta" % urllib.parse.quote(query)
    seqs, cur = [], None
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=180) as r:
            for raw in r:
                line = raw.decode("utf-8", "ignore").strip()
                if line.startswith(">"):
                    if cur is not None and 40 <= len(cur) <= MAXLEN and set(cur) <= set(AA):
                        seqs.append(cur)
                    cur = ""
                    if len(seqs) >= cap:
                        break
                elif cur is not None:
                    cur += line
    except Exception as ex:
        print("fetch warning:", ex)
    return seqs


D = device()
print("device:", D, "|", (torch.xpu.get_device_name(0) if D == "xpu" else D))

print("Fetching training corpus...")
corpus = fetch("reviewed:true AND length:[40 TO 200]", 16000)
print(f"  corpus: {len(corpus)} proteins")

print("Fetching labeled families for validation...")
FAMS = {
    "globin": 'protein_name:globin',
    "cytochrome c": 'protein_name:"cytochrome c"',
    "protein kinase": 'protein_name:"protein kinase"',
    "histone": 'protein_name:histone',
    "ribonuclease": 'protein_name:ribonuclease',
}
labeled = []
for fam, q in FAMS.items():
    fs = fetch("reviewed:true AND (%s) AND length:[40 TO 200]" % q, 130)
    for s in fs[:120]:
        labeled.append((s, fam))
    print(f"  {fam}: {len([1 for _, f in labeled if f == fam])}")


class PLM(nn.Module):
    def __init__(self, d=192, nhead=6, layers=4, ff=512):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, d)
        self.pos = nn.Embedding(MAXLEN + 1, d)
        layer = nn.TransformerEncoderLayer(d, nhead, ff, batch_first=True, dropout=0.1)
        self.enc = nn.TransformerEncoder(layer, layers)
        self.head = nn.Linear(d, VOCAB)

    def hidden(self, x, pad):
        L = x.size(1)
        h = self.emb(x) + self.pos(torch.arange(L, device=x.device))[None]
        return self.enc(h, src_key_padding_mask=pad)

    def forward(self, x, pad):
        return self.head(self.hidden(x, pad))

    def embed(self, x, pad):
        h = self.hidden(x, pad)
        w = (~pad).unsqueeze(-1).float()
        return (h * w).sum(1) / w.sum(1).clamp(min=1)


def batchify(seqs, bs=32):
    order = np.random.permutation(len(seqs))
    for i in range(0, len(seqs), bs):
        chunk = [seqs[k] for k in order[i:i + bs]]
        L = max(len(s) for s in chunk)
        x = torch.full((len(chunk), L), PAD, dtype=torch.long)
        pad = torch.ones((len(chunk), L), dtype=torch.bool)
        for j, s in enumerate(chunk):
            ids = encode(s)
            x[j, :len(ids)] = torch.tensor(ids)
            pad[j, :len(ids)] = False
        yield x, pad


np.random.seed(0)
torch.manual_seed(0)
model = PLM().to(D)
opt = torch.optim.AdamW(model.parameters(), lr=3e-4)
lossf = nn.CrossEntropyLoss(ignore_index=-100)
print(f"model: {sum(p.numel() for p in model.parameters())/1e6:.1f}M params")

split = int(0.97 * len(corpus))
train, val = corpus[:split], corpus[split:]

TIME_BUDGET = 720  # seconds of training
t0 = time.time()
step = 0
epoch = 0
sync(D)
while time.time() - t0 < TIME_BUDGET:
    epoch += 1
    model.train()
    for x, pad in batchify(train, 32):
        x, pad = x.to(D), pad.to(D)
        # mask 15% of real (non-pad) positions
        real = ~pad
        m = (torch.rand_like(x, dtype=torch.float) < 0.15) & real
        target = torch.where(m, x, torch.full_like(x, -100))
        xin = torch.where(m, torch.full_like(x, MASK), x)
        logits = model(xin, pad)
        loss = lossf(logits.reshape(-1, VOCAB), target.reshape(-1))
        opt.zero_grad(); loss.backward(); opt.step()
        step += 1
        if time.time() - t0 >= TIME_BUDGET:
            break
    # quick val accuracy on masked positions
    model.eval()
    cor = tot = 0
    with torch.no_grad():
        for x, pad in batchify(val, 32):
            x, pad = x.to(D), pad.to(D)
            real = ~pad
            m = (torch.rand_like(x, dtype=torch.float) < 0.15) & real
            xin = torch.where(m, torch.full_like(x, MASK), x)
            pred = model(xin, pad).argmax(-1)
            cor += (pred[m] == x[m]).sum().item(); tot += int(m.sum().item())
    print(f"epoch {epoch} (step {step}, {time.time()-t0:.0f}s): masked-residue accuracy {100*cor/max(tot,1):.1f}%  (random = {100/20:.0f}%)")

sync(D)
print(f"trained {step} steps in {time.time()-t0:.0f}s on the {D}")

# ---- validate: do family embeddings cluster? ----
model.eval()
embs, labs = [], []
with torch.no_grad():
    for s, fam in labeled:
        ids = torch.tensor([encode(s)], device=D)
        pad = torch.zeros_like(ids, dtype=torch.bool)
        embs.append(model.embed(ids, pad)[0].cpu().numpy()); labs.append(fam)
E = np.array(embs)
E = (E - E.mean(0)) / (E.std(0) + 1e-6)
# PCA to 2D
U, S, Vt = np.linalg.svd(E - E.mean(0), full_matrices=False)
P = (E - E.mean(0)) @ Vt[:2].T

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.5, 6))
colors = {"globin": "#c43d2e", "cytochrome c": "#2E7D46", "protein kinase": "#1565c0",
          "histone": "#9c27b0", "ribonuclease": "#d98a26"}
for fam in FAMS:
    idx = [i for i, l in enumerate(labs) if l == fam]
    ax.scatter(P[idx, 0], P[idx, 1], s=18, color=colors[fam], label=fam, alpha=0.8)
ax.set_title("protein families in the language model's embedding space (PCA)")
ax.legend(); ax.set_xticks([]); ax.set_yticks([])
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.savefig("../figures/fig11_plm_embeddings.png", dpi=150, bbox_inches="tight")
print("saved figures/fig11_plm_embeddings.png")

# simple cluster-quality score: nearest-neighbour same-family rate
from math import inf
same = 0
for i in range(len(P)):
    best, bj = inf, -1
    for j in range(len(P)):
        if i == j:
            continue
        dsq = float(((E[i] - E[j]) ** 2).sum())
        if dsq < best:
            best, bj = dsq, j
    if labs[bj] == labs[i]:
        same += 1
print(f"nearest-neighbour same-family rate: {100*same/len(P):.0f}% (chance ~ {100//len(FAMS)}%)")

torch.save(model.state_dict(), "plm.pt")
print("saved plm.pt")
