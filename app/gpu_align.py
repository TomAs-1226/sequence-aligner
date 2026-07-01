"""
GPU-accelerated batched alignment scores using PyTorch.

The normal aligner does one pair at a time on the CPU. Some tasks need the SAME
query aligned against thousands of other sequences (for example the significance
test, which aligns against hundreds of shuffled sequences). Those are all
independent, so they run in parallel on the GPU.

This computes global (Needleman-Wunsch) alignment SCORES for one query against a
batch of equal-length targets, using the anti-diagonal wavefront: every cell on
one diagonal is independent, so a whole diagonal is one vectorized GPU step.

Falls back to the CPU automatically if no GPU is found.
"""
import torch


def best_device():
    if hasattr(torch, "xpu") and torch.xpu.is_available():
        return "xpu"            # Intel GPU
    if torch.cuda.is_available():
        return "cuda"           # NVIDIA GPU
    return "cpu"


def sync(device):
    if device == "xpu":
        torch.xpu.synchronize()
    elif device == "cuda":
        torch.cuda.synchronize()


def build_matrix_tensor(matrix, alphabet):
    """Turn a {(a,b): score} dict into a [A, A] tensor and a char->index map."""
    idx = {c: i for i, c in enumerate(alphabet)}
    A = len(alphabet)
    M = torch.zeros((A, A), dtype=torch.float32)
    for a in alphabet:
        for b in alphabet:
            M[idx[a], idx[b]] = matrix.get((a, b), matrix.get((a, "X"), 0))
    return M, idx


def encode(seq, idx):
    x = idx.get("X", 0)
    return [idx.get(c, x) for c in seq]


def batched_global_scores(query, targets, matrix, gap, alphabet, device=None):
    """Global-alignment scores of `query` against each target (all targets must
    be the same length). Returns a list of scores. `matrix` is a {(a,b):int} dict."""
    if device is None:
        device = best_device()
    M, idx = build_matrix_tensor(matrix, alphabet)
    SM = M.to(device)
    n = len(query)
    m = len(targets[0])
    B = len(targets)
    qi = torch.tensor(encode(query, idx), device=device)                 # [n]
    ti = torch.tensor([encode(t, idx) for t in targets], device=device)  # [B, m]

    S = torch.empty((B, n + 1, m + 1), device=device, dtype=torch.float32)
    S[:, 0, :] = -gap * torch.arange(m + 1, device=device, dtype=torch.float32)
    S[:, :, 0] = -gap * torch.arange(n + 1, device=device, dtype=torch.float32).unsqueeze(0)

    for d in range(2, n + m + 1):
        i_lo = max(1, d - m)
        i_hi = min(n, d - 1)
        if i_lo > i_hi:
            continue
        i_idx = torch.arange(i_lo, i_hi + 1, device=device)   # [L]
        j_idx = d - i_idx                                      # [L]
        qrow = qi[i_idx - 1]                                   # [L]
        tcol = ti[:, j_idx - 1]                                # [B, L]
        sub = SM[qrow.unsqueeze(0).expand(B, -1), tcol]        # [B, L]
        diag = S[:, i_idx - 1, j_idx - 1] + sub
        up = S[:, i_idx - 1, j_idx] - gap
        left = S[:, i_idx, j_idx - 1] - gap
        S[:, i_idx, j_idx] = torch.maximum(torch.maximum(diag, up), left)

    sync(device)
    return S[:, n, m].tolist()


# protein alphabet that matches the BLOSUM62 / PAM250 matrices
PROTEIN_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZX*"
DNA_ALPHABET = "ACGTN"
