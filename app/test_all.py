"""Comprehensive bug check for the alignment engine and the GPU aligner.
Cross-checks against brute-force references and known values. Prints PASS/FAIL."""
import random
import sys
from functools import lru_cache

import engine as e
import gpu_align as g

sys.path.insert(0, "../Alignment")

fails = []


def check(name, cond):
    if not cond:
        fails.append(name)
        print("  FAIL:", name)


# ---------- brute-force references ----------
def bf_global(s1, s2, match, mismatch, gap):
    @lru_cache(None)
    def f(i, j):
        if i == 0 and j == 0:
            return 0.0
        best = float("-inf")
        if i > 0 and j > 0:
            best = max(best, f(i - 1, j - 1) + (match if s1[i - 1] == s2[j - 1] else -mismatch))
        if i > 0:
            best = max(best, f(i - 1, j) - gap)
        if j > 0:
            best = max(best, f(i, j - 1) - gap)
        return best
    return f(len(s1), len(s2))


def bf_local(s1, s2, match, mismatch, gap):
    best = 0.0
    for a in range(len(s1) + 1):
        for b in range(a, len(s1) + 1):
            for c in range(len(s2) + 1):
                for d in range(c, len(s2) + 1):
                    best = max(best, bf_global(s1[a:b], s2[c:d], match, mismatch, gap))
    return best


def rand_dna(n):
    return "".join(random.choice("ACGT") for _ in range(n))


random.seed(7)

# ---------- 1. global vs brute force ----------
print("[1] global alignment vs brute force")
for _ in range(120):
    n1, n2 = random.randint(0, 6), random.randint(0, 6)
    if n1 == 0 or n2 == 0:
        continue
    s1, s2 = rand_dna(n1), rand_dna(n2)
    m, mm, gp = random.choice([1, 2]), random.choice([1, 2, 3]), random.choice([1, 2])
    r = e.align_global(s1, s2, e.dna_scorer(m, mm), gp)
    check(f"global {s1}/{s2} m{m}mm{mm}g{gp}", abs(r[2] - bf_global(s1, s2, m, mm, gp)) < 1e-6)
    # score of returned alignment must equal reported score
    check(f"global-score-consistent {s1}/{s2}", abs(e.alignment_score(r[0], r[1], e.dna_scorer(m, mm), gp) - r[2]) < 1e-6)
    # rows equal length, removing gaps gives back originals
    check(f"global-rows {s1}/{s2}", len(r[0]) == len(r[1]) and r[0].replace("-", "") == s1 and r[1].replace("-", "") == s2)

# ---------- 2. local vs brute force ----------
print("[2] local alignment vs brute force")
for _ in range(80):
    s1, s2 = rand_dna(random.randint(1, 6)), rand_dna(random.randint(1, 6))
    m, mm, gp = random.choice([1, 2]), random.choice([1, 2]), random.choice([1, 2])
    r = e.align_local(s1, s2, e.dna_scorer(m, mm), gp)
    check(f"local {s1}/{s2}", abs(r[2] - bf_local(s1, s2, m, mm, gp)) < 1e-6)

# ---------- 3. banded == full global (band large) ----------
print("[3] banded vs full global")
for _ in range(60):
    s1 = rand_dna(random.randint(20, 60))
    s2 = "".join(c for c in s1 if random.random() > 0.1)  # deletions
    s2 = "".join((c if random.random() > 0.1 else random.choice("ACGT")) for c in s2)  # subs
    full = e.align_global(s1, s2, e.dna_scorer(2, 3), 2)[2]
    band = e.align_banded(s1, s2, e.dna_scorer(2, 3), 2, band=max(len(s1), len(s2)))[2]
    check("banded==full", abs(full - band) < 1e-6)

# ---------- 4. affine with open==extend == linear global ----------
print("[4] affine (open==extend) == linear global")
for _ in range(60):
    s1, s2 = rand_dna(random.randint(1, 7)), rand_dna(random.randint(1, 7))
    gp = random.choice([1, 2])
    lin = e.align_global(s1, s2, e.dna_scorer(2, 3), gp)[2]
    aff = e.align_global_affine(s1, s2, e.dna_scorer(2, 3), gp, gp)[2]
    check("affine==linear when open==extend", abs(lin - aff) < 1e-6)

# ---------- 5. semi-global properties ----------
print("[5] semi-global (free end gaps)")
for _ in range(60):
    core = rand_dna(random.randint(3, 8))
    pre, suf = rand_dna(random.randint(0, 5)), rand_dna(random.randint(0, 5))
    s2 = pre + core + suf
    # core vs (pre+core+suf): free end gaps => score is the perfect match of core
    r = e.align_semiglobal(core, s2, e.dna_scorer(2, 3), 4)
    check("semiglobal free ends", abs(r[2] - 2 * len(core)) < 1e-6)
    # semi-global score >= global score
    a, b = rand_dna(random.randint(1, 7)), rand_dna(random.randint(1, 7))
    sg = e.align_semiglobal(a, b, e.dna_scorer(2, 3), 2)[2]
    gl = e.align_global(a, b, e.dna_scorer(2, 3), 2)[2]
    check("semiglobal>=global", sg >= gl - 1e-6)

# ---------- 6. readouts, translate, reverse complement ----------
print("[6] readouts / translate / revcomp")
check("identity", abs(e.percent_identity("AC-GT", "ACTGT") - 100.0) < 1e-9)  # gap col excluded -> 4/4
check("identity2", abs(e.percent_identity("ACGT", "ACCT") - 75.0) < 1e-9)
check("gapstats", e.gap_stats("A--GT-", "ABCGTX") == (2, [2, 1]))
check("translate", e.translate("ATGGCCTGA") == "MA")
check("translate-stop", e.translate("ATGTAAGGG") == "M")
check("translate-lower/u", e.translate("augGCC".replace("a","A")) == e.translate("ATGGCC"))
check("revcomp", e.reverse_complement("ATTTAC") == "GTAAAT")
check("revcomp-involution", all(e.reverse_complement(e.reverse_complement(x := rand_dna(20))) == x for _ in range(5)))
# translate reproduces the provided spike proteins
for nt, aa in [("SARS-CoV_genome_spike_protein.fasta", "SARS-CoV_spike_protein_aa.fasta"),
               ("SARS-CoV-2_genome_spike_protein.fasta", "SARS-CoV-2_spike_protein_aa.fasta")]:
    d = "../Alignment/Data/Coronaviruses/"
    check(f"translate-real {nt}", e.translate(e.read_fasta(d + nt)) == e.read_fasta(d + aa))

# ---------- 7. GPU aligner == CPU, and CPU fallback works ----------
print("[7] GPU batched == CPU + fallback")
q = e.read_fasta("../Alignment/Data/Hemoglobin/Homo_sapiens_hemoglobin.fasta")
ts = [e.read_fasta("../Alignment/Data/Hemoglobin/Danio_rerio_hemoglobin.fasta")] * 5
cpu = [e.align_global(q, t, e.matrix_scorer(e.BLOSUM62), 10)[2] for t in ts]
for dev in ("cpu", g.best_device()):
    gs = g.batched_global_scores(q, ts, e.BLOSUM62, 10, g.PROTEIN_ALPHABET, device=dev)
    check(f"gpu batched=={dev}", all(abs(a - b) < 1e-4 for a, b in zip(gs, cpu)))
# forced-bad device must fall back to cpu (not crash)
try:
    gs = g.batched_global_scores(q, ts, e.BLOSUM62, 10, g.PROTEIN_ALPHABET, device="nonsense")
    check("gpu fallback on bad device", all(abs(a - b) < 1e-4 for a, b in zip(gs, cpu)))
except Exception as ex:
    check("gpu fallback on bad device", False)

# ---------- 8. assignment stub functions ----------
print("[8] assignment functions (globalAlignment.py etc.)")
try:
    from globalAlignment import global_alignment
    from localAlignment import local_alignment
    from sharedKMers import count_shared_kmers
    ga = global_alignment("GCATGCU", "GATTACA", 1, 1, 1)
    check("stub global_alignment score", abs(e.alignment_score(ga[0], ga[1], e.dna_scorer(1, 1), 1) - 0.0) < 1e-9)
    la = local_alignment("GGTTGACTA", "TGTTACGG", 1, 1, 1)
    check("stub local_alignment score", abs(e.alignment_score(la[0][0], la[0][1], e.dna_scorer(1, 1), 1) - 4.0) < 1e-9)
    check("stub shared kmers", count_shared_kmers("GATTACA", "GATTACA", 3) == 5)
except Exception as ex:
    check("assignment stubs import/run: " + str(ex), False)

print()
if fails:
    print(f"RESULT: {len(fails)} FAILURES")
    for f in fails:
        print("  -", f)
    sys.exit(1)
else:
    print("RESULT: ALL CHECKS PASSED")
