"""Thin wrapper the web page calls. It takes plain values from the browser,
runs the real engine, and returns a JSON string. No changes to engine.py."""
import json
import random
import re
import engine as e


def _clean(seq, seqtype):
    seq = seq.upper()
    if seqtype == "dna":
        return re.sub(r"[^ACGTUN]", "", seq)
    return re.sub(r"[^A-Z*]", "", seq)


def _scorer_for(p):
    if p["seqtype"] == "dna":
        return e.dna_scorer(float(p["match"]), float(p["mismatch"]))
    matrix = e.BLOSUM62 if p.get("matrix") == "BLOSUM62" else e.PAM250
    return e.matrix_scorer(matrix)


def _align(seq1, seq2, scorer, gap, mode):
    if mode == "local":
        res = e.align_local(seq1, seq2, scorer, gap)
        return res[0], res[1], res[2]
    if mode == "semiglobal":
        return e.align_semiglobal(seq1, seq2, scorer, gap)
    return e.align_global(seq1, seq2, scorer, gap)


def web_align_json(params_json):
    p = json.loads(params_json)
    seqtype = p["seqtype"]
    seq1 = _clean(p["seq1"], seqtype)
    seq2 = _clean(p["seq2"], seqtype)
    if not seq1 or not seq2:
        return json.dumps({"error": "Enter two sequences to align."})

    mode = p["mode"]
    gapmodel = p["gapmodel"]
    gap = float(p["gap"])

    if seqtype == "dna":
        scorer = e.dna_scorer(float(p["match"]), float(p["mismatch"]))
    else:
        matrix = e.BLOSUM62 if p.get("matrix") == "BLOSUM62" else e.PAM250
        scorer = e.matrix_scorer(matrix)

    # Keep the browser responsive: very large pairs use the banded method.
    if len(seq1) * len(seq2) > 4_000_000:
        row1, row2, score = e.align_banded(seq1, seq2, scorer, gap)
        used = "banded (long sequences)"
    elif gapmodel == "affine" and mode == "global":
        r = e.align_global_affine(seq1, seq2, scorer,
                                  float(p["gap_open"]), float(p["gap_extend"]))
        row1, row2, score = r[0], r[1], r[2]
        used = "global, affine gaps"
    elif mode == "local":
        res = e.align_local(seq1, seq2, scorer, gap)
        row1, row2, score = res[0], res[1], res[2]
        used = "local (Smith-Waterman)"
    elif mode == "semiglobal":
        row1, row2, score = e.align_semiglobal(seq1, seq2, scorer, gap)
        used = "semi-global (free end gaps)"
    else:
        row1, row2, score = e.align_global(seq1, seq2, scorer, gap)
        used = "global (Needleman-Wunsch)"

    ident = e.percent_identity(row1, row2)
    ngaps, glens = e.gap_stats(row1, row2)
    match = "".join(
        "|" if (a == b and a != "-") else (" " if (a == "-" or b == "-") else ".")
        for a, b in zip(row1, row2)
    )
    return json.dumps({
        "row1": row1, "row2": row2, "match": match,
        "score": round(float(score), 2),
        "identity": round(float(ident), 1),
        "gaps": ngaps, "gap_lengths": glens, "used": used,
        "len1": len(seq1), "len2": len(seq2), "cols": len(row1),
    })


def sweep_json(params_json):
    """Align the same pair across a range of gap penalties, so the page can
    show the score climbing as the gap count climbs (the whole point)."""
    p = json.loads(params_json)
    seqtype = p["seqtype"]
    seq1 = _clean(p["seq1"], seqtype)
    seq2 = _clean(p["seq2"], seqtype)
    if not seq1 or not seq2:
        return json.dumps({"error": "Enter two sequences first."})
    mode = p["mode"]
    gmax = 16.0 if seqtype == "dna" else 24.0
    product = len(seq1) * len(seq2)
    steps = 12 if product < 300_000 else 6
    points = []
    for i in range(steps + 1):
        gap = gmax * i / steps
        row1, row2, score = _align(seq1, seq2, _scorer_for(p), gap, mode)
        ng, _ = e.gap_stats(row1, row2)
        points.append({
            "gap": round(gap, 2),
            "score": round(float(score), 2),
            "identity": round(e.percent_identity(row1, row2), 1),
            "gaps": ng,
        })
    return json.dumps({"points": points, "gmax": gmax, "gap_now": float(p["gap"])})


def significance_json(params_json):
    """Align the real pair, then align many shuffles of sequence 2, so the page
    can show whether the real score is better than chance."""
    p = json.loads(params_json)
    seqtype = p["seqtype"]
    seq1 = _clean(p["seq1"], seqtype)
    seq2 = _clean(p["seq2"], seqtype)
    if not seq1 or not seq2:
        return json.dumps({"error": "Enter two sequences first."})
    product = len(seq1) * len(seq2)
    if product > 600_000:
        return json.dumps({"error": "These sequences are long. Try hemoglobin or "
                                    "cytochrome c for the live shuffle test."})
    mode = p["mode"]
    gap = float(p["gap"])
    scorer = _scorer_for(p)
    real = _align(seq1, seq2, scorer, gap, mode)[2]
    n = 300 if product < 60_000 else (150 if product < 300_000 else 60)
    rnd = random.Random(0)
    chars = list(seq2)
    nulls = []
    for _ in range(n):
        rnd.shuffle(chars)
        nulls.append(_align(seq1, "".join(chars), scorer, gap, mode)[2])
    mu = sum(nulls) / len(nulls)
    sd = (sum((x - mu) ** 2 for x in nulls) / len(nulls)) ** 0.5
    z = (real - mu) / sd if sd > 0 else 0.0
    return json.dumps({
        "nulls": [round(float(x), 2) for x in nulls],
        "real": round(float(real), 2), "mu": round(mu, 2),
        "sd": round(sd, 2), "z": round(z, 1), "n": n,
    })


def trace_align_json(params_json):
    """Run the alignment under a line tracer so the page can show exactly which
    lines of engine.py executed. Uses a short slice of the sequences so the
    trace is fast; the code path is the same as the full alignment."""
    import sys
    p = json.loads(params_json)
    seqtype = p["seqtype"]
    seq1 = _clean(p["seq1"], seqtype)[:70]
    seq2 = _clean(p["seq2"], seqtype)[:70]
    if not seq1 or not seq2:
        return json.dumps({"error": "Enter two sequences first."})
    scorer = _scorer_for(p)
    gap = float(p["gap"])
    mode = p["mode"]
    gapmodel = p["gapmodel"]
    engine_file = e.__file__

    seen = set()
    order = []

    def tracer(frame, event, arg):
        if event == "line" and frame.f_code.co_filename == engine_file:
            ln = frame.f_lineno
            if ln not in seen:
                seen.add(ln)
                order.append(ln)
        return tracer

    sys.settrace(tracer)
    try:
        if gapmodel == "affine" and mode == "global":
            entry = "align_global_affine"
            e.align_global_affine(seq1, seq2, scorer, float(p["gap_open"]), float(p["gap_extend"]))
        else:
            entry = {"local": "align_local", "semiglobal": "align_semiglobal"}.get(mode, "align_global")
            _align(seq1, seq2, scorer, gap, mode)
    finally:
        sys.settrace(None)

    return json.dumps({"lines": sorted(seen), "order": order, "entry": entry,
                       "total": len(seen)})


def do_translate(payload):
    # payload is "seq1||seq2"; translate both, return a JSON list.
    parts = payload.split("||")
    outs = [e.translate(re.sub(r"[^ACGTUNacgtun]", "", x)) for x in parts]
    return json.dumps(outs)


def do_revcomp(dna):
    return e.reverse_complement(re.sub(r"[^ACGTUNacgtun]", "", dna))
