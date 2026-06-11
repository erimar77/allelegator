"""Pairwise Needleman-Wunsch alignment and center-star multiple sequence alignment."""


def nucleotide_score(a, b):
    return 1 if a == b else -1


# Minimal protein scoring: identity-weighted with a mild mismatch penalty.
# A full BLOSUM62 matrix can replace this later; for near-identical alleles the
# distinction rarely changes the alignment.
def protein_score(a, b):
    return 2 if a == b else -1


def _score_fn(seq_type):
    return nucleotide_score if seq_type == "nucleotide" else protein_score


def needleman_wunsch(seq_a, seq_b, score, gap=-2):
    """Global alignment. Returns (aligned_a, aligned_b, score).

    The gap penalty is deliberately harsher than a mismatch so that
    near-identical, equal-length sequences (alleles) align by substitution
    rather than sprouting spurious gaps; genuine indels still score best.
    """
    n, m = len(seq_a), len(seq_b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + score(seq_a[i - 1], seq_b[j - 1])
            up = dp[i - 1][j] + gap
            left = dp[i][j - 1] + gap
            dp[i][j] = max(diag, up, left)
    ai, bi = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + score(seq_a[i - 1], seq_b[j - 1]):
            ai.append(seq_a[i - 1]); bi.append(seq_b[j - 1]); i -= 1; j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + gap:
            ai.append(seq_a[i - 1]); bi.append("-"); i -= 1
        else:
            ai.append("-"); bi.append(seq_b[j - 1]); j -= 1
    return "".join(reversed(ai)), "".join(reversed(bi)), dp[n][m]


def center_star_msa(seqs, seq_type):
    """Center-star MSA. Returns equal-length aligned rows in input order."""
    if len(seqs) <= 1:
        return list(seqs)
    score = _score_fn(seq_type)

    # Pick the center sequence: highest summed pairwise alignment score.
    best_idx, best_total = 0, None
    for i in range(len(seqs)):
        total = 0
        for j in range(len(seqs)):
            if i == j:
                continue
            _, _, s = needleman_wunsch(seqs[i], seqs[j], score)
            total += s
        if best_total is None or total > best_total:
            best_total, best_idx = total, i

    # Align every other sequence to the center, merging induced gaps so columns
    # stay consistent ("once a gap, always a gap").
    aligned_center = seqs[best_idx]
    others = {}
    for j, s in enumerate(seqs):
        if j == best_idx:
            continue
        ac, aj, _ = needleman_wunsch(aligned_center, s, score)
        others = _propagate_gaps(others, aligned_center, ac)
        aligned_center = ac
        others[j] = aj

    result = [aligned_center if idx == best_idx else others[idx] for idx in range(len(seqs))]
    width = max(len(r) for r in result)
    return [r + "-" * (width - len(r)) for r in result]


def _propagate_gaps(others, old_center, new_center):
    """When the center gains gaps (old_center -> new_center), insert matching
    gaps into previously aligned sequences so columns stay aligned."""
    if not others or old_center == new_center:
        return others
    updated = {k: [] for k in others}
    oi = 0
    for c in new_center:
        if oi < len(old_center) and old_center[oi] == c:
            for k in others:
                updated[k].append(others[k][oi])
            oi += 1
        else:
            for k in others:
                updated[k].append("-")
    return {k: "".join(v) for k, v in updated.items()}
