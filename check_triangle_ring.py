"""Exhaustive check of the least energy of the graph H built from a ring of weighted triangles.

For each ring, H is built as in the paper: an inner pole, an outer pole, one vertex per triangle,
a path of length 2A from each triangle vertex to the inner pole and two paths of length 2(B-2A)
to the outer pole. The energy of an orientation is d^T L^+ d with d = out-degree minus in-degree
and L the Laplacian of H, computed with exact rational arithmetic. The least energy over strongly
connected orientations is compared with the closed form W - 4U + 4 rho^2 / S.

Standard library only. Run: python check_triangle_ring.py
"""
import itertools
from fractions import Fraction as F
from math import lcm


def build(pairs):
    n = 2 + len(pairs)          # 0 = inner pole, 1 = outer pole, 2.. = triangle vertices
    edges, paths = [], []

    def path(u, v, length):
        nonlocal n
        seq = [u] + list(range(n, n + length - 1)) + [v]
        n += length - 1
        ids = []
        for a, b in zip(seq, seq[1:]):
            ids.append(len(edges))
            edges.append((a, b))
        paths.append(ids)

    for i, (A, B) in enumerate(pairs):
        path(2 + i, 0, 2 * A)
        path(2 + i, 1, 2 * (B - 2 * A))
        path(2 + i, 1, 2 * (B - 2 * A))
    return n, edges, paths


def reduced_inverse(n, edges):
    L = [[F(0)] * n for _ in range(n)]
    for a, b in edges:
        L[a][a] += 1; L[b][b] += 1; L[a][b] -= 1; L[b][a] -= 1
    m = n - 1
    aug = [L[i + 1][1:] + [F(int(i == j)) for j in range(m)] for i in range(m)]
    for c in range(m):
        p = next(r for r in range(c, m) if aug[r][c] != 0)
        aug[c], aug[p] = aug[p], aug[c]
        piv = aug[c][c]
        aug[c] = [x / piv for x in aug[c]]
        for r in range(m):
            if r != c and aug[r][c] != 0:
                fac = aug[r][c]
                aug[r] = [x - fac * y for x, y in zip(aug[r], aug[c])]
    inv = [row[m:] for row in aug]
    D = lcm(*[x.denominator for row in inv for x in row])
    return [[int(x * D) for x in row] for row in inv], D


def strong(n, edges, signs):
    out = [[] for _ in range(n)]; inn = [[] for _ in range(n)]
    for (a, b), s in zip(edges, signs):
        if s < 0:
            a, b = b, a
        out[a].append(b); inn[b].append(a)

    def reach(adj):
        seen = {0}; st = [0]
        while st:
            u = st.pop()
            for v in adj[u]:
                if v not in seen:
                    seen.add(v); st.append(v)
        return len(seen) == n
    return reach(out) and reach(inn)


def least_energy(pairs, coherent):
    n, edges, paths = build(pairs)
    Inv, Dn = reduced_inverse(n, edges)
    m = len(edges)
    if coherent:
        def gen():
            for ps in itertools.product((1, -1), repeat=len(paths)):
                s = [0] * m
                for pid, sg in zip(paths, ps):
                    for e in pid:
                        s[e] = sg
                yield tuple(s)
        signs = gen()
    else:
        signs = itertools.product((1, -1), repeat=m)
    best_all, best_strong, count = None, None, 0
    for s in signs:
        count += 1
        d = [0] * n
        for (a, b), sg in zip(edges, s):
            if sg > 0:
                d[a] += 1; d[b] -= 1
            else:
                d[b] += 1; d[a] -= 1
        dr = d[1:]
        e = F(sum(dr[i] * sum(Inv[i][j] * dr[j] for j in range(len(dr))) for i in range(len(dr))), Dn)
        if best_all is None or e < best_all:
            best_all = e
        if (best_strong is None or e < best_strong) and strong(n, edges, s):
            best_strong = e
    return m, count, best_all, best_strong


def closed_form(pairs):
    W = sum(2 * A + 4 * (B - 2 * A) for A, B in pairs)
    U = sum(F(B - 2 * A) + F(A * A, B) for A, B in pairs)
    S = sum(F(1, B) for A, B in pairs)
    a = [F(A, B) for A, B in pairs]
    rho = min(abs(t - round(t)) for t in (sum(sg * x for sg, x in zip(sig, a))
                                         for sig in itertools.product((1, -1), repeat=len(a))))
    return W, W - 4 * U + 4 * rho * rho / S, rho


if __name__ == "__main__":
    for pairs, coherent in [([(1, 3)] * 3, False), ([(1, 3), (1, 3), (1, 5)], True)]:
        m, count, emin_all, emin_strong = least_energy(pairs, coherent)
        W, pred, rho = closed_form(pairs)
        print("ring", pairs, "edges", m, "orientations", count,
              "(constant on each path)" if coherent else "(all)")
        scope = "orientations constant on each path" if coherent else "all orientations"
        print("  least energy, %s:" % scope, emin_all)
        print("  least energy, strongly connected among %s:" % scope, emin_strong)
        print("  closed form W - 4U + 4 rho^2/S:", pred, " rho =", rho)
        print("  match:", emin_strong == pred == emin_all)
