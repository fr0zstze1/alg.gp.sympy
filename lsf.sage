################################################################################
# ######################### Lsf(Loop symmetric function)  cyinv(crystal invariant)  ############################################

from sage.combinat.partitions import *
from sage.combinat.tableau import *
from itertools import combinations, combinations_with_replacement
mat=matrix
ssytx=SemistandardTableaux
ssyt=SemistandardTableau
ssktx=SemistandardSkewTableaux

# ---------some fast helpers for Lsf ----------------------------------------
def _cells_of_shape(lam, mu=None):
    lam = list(lam)
    mu  = list(mu) if mu is not None else []
    L = max(len(lam), len(mu))
    lam += [0]*(L - len(lam))
    mu  += [0]*(L - len(mu))
    if any(lam[i] < mu[i] for i in range(L)):
        raise ValueError("Invalid skew shape: some λ_i < μ_i.")
    cells = []
    for i in range(1, L+1):
        for j in range(mu[i-1] + 1, lam[i-1] + 1):
            cells.append((i, j))
    return cells, lam, mu

def ring_vars(m, n, base_ring=ZZ, varprefix='x'):
    # Variables named x{i}_{a} (e.g. x1_2) corresponding to x_i^{(a)} with a in {1,...,n}
    names = [f"{varprefix}{i}_{a}" for i in range(1, m+1) for a in range(1, n+1)]
    R = PolynomialRing(base_ring, names)
    gens = R.gens()
    X = {}
    k = 0
    for i in range(1, m+1):
        for a in range(1, n+1):
            X[(i, a)] = gens[k]
            k += 1
    return R, X

def _res1(r, shift, n):
    """Return the 1-based residue ((r + shift) mod n) in {1,...,n}."""
    return int(((Integer(r) - 1 + Integer(shift)) % Integer(n)) + 1)

# inject vars x1_1..xm_n to global env.
def inject(m,n):
    ring_vars(m,n)[0].inject_variables()

def _conjugate_list(parts):
    """Conjugate partition as a Python list."""
    if not parts:
        return []
    return list(Partition(parts).conjugate())




# --------- loop schur function s_λ\μ^(r)(x1..xm) ------------------------------
def sr(lam,mu,r,m,n,base_ring=ZZ,varprefix='x'):
    """
    s_λ\μ^(r)(x1..xm)
    by summing weights over SSYTs of shape λ/μ with entries in {1,...,m}.
    Weight for a tableau T (matrix coords):
        Π_{(i,j)∈λ/μ} x_{T(i,j)}^{( ((r + i - j - 1) % n) + 1 )}.
    """
    R,X=ring_vars(m,n,base_ring=base_ring,varprefix=varprefix)
    
    lam_f = list(lam); mu_f = [] if mu is None else list(mu)
    it = (SemistandardSkewTableaux([lam_f, mu_f], max_entry=m)
        if (mu is not None and any(mu_f))
        else SemistandardTableaux(lam_f, max_entry=m))
    tot=R(0)
    rz=Integer(r)
    
    for T in it:
        w=R(1)
        for i_idx,row in enumerate(T,start=1):
            for j_idx,entry in enumerate(row,start=1):
                if entry is None:
                    continue
                a = int(((rz+i_idx-j_idx-1)%n)+1)
                w *= X[(int(entry),a)]
        tot +=w
    return tot


# --------- loop schur by h-jt matrix ------------------------------------------
def s_jt_h(lam, mu=None, r=0, m=2, n=3):
    """
    Construct the symbolic Jacobi-Trudi matrix for the loop Schur function s^{(r)}_{λ/μ}
    using loop homogeneous symmetric functions H_k^{(r)}, displayed as h{k}_{r}.

    Formula (h-version):
        s^{(r)}_{λ/μ} = det( H^{(r - μ_j + j - 1)}_{λ_i - μ_j - i + j} )_{i,j}
    with H_0^{(*)}=1 and H_k^{(*)}=0 for k<0 or k>m.
    The color index (r part) is reduced mod n into {1,...,n} (no zero).
    """
    r = Integer(r)
    lam = list(lam)
    mu  = [] if mu is None else list(mu)
    L = max(len(lam), len(mu))
    lam += [0]*(L - len(lam))
    mu  += [0]*(L - len(mu))

    M = matrix(SR, L, L)
    for i in range(1, L+1):
        for j in range(1, L+1):
            k  = lam[i-1] - mu[j-1] - i + j
            rshift = r - mu[j-1] + j - 1  # superscript shift

            if k < 0:
                M[i-1, j-1] = SR(0)
                continue
            if k == 0:
                M[i-1, j-1] = SR(1)
                continue
            if k > m:
                M[i-1, j-1] = SR(0)
                continue

            jcolor = int(((Integer(rshift) - 1) % Integer(n)) + 1)
            symname = f"h{k}_{jcolor}"
            M[i-1, j-1] = var(symname)

    return M

# ---------- loop schur by e-jt matrix
def s_jt(lam, mu=None, r=0, m=2, n=3):
    r = Integer(r)
    lam = list(lam)
    mu  = [] if mu is None else list(mu)

    # Conjugates (e-version uses λ', μ')
    lamc = _conjugate_list(lam)
    muc  = _conjugate_list(mu)

    Lp = max(len(lamc), len(muc))
    lamc += [0]*(Lp - len(lamc))
    muc  += [0]*(Lp - len(muc))

    M = matrix(SR, Lp, Lp)

    for i in range(1, Lp+1):
        for j in range(1, Lp+1):
            # Size (row) index for e_k and color (column) index before mod
            k  = lamc[i-1] - muc[j-1] - i + j         # k can be negative/zero/positive
            rshift = r + muc[j-1] - j + 1             # raw color shift (any integer)

            if k < 0:
                M[i-1, j-1] = SR(0)
                continue

            if k == 0:
                # By convention E_0^{(⋅)} = 1
                M[i-1, j-1] = SR(1)
                continue

            # Enforce your constraints:
            #  - vanish when k > m (elementary e_k is 0 for k>m)
            if k > m:
                M[i-1, j-1] = SR(0)
                continue

            # Map color to {1,...,n} (no zero)
            jcolor = int(((Integer(rshift) - 1) % Integer(n)) + 1)

            # Symbol name of the form e{k}_{j}
            symname = f"e{k}_{jcolor}"
            M[i-1, j-1] = var(symname)

    return M

# ---------- loop elementary E_k^(r)(x1..xm) ------------------------------------
def er(k, r, m, n, base_ring=ZZ, varprefix='x', X=None):
    """
    E_k^{(r)}(x_1,...,x_m) with 1-based residue labels in {1,...,n} (mod n):
        product has residues r, r+1, ..., r+k-1.
    Returns (R, polynomial, X)[1].
    If X is provided (a dict (i,a) -> variable), uses its parent ring.
    """
    if X is None:
        R, X = ring_vars(m, n, base_ring=base_ring, varprefix=varprefix)
    else:
        # infer ring from provided variables
        R = next(iter(X.values())).parent()

    if k == 0:
        return R, R(1), X
    if k > m:
        return R, R(0), X

    total = R(0)
    for idxs in combinations(range(1, m+1), k):   # strict i1<...<ik
        term = R(1)
        for t, i in enumerate(idxs):              # t = 0..k-1
            a = _res1(r, t, n)                    # r+t
            term *= X[(i, a)]
        total += term
    return total

# ---------- loop homogeneous H_k^(r)(x1..xm) --------------------------------------------
def hr(k, r, m, n, base_ring=ZZ, varprefix='x', X=None):
    """
    H_k^{(r)}(x_1,...,x_m) with 1-based residue labels in {1,...,n} (mod n):
        product has residues r, r-1, ..., r-k+1.
    Returns (R, polynomial, X)[1].
    If X is provided (a dict (i,a) -> variable), uses its parent ring.
    """
    if X is None:
        R, X = ring_vars(m, n, base_ring=base_ring, varprefix=varprefix)
    else:
        R = next(iter(X.values())).parent()

    if k == 0:
        return R, R(1), X

    total = R(0)
    for idxs in combinations_with_replacement(range(1, m+1), k):  # weak i1<=...<=ik
        term = R(1)
        for t, i in enumerate(idxs):                               # t = 0..k-1
            a = _res1(r, -t, n)                                    # r - t
            term *= X[(i, a)]
        total += term
    return total
