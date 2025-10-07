######################################################################alg.ca
from sage.rings.polynomial.laurent_polynomial import *
from sage.algebras.cluster_algebra import *
from sage.combinat.cluster_algebra_quiver.cluster_seed import *
from quiver import *

mat=matrix
ca=ClusterAlgebra


def etd(self,*args):
    self.explore_to_depth(*args)
ClusterAlgebra.etd = etd

def gv(self,*args):
    return self.g_vector(*args)
ClusterAlgebraSeed.gv = gv
ClusterAlgebra.gv = gv
ClusterSeed.gv = gv


def gvs(self):
    return self.g_vectors()
ClusterAlgebraSeed.gvs = gvs
ClusterAlgebra.gvs = gvs

def gvsf(self):
    return self.g_vectors_so_far()
ClusterAlgebra.gvsf = gvsf

def fp(self,*args):
    return self.F_polynomial(*args)
ClusterAlgebraSeed.fp = fp
ClusterAlgebra.fp = fp


def fps(self):
    return self.F_polynomials()
ClusterAlgebraSeed.fps = fps
ClusterAlgebra.fps = fps

def fpsf(self):
    return self.F_polynomials_so_far()
ClusterAlgebra.fpsf = fpsf


def c_v(self,*args):
    return self.c_vector(*args)
ClusterAlgebraSeed.c_v = c_v


def c_vs(self):
    return self.c_vectors()
ClusterAlgebraSeed.c_vs = c_vs

def cv(self,*args):
    return self.cluster_variable(*args)
ClusterAlgebraSeed.cv = cv
ClusterSeed.cv = cv
ClusterAlgebra.cv = cv

def cvs(self):
    return self.cluster_variables()
ClusterAlgebraSeed.cvs = cvs

def cvs(self):
    return list(self.cluster_variables())
ClusterAlgebra.cvs = cvs

def cvsf(self):
    return self.cluster_variables_so_far()
ClusterAlgebra.cvsf = cvsf

def d_g(self,*args):
    return self.d_vector_to_g_vector(*args)
ClusterAlgebra.d_g = d_g

def g_d(self,*args):
    return self.g_vector_to_d_vector(*args)
ClusterAlgebra.g_d = g_d

def c_sd(self):
    return self.current_seed()
ClusterAlgebra.c_sd = c_sd

def i_sd(self):
    return self.initial_seed()
ClusterAlgebra.i_sd = i_sd

def c_sd(self):
    return self.current_seed()
ClusterAlgebra.c_sd = c_sd

#####################################################################################combinat.ca
from sage.combinat.cluster_algebra_quiver.cluster_seed import *
from sage.combinat.cluster_algebra_quiver.quiver import *

cq=ClusterQuiver
sd = ClusterSeed

def cls(self,*args):
    return self.cluster_class(*args)
ClusterSeed.cls = cls

def cls_it(self,*args):
    return self.cluster_class_iter(*args)
ClusterSeed.cls_it = cls_it

def fp(self,*args):
    return self.f_polynomial(*args)
ClusterSeed.fp = fp

def fps(self):
    return self.f_polynomials()
ClusterSeed.fps = fps

def gm(self):
    return self.g_matrix()
ClusterSeed.gm =gm


########################################

## LL = lower left minor
def F(i, j, g):
    """
    返回以 (i,j) 为左下角的最大方形子矩阵
    参数:
        i, j: 1-based 索引 (左下角位置)
        g: m×n 矩阵
    返回:
        以 (i,j) 为左下角的 k×k 子矩阵，k = min(i, n-j+1)
    """
    m, n = g.nrows(), g.ncols()
    # 计算最大可能的子矩阵大小
    k = min(i, n - j + 1)
    
    # 转换为 0-based 索引
    # 子矩阵的行范围: [i-k, i-1] (0-based)
    start_row = i - k
    end_row = i  # 切片结束索引 (不包含)
    start_col = j - 1
    end_col = j - 1 + k  # 切片结束索引 (不包含)
    
    return g[start_row:end_row, start_col:end_col]

def f(i, j, g):
    """
    返回以 (i,j) 为左下角的最大方形子矩阵的行列式
    参数:
        i, j: 1-based 索引 (左下角位置)
        g: m×n 矩阵
    返回:
        Fij(g) 的行列式（如果子矩阵是方阵）
    """
    sub = F(i, j, g)
    # 确保子矩阵是方阵
    if sub.nrows() != sub.ncols():
        raise ValueError(f"子矩阵不是方阵：{sub.nrows()}×{sub.ncols()}")
    return sub.det().factor()

## UL = upper left minor

def D(i, j, g):
    """
    返回以 (i,j) 为左上角的最大方形子矩阵
    参数:
        i, j: 1-based 索引 (左上角位置)
        g: m×n 矩阵
    返回:
        以 (i,j) 为左上角的 k×k 子矩阵，k = min(m-i+1, n-j+1)
    """
    m, n = g.nrows(), g.ncols()
    # 计算最大可能的子矩阵大小
    k = min(m - i + 1, n - j + 1)
    
    # 转换为 0-based 索引
    # 子矩阵的行范围: [i-1, i-1+k-1] (0-based)
    start_row = i - 1
    end_row = i - 1 + k  # 切片结束索引 (不包含)
    start_col = j - 1
    end_col = j - 1 + k  # 切片结束索引 (不包含)
    
    return g[start_row:end_row, start_col:end_col]

def d(i, j, g):
    """
    返回以 (i,j) 为左上角的最大方形子矩阵的行列式
    参数:
        i, j: 1-based 索引 (左上角位置)
        g: m×n 矩阵
    返回:
        Dij(g) 的行列式（如果子矩阵是方阵）
    """
    sub = D(i, j, g)
    # 确保子矩阵是方阵
    if sub.nrows() != sub.ncols():
        raise ValueError(f"子矩阵不是方阵：{sub.nrows()}×{sub.ncols()}")
    return sub.det().factor()


################################################################################
# ######################### cy.Inv  ############################################
from sage.combinat import *
from itertools import combinations, combinations_with_replacement
# Loop Schur s_{λ/μ}^{(r)} with 1-based residues (r ∈ {1,...,n} mod n)
# Uses sage.combinat (Semistandard[Skew]Tableaux) when available; otherwise falls back to a manual enumerator.

# ---------------- helpers ----------------
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

def _neighbors(cells):
    S = set(cells)
    left_of, up_of = {}, {}
    for (i,j) in cells:
        left_of[(i,j)] = (i, j-1) if (i, j-1) in S else None
        up_of[(i,j)]   = (i-1, j) if (i-1, j) in S else None
    return left_of, up_of




# ---------------- loop Schur s_λ\μ^(r)(x1..xm) ----------------
def sr(lam, mu=None, r=0, m=2, n=3, base_ring=ZZ, varprefix='x', use_combinat=True):
    """
    Loop Schur s_{λ/μ}^{(r)}(x_1,...,x_m) with 1-based colors and c(i,j)=i-j (mod n).

    Weight of a tableau T:
        ∏_{(i,j)∈λ/μ} x_{T(i,j)}^{(a(i,j))},  where  a(i,j) = (r + (i-j)mod n).
    """
    R, X = ring_vars(m, n, base_ring=base_ring, varprefix=varprefix)
    cells, lam_f, mu_f = _cells_of_shape(lam, mu)
    r = Integer(r)

    # try Sage's enumerators
    if use_combinat:
        try:
            from sage.combinat.tableau import SemistandardTableaux, SemistandardSkewTableaux
            iterator = (SemistandardSkewTableaux([lam_f, mu_f], max_entry=m)
                        if (mu is not None and any(mu_f))
                        else SemistandardTableaux(lam_f, max_entry=m))
            total = R(0)
            for T in iterator:
                w = R(1)
                for i_idx, row in enumerate(T, start=1):
                    mu_i = mu_f[i_idx-1] if i_idx-1 < len(mu_f) else 0
                    for local_j, entry in enumerate(row, start=1):
                        j = mu_i + local_j                 
                        a = int(((r + i_idx - j - 1) % n) + 1) 
                        w *= X[(int(entry), a)]
                total += w
            return R, total, X
        except Exception:
            pass

    # fallback: manual SSYT enumeration
    left_of, up_of = _neighbors(cells)
    cells_sorted = sorted(cells)
    assignment = {}
    total = R(0)

    def backtrack(k, weight):
        nonlocal total
        if k == len(cells_sorted):
            total += weight; return
        (i, j) = cells_sorted[k]
        lb = 1
        L = left_of[(i, j)]
        U = up_of[(i, j)]
        if L is not None and L in assignment: lb = max(lb, assignment[L])      # rows weak
        if U is not None and U in assignment: lb = max(lb, assignment[U] + 1)  # cols strict
        a = int(((r + i - j - 1) % n) + 1)
        for val in range(lb, m+1):
            assignment[(i, j)] = val
            backtrack(k+1, weight * X[(val, a)])
            del assignment[(i, j)]

    backtrack(0, R(1))
    return  total

def sr_jt(lam, mu=None, r=0, m=2, n=3, base_ring=ZZ, varprefix='x', X=None):
    """
    Loop Schur s^{(r)}_{λ/μ} via the Jacobi-Trudi determinant:
        s^{(r)}_{λ/μ} = det( h^{(r - μ_j + j - 1)}_{λ_i - μ_j - i + j} )_{1<=i,j<=L}
    where L = max(len(λ), len(μ)).
    Conventions: h_0^{(*)}=1, h_k^{(*)}=0 for k<0. Colors are 1..n (mod n).

    Returns (M, M.det()) where M is the jt_matrix of sr
    """
    
    lam = list(lam)
    mu  = [] if mu is None else list(mu)
    L = max(len(lam), len(mu))
    lam += [0]*(L - len(lam))
    mu  += [0]*(L - len(mu))

    if X is None:
        R, X = ring_vars(m, n, base_ring=base_ring, varprefix=varprefix)
    else:
        R = next(iter(X.values())).parent()

    # Build the L x L matrix of loop-complete functions
    M = matrix(R, L, L)
    for i in range(1, L+1):
        for j in range(1, L+1):
            k  = lam[i-1] - mu[j-1] - i + j               
            sr = r - mu[j-1] + j - 1                      
            if k < 0:
                M[i-1, j-1] = R(0)
            else:
                Hij = hr(k, sr, m=m, n=n, base_ring=base_ring,
                                           varprefix=varprefix, X=X)
                M[i-1, j-1] = Hij

    return M,M.det()

# ---------- loop elementary E_k^(r)(x1..xm) ----------
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

# ---------- loop complete H_k^(r)(x1..xm) ----------
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
