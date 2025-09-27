
from sympy import Matrix, eye, symbols, IndexedBase, diag, Function
import os
from sympy.abc import x,y,z,i,j,t,p,q

mat = Matrix

a=IndexedBase('a')
x=IndexedBase('x')
y=IndexedBase('y')
z=IndexedBase('z')

class R:
    def __getitem__(self, slice_spec):
        if isinstance(slice_spec, slice):
            start = slice_spec.start or 1
            stop = slice_spec.stop
            step = slice_spec.step or 1
            return list(range(start, stop+1, step))
        return list(range(1, slice_spec+1))

r = R()

def b(n,prefix='b'):
    b = IndexedBase(prefix)
    M = Matrix.zeros(n, n) 
    for i in range(n):
        for j in range(i,n):  # 从对角线开始到最后一列
            M[i, j] = b[i+1, j+1]  # 索引从1开始
    return M

def b_(n,prefix='b'):
    b = IndexedBase(prefix)
    M = Matrix.zeros(n,n)
    for i in range(n):
        for j in range(i+1):
            M[i,j] = b[i+1,j+1]
    return M

def u(n,prefix='u'):
    u = IndexedBase(prefix)
    return Matrix(n, n, lambda i,j: 1 if i == j else (u[i+1, j+1] if j > i else 0))

def u_(n,prefix='u'):
    u = IndexedBase(prefix)
    return Matrix(n, n, lambda i,j: 1 if i == j else (u[i+1,j+1] if i > j else 0))

def xij(n,i,j,t):
    I = eye(n)
    E = zeros(n)
    E[i-1,j-1] = 1
    return I + t * E

def sij(n,i,j):
    S = eye(n)
    S[i-1,j-1] = 1
    S[j-1,i-1] = 1

    S[i-1,i-1] = 0
    S[j-1,j-1] = 0

    return S

def n_cell_w(n,w):
    p1=b_(n) * w * b_(n)
    p2=p1.upper_triangular()
    for i in range(p2.rows):
        p2[i,i] = 1
    return p2

def n_sg_w(n,w):
    p1=w**-1 * u_(n) * w
    p2=p1.upper_triangular()
    for i in range(p2.rows):
        p2[i,i] = 1
    return p2




