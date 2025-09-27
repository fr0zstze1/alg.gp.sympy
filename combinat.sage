######################################################################alg.ca
from sage.rings.polynomial.laurent_polynomial import *
from sage.algebras.cluster_algebra import *
from sage.combinat.cluster_algebra_quiver.cluster_seed import *
from quiver import *

ca=ClusterAlgebra


def gv(self,*args):
    return self.g_vector(*args)
ClusterAlgebraSeed.gv = gv
ClusterAlgebra.gv = gv
ClusterSeed.gv = gv


def gvs(self):
    return self.g_vectors()
ClusterAlgebraSeed.gvs = gvs
ClusterAlgebra.gvs = gvs

def fp(self,*args):
    return self.F_polynomial(*args)
ClusterAlgebraSeed.fp = fp
ClusterAlgebra.fp = fp


def fps(self):
    return self.F_polynomials()
ClusterAlgebraSeed.fps = fps
ClusterAlgebra.fps = fps


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

def cvs(self):
    return self.cluster_variables()
ClusterAlgebraSeed.cvs = cvs

def cvs(self):
    return list(self.cluster_variables())
ClusterAlgebra.cvs = cvs

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

############################################################################# alg.qg
from sage.algebras.quantum_groups.quantum_group_gap import *
from sage.algebras.quantum_groups.q_numbers import *

q_bin = q_binomial
q_fac = q_factorial

qgr=QuantumGroup

def e_(self,i):
    return self.e_tilde(i)
QuaGroupModuleElement.e_ = e_

def f_(self,i):
    return self.f_tilde(i)
QuaGroupModuleElement.f_ = f_

def d(self,k):
    return self^k/q_factorial(k)
QuaGroupModuleElement.d = d

def Tw(self,*args):
    return self.braid_group_action(*args)
QuaGroupModuleElement.Tw = Tw


def hwm(self,wt):
    return self.highest_weight_module(wt)
QuantumGroup.hwm = hwm

def hwv(self):
    return self.highest_weight_vector()
QuantumGroupModule.hwv = hwv

def hw_dec(self):
    return self.highest_weight_decomposition()
QuantumGroupModule.hw_dec = hw_dec

def cyg(self,*args):
    return self.crystal_graph(*args)
QuantumGroupModule.cyg = cyg

def cyb(self,*args):
    return self.crystal_basis(*args)
QuantumGroupModule.cyb = cyb




############################################################################# combinat.sf
def tm(self, basis, *args, **kwargs):
    return self.transition_matrix(basis, *args, **kwargs)

from sage.combinat.sf.macdonald import *
from sage.combinat.sf.sfa import *
from sage.combinat.sf.new_kschur import *


def ex(self,n,**kwargs):
    return self.expand(n,[f'x{i}' for i in [1..n]])
SymmetricFunctionAlgebra_generic_Element.ex = ex

MacdonaldPolynomials_generic.tm = tm
SymmetricFunctionAlgebra_generic.tm=tm
kSchur.tm = tm

kk=FractionField(QQ['q','t'])
kk.inject_variables()

sf=SymmetricFunctions(kk)
sf.inject_shorthands()

HLQp=sf.hall_littlewood().Qp()
Qp1=sf.hall_littlewood(t=1).Qp()
HLP=sf.hall_littlewood().P()
HLQ=sf.hall_littlewood().Q()

P=sf.macdonald().P()
Q=sf.macdonald().Q()
J=sf.macdonald().J()
H=sf.macdonald().H()
Ht=sf.macdonald().Ht()

## H/s_y,x = K(q,t)_x,y
Kqt=lambda x,y:qt_kostka(x,y)

##
## H(q=0) = Qp

#########################################################################combinat.rs-wg
from sage.combinat.root_system.weyl_group import *
from sage.combinat.affine_permutation import *
from sage.combinat.root_system.cartan_type import *

wg=WeylGroup
rs=RootSystem


from sage.combinat.root_system import *
from sage.modules.with_basis.indexed_element import *

ct=CartanType

def cm(self):
    return self.cartan_matrix()
CartanType_abstract.cm = cm
RootSystem.cm = cm

def bm(self):
    return self.cartan_matrix().symmetrized_matrix()
CartanType_abstract.bm = bm
RootSystem.bm = bm
CartanMatrix.bm = bm

def dd(self):
    return DynkinDiagram(self)
CartanType_abstract.dd = dd

def coxeterd(self):
    return self.coxeter_diagram()
CartanType_abstract.coxeterd = coxeterd

def coxeterm(self):
    return self.coxeter_matrix()
CartanType_abstract.coxeterm = coxeterm


def dd_a(self):
    print(self.ascii_art(label=self.a().__getitem__))
CartanType_abstract.dd_a = dd_a

def dd_c(self):
    print(self.ascii_art(label=self.acheck().__getitem__))
CartanType_abstract.dd_c = dd_c

#######################################P, Q, A
#####################rs
from sage.combinat.root_system.ambient_space import *
from sage.combinat.root_system.weight_space import *

def rs_wg(self,*args):
    return WeylGroup(self,*args)
RootSystem.rs_wg = rs_wg

## for A, P, Q, use .weyl_group()

def w_act(self,*args,**kwargs):
    return self.weyl_action(*args,**kwargs)
weight_space.WeightSpaceElement.w_act = w_act
root_space.RootSpaceElement.w_act = w_act
ambinet_space.AmbientSpaceElement.w_act = w_act

def r_sp(self,*args):
    return self.root_space(*args)
RootSystem.r_sp = r_sp

def r_l(self):
    return self.root_lattice()
RootSystem.r_l = r_l

def w_sp(self,*args):
    return self.weight_space(*args)
RootSystem.w_sp = w_sp

def w_l(self):
    return self.weight_lattice()
RootSystem.w_l = w_l

def cr_sp(self,*args):
    return self.coroot_space(*args)
RootSystem.cr_sp = cr_sp
root_space.RootSpace.cr_sp = cr_sp

def cr_l(self):
    return self.coroot_lattice()
RootSystem.cr_l = cr_l
root_space.RootSpace.cr_l = cr_l
ambient_space.AmbientSpace.cr_l = cr_l

def cw_sp(self,*args):
    return self.coweight_space(*args)
RootSystem.cw_sp = cw_sp
root_space.RootSpace.cw_sp = cw_sp

def cw_l(self):
    return self.coweight_lattice()
RootSystem.cw_l = cw_l
root_space.RootSpace.cw_l = cw_l

def a_sp(self,*args):
    return self.ambient_space(*args)
RootSystem.a_sp = a_sp

def a_l(self):
    return self.ambient_lattice()
RootSystem.a_l = a_l

def fw(self):
    return self.fundamental_weights()
weight_space.WeightSpace.fw = fw
ambient_space.AmbientSpace.fw = fw

def al(self):
    return self.alpha()
root_space.RootSpace.al = al
weight_space.WeightSpace.al = al
ambient_space.AmbientSpace.al = al

def alc(self):
    return self.alphacheck()
root_space.RootSpace.alc = alc
weight_space.WeightSpace.alc = alc
ambient_space.AmbientSpace.alc = alc





## class WeylGroup_gens and AffinePermutationGroupGeneric
def s_(self, word):
    return self.from_reduced_word(word)
WeylGroup_gens.s_ = s_

def s_(self, word):
    return self.from_word(word)
AffinePermutationGroupGeneric.s_ = s_

def sref(self):
    return self.simple_reflections()
WeylGroup_gens.sref = sref

def ref(self):
    return self.reflections()
WeylGroup_gens.ref = ref

def chartable(self):  ## character table of finite Weyl gps.
    return self.character_table()
WeylGroup_gens.chartable = chartable

def afgr(self,k):
    return self.affine_grassmannian_elements_of_given_length(k)
WeylGroup_gens.afgr = afgr

def afgr_wd(self,k):
    return [x.reduced_word() for x in self.affine_grassmannian_elements_of_given_length(k)]
WeylGroup_gens.afgr_wd = afgr_wd



## class WeylGroupElement
def wd(self):
    return self.reduced_word()
WeylGroupElement.wd = wd
AffinePermutation.wd = wd

def on(self,v):
    return self.action(v)
WeylGroupElement.on = on

def afgr_cr(self):
    return self.affine_grassmannian_to_core()
WeylGroupElement.afgr_cr = afgr_cr

def afgr_pk(self):
    return self.affine_grassmannian_to_partition()
WeylGroupElement.afgr_pk = afgr_pk


#############################################################################combinat.par-tab
## class Core
from sage.combinat.core import Core
def cr_afgr(self):
    return self.to_grassmannian()
Core.cr_afgr = cr_afgr

def cr_pk(self):
    return self.to_bounded_partition()
Core.cr_pk = cr_pk




from sage.combinat.partition import *

pas=Partitions
pas.options(convention='french')
pa=Partition
skp=SkewPartition
skps=SkewPartitions
cy=crystals
cr=Core
crs=Cores

tx=Tableaux
tab=Tableau

s_tx = StandardTableaux
s_t = StandardTableau
ss_tx = SemistandardTableaux
ss_t = SemistandardTableau

def content_mod(self,*args):
    return self.contents_tableau([Zmod(*args)(0)]).pp()
Partition.content_mod = content_mod

def pk_afgr(self,k):
    return self.from_kbounded_to_grassmannian(k)
Partition.pk_afgr = pk_afgr

def pk_cr_skp(self,k):
    c=self.to_core(k)
    return pa(c).k_skew(k-1)
Partition.pk_cr_skp = pk_cr_skp

def pk_cr(self,k):
    return self.to_core(k)
Partition.pk_cr = pk_cr


def jt(self):
    return self.jacobi_trudi()
Partition.jt = jt


def pk(n,k):
    return Partitions(n,max_part=k)










################################################################################graph

from sage.graphs.graph_plot import *
sage.graphs.graph_plot.DEFAULT_PLOT_OPTIONS['figsize'] = (4)
sage.graphs.graph_plot.DEFAULT_PLOT_OPTIONS['edge_labels']=true

sage.graphs.graph_plot.DEFAULT_PLOT_OPTIONS['layout']='graphviz'
sage.graphs.graph_plot.DEFAULT_PLOT_OPTIONS['color_by_label']=true

sage.graphs.graph_plot.DEFAULT_SHOW_OPTIONS['figsize']= (3)