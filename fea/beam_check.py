"""FEA smoke test: cantilever beam (tet P2) vs Euler-Bernoulli / Timoshenko.

PITFALL: ElementHex2 assembly in scikit-fem 12 is pathologically slow
(56 s for 120 elements vs 0.5 s for hex1) -> use ElementTetP2.

Geometry: L x b x h, clamped at x=0, uniform traction on x=L face = -P/A in z.
Verification targets:
  tip deflection  d_EB   = P*L^3 / (3*E*I),  I = b*h^3/12
  shear add-on    d_sh   = (6/5) * P*L / (G*A),  G = E/(2*(1+nu))
  midspan stress  sig(x) = P*(L-x)*(h/2)/I   (Saint-Venant, away from root)

Run:  python beam_check.py
"""
import json
import numpy as np

from skfem import MeshTet, Basis, FacetBasis, LinearForm, asm, condense, solve
from skfem.element import ElementVector, ElementTetP2
from skfem.models.elasticity import linear_elasticity

# ---- parameters (PETG scraper-arm-like) ------------------------------------
L, b, h = 200.0, 30.0, 6.0      # mm
P = 15.0                        # N, -z on tip face
E, nu = 2100.0, 0.38            # MPa (PETG)

# ---- theory ----------------------------------------------------------------
I = b * h**3 / 12.0
A = b * h
G = E / (2 * (1 + nu))
d_eb = P * L**3 / (3 * E * I)
d_sh = (6.0 / 5.0) * P * L / (G * A)
d_theory = d_eb + d_sh
sig_mid_theory = P * (L - 0.25 * L) * (h / 2) / I

# ---- mesh + basis ----------------------------------------------------------
NX, NY, NZ = 48, 8, 4
m = MeshTet.init_tensor(
    np.linspace(0, L, NX + 1), np.linspace(0, b, NY + 1), np.linspace(0, h, NZ + 1)
)
e = ElementVector(ElementTetP2())
basis = Basis(m, e)

lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
mu = G

K = asm(linear_elasticity(lmbda, mu), basis)

A_t = A
@LinearForm
def tip_traction(v, w):
    return v[2] * (-P / A_t)

fbasis = FacetBasis(m, e, facets=m.facets_satisfying(lambda x: x[0] > L - 1e-6))
f = asm(tip_traction, fbasis)

D = basis.get_dofs(lambda x: x[0] < 1e-6).all()
u = solve(*condense(K, f, D=D))

# ---- post: von Mises at quadrature points ----------------------------------
uh = basis.interpolate(u)
Gm = uh.grad                                   # (3,3,nelems,nqp): [i,j]=du_i/dx_j
eps = 0.5 * (Gm + np.einsum('ij...->ji...', Gm))
tr = np.trace(eps, axis1=0, axis2=1)
sig = lmbda * tr * np.eye(3)[:, :, None, None] + 2.0 * mu * eps
s_tr = np.trace(sig, axis1=0, axis2=1)         # STRESS trace = (3*lmbda+2*mu)*tr_eps
dev = sig - (s_tr / 3.0) * np.eye(3)[:, :, None, None]
vm = np.sqrt(1.5 * np.sum(dev * dev, axis=(0, 1)))   # (nelems,nqp)
vm_max = float(vm.max())
vm_elem = vm.max(axis=-1)

# midspan VM: elements with centroid x in [0.24L, 0.26L]
cx = m.p[0][m.t].mean(axis=0)
band = np.abs(cx - 0.25 * L) < 0.02 * L
vm_mid = float(vm_elem[band].max()) if band.any() else float('nan')

# ---- tip deflection via interpolator --------------------------------------
interp = basis.interpolator(u)
tip_pts = np.array([
    [L, L, L, L], [0, b, 0, b], [0, 0, h, h]
], dtype=float)
u_tip = interp(tip_pts)                        # (3, 4)
d_fem = float(-u_tip[2].mean())                # +z load is -z; report positive drop

err = abs(d_fem - d_theory) / d_theory * 100
err_s = abs(vm_mid - sig_mid_theory) / sig_mid_theory * 100 if band.any() else float('nan')

report = dict(
    geometry_mm=dict(L=L, b=b, h=h), load_N=P, material=dict(E_MPa=E, nu=nu),
    tip_deflection_mm=dict(fem=d_fem, theory=d_theory,
                           eber=d_eb, shear_add=d_sh, error_pct=err),
    midspan_VM_MPa=dict(fem=vm_mid, theory=sig_mid_theory, error_pct=err_s),
    max_VM_MPa=vm_max,
    root_VM_theory_MPa=P * L * (h / 2) / I,
    safety_factor_vs_yield50=50.0 / vm_max,
    mesh=dict(elements=int(m.t.shape[1]), dofs=int(len(u))),
)

print(json.dumps(report, indent=2))

# quick gate: FEA must match theory
assert err < 5.0, f"tip deflection off by {err:.2f}%"
assert err_s < 10.0, f"midspan VM off by {err_s:.2f}%"
print("PASS: FEA matches beam theory within tolerance")
