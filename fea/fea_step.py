"""Generic STEP -> FEA: von Mises + deflection report (linear static elasticity).

Pipeline: gmsh (OCC, STEP import) -> linear tets -> scikit-fem ElementTetP2.
Units: mm, N, MPa.  BCs:
  --fix      : all dofs clamped for nodes inside the box (box hugs the fixed face)
  --load-box : boundary facets whose centroid is inside the loaded face box
  --load      : total force in N by default; use --load-mode traction for MPa
                surface traction; --load-area makes total-force normalization exact


Usage:
  python fea_step.py MODEL.step --E 2100 --nu 0.38 \
      --fix=-100,-15,-3:-98,15,3 --load-box=99.99,-15,-3:100.01,15,3 --load 0,0,-15 \
      --probe 100,0,3 --elem-size 4 --yield-mpa 50 --out report.json --png def.png

(load box must hug ONLY the loaded face — a deeper box also captures adjacent
side facets and shifts load onto them.)

Verification anchor: beam_check.py (analytic cantilever; tip deflection
FEM 34.51 mm vs theory 35.30 mm). Cross-check of THIS script on beam.step must
agree with beam_check's FEM value within ~5%.

PITFALLS (learned 2026-10-25):
- scikit-fem 12 ElementHex2 assembly is pathologically slow (56 s for 120
  elements vs 0.5 s for hex1) -> use ElementTetP2.
- von Mises deviator must use the STRESS trace (3*lmbda+2*mu)*tr(eps), NOT the
  strain trace; the wrong trace inflates VM by exactly sqrt(1.5).
- pip installs on this box must use a mirror (pypi.org times out):
  --index-url https://pypi.tuna.tsinghua.edu.cn/simple
"""
import argparse
import json
import time

import numpy as np
import gmsh
import meshio
from skfem import MeshTet, Basis, FacetBasis, LinearForm, asm, condense, solve
from skfem.element import ElementVector, ElementTetP2
from skfem.models.elasticity import linear_elasticity


def box(s):
    lo, hi = s.split(":")
    return (np.array([float(v) for v in lo.split(",")]),
            np.array([float(v) for v in hi.split(",")]))


def in_box_points(x, lo, hi):
    return np.all((x >= (lo - 1e-6)[:, None]) & (x <= (hi + 1e-6)[:, None]), axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step_file")
    ap.add_argument("--E", type=float, required=True, help="Young's modulus, MPa")
    ap.add_argument("--nu", type=float, required=True, help="Poisson ratio")
    ap.add_argument("--fix", type=box, required=True, help="lo_x,lo_y,lo_z:hi_x,hi_y,hi_z")
    ap.add_argument("--load-box", type=box, required=True)
    ap.add_argument("--load", required=True, help="Fx,Fy,Fz in N")
    ap.add_argument("--load-mode", choices=("traction", "total"), default="total",
                    help="surface load interpretation: total force, or constant traction")
    ap.add_argument("--load-area", type=float, default=None,
                    help="loaded face area in mm^2 when --load-mode=total; otherwise estimate from selected facets")
    ap.add_argument("--probe", default=None, help="x,y,z point to report displacement")
    ap.add_argument("--elem-size", type=float, default=4.0, help="gmsh max element size, mm")
    ap.add_argument("--scale", type=float, default=1.0, help="deformed-plot exaggeration")
    ap.add_argument("--yield-mpa", type=float, default=None, help="yield strength for SF")
    ap.add_argument("--symmetry-axis", choices=("x", "y", "z"), default=None,
                    help="constrain displacement normal to a symmetry plane")
    ap.add_argument("--symmetry-at", type=float, default=0.0,
                    help="coordinate of the symmetry plane, mm")
    ap.add_argument("--symmetry-both-faces", action="store_true",
                    help="also constrain the normal displacement on the opposite outer face")
    ap.add_argument("--traction-axis", choices=("x", "y", "z"), default=None,
                    help="optional single axis of the applied surface traction")
    ap.add_argument("--out", default="report.json")
    ap.add_argument("--png", default=None)
    a = ap.parse_args()

    t0 = time.time()
    F = np.array([float(v) for v in a.load.split(",")])

    # ---- mesh via gmsh (OCC STEP import) ------------------------------------
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(a.step_file)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", a.elem_size)
    gmsh.model.mesh.generate(3)
    msh = a.step_file + ".msh"
    gmsh.write(msh)
    gmsh.finalize()

    mesh = meshio.read(msh)
    pts = np.asarray(mesh.points, dtype=float).T
    tets = mesh.cells_dict["tetra"].T.astype(np.int64)
    m = MeshTet(pts, tets)
    print(f"mesh: {tets.shape[1]} tets from {a.step_file} ({time.time()-t0:.1f}s)", flush=True)

    # ---- FE setup -----------------------------------------------------------
    E, nu = a.E, a.nu
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    e = ElementVector(ElementTetP2())
    basis = Basis(m, e)
    K = asm(linear_elasticity(lmbda, mu), basis)
    print(f"assembled K, {basis.N} dofs ({time.time()-t0:.1f}s)", flush=True)

    # ---- load: uniform traction on BOUNDARY facets inside load box ----------
    lo_l, hi_l = a.load_box
    # Boundary facets = facets belonging to exactly one element.
    fcount = np.bincount(m.t2f.reshape(-1), minlength=m.facets.shape[1])
    bnd = np.where(fcount == 1)[0]
    sel = np.intersect1d(m.facets_satisfying(lambda x: in_box_points(x, lo_l, hi_l)), bnd)
    if len(sel) == 0:
        raise SystemExit("load box selects no boundary facets")
    fb = FacetBasis(m, e, facets=sel)

    # A single STEP face may be triangulated into many mesh facets. Assemble
    # the surface integral using a constant traction. For a total force, use
    # the exact planar face area supplied separately by --load-area.
    if a.traction_axis is None:
        nonzero = np.flatnonzero(np.abs(F) > 1e-12)
        if len(nonzero) != 1:
            raise SystemExit("provide --traction-axis for a multi-axis surface load")
        traction_axis = int(nonzero[0])
    else:
        traction_axis = {"x": 0, "y": 1, "z": 2}[a.traction_axis]
    axis_name = "xyz"[traction_axis]
    if a.load_mode == "total":
        if not np.allclose(np.delete(F, traction_axis), 0.0):
            raise SystemExit(
                "--load must be zero in non-selected components; "
                f"got F={F.tolist()} for traction axis {axis_name}"
            )
        if a.load_area is None:
            # Fallback for generic STEP runs: estimate the selected planar
            # face area from its boundary mesh facets. A benchmark can pass
            # an exact CAD face area through --load-area instead.
            fv = m.facets[:, sel]
            q0, q1, q2 = m.p[:, fv[0]], m.p[:, fv[1]], m.p[:, fv[2]]
            area = float(0.5 * np.linalg.norm(
                np.cross(q1 - q0, q2 - q0, axis=0), axis=0
            ).sum())
        else:
            area = a.load_area
        traction_value = F[traction_axis] / area
    else:
        area = None
        traction_value = F[traction_axis]
    tvec = np.zeros(3)
    tvec[traction_axis] = traction_value
    area_text = f"{area:.2f} mm^2" if area is not None else "not normalized"
    print(f"load facets: {len(sel)}, area={area_text}, traction={tvec} MPa, axis={axis_name}, mode={a.load_mode}", flush=True)

    @LinearForm
    def traction(v, w):
        return v[traction_axis] * tvec[traction_axis]

    f = asm(traction, fb)

    # ---- fix: clamp all displacement components on the fixed face -----------
    lo_f, hi_f = a.fix
    fixed = basis.get_dofs(lambda x: in_box_points(x, lo_f, hi_f))
    D = fixed.all()
    print(f"fixed dofs: {len(D)} ({time.time()-t0:.1f}s)", flush=True)

    # Optional symmetry constraints. For a solid model, constrain the
    # displacement normal to a symmetry plane. The selected component includes
    # nodal and edge DOFs for ElementTetP2; facet DOFs may not exist.
    # Syntax: --symmetry-axis z --symmetry-at 0 (axis x/y/z, coordinate mm).
    if a.symmetry_axis is not None:
        axis = {"x": 0, "y": 1, "z": 2}[a.symmetry_axis]
        plane = a.symmetry_at
        planes = [plane]
        if a.symmetry_both_faces:
            planes.append(float(np.max(m.p[axis])))
        component = {0: "u^1", 1: "u^2", 2: "u^3"}[axis]
        symmetry_dofs = []
        for value in planes:
            sym = basis.get_dofs(lambda x, value=value: np.abs(x[axis] - value) <= 1e-6)
            for dof_group in (sym.nodal, sym.edge, sym.facet):
                symmetry_dofs.extend(dof_group.get(component, []))
        D = np.unique(np.concatenate([D, np.asarray(symmetry_dofs, dtype=int)]))
        face_text = ",".join(f"{value:g}" for value in planes)
        print(f"symmetry dofs: {len(D)} total ({a.symmetry_axis}={face_text})", flush=True)

    u = solve(*condense(K, f, D=D))
    # Reaction check: prescribed-DOF reactions should balance the applied load.
    residual = K @ u - f
    reaction = residual[D]
    applied = float(np.sum(f))
    constrained_reaction = float(np.sum(reaction))
    equilibrium_error = (abs(constrained_reaction + applied) /
                         max(abs(applied), 1e-12))
    print(f"equilibrium: applied={applied:.6g}, reactions={constrained_reaction:.6g}, rel_error={equilibrium_error:.3e}", flush=True)
    print(f"solved ({time.time()-t0:.1f}s)", flush=True)

    # ---- post: von Mises (STRESS trace for the hydrostatic part!) ------------
    uh = basis.interpolate(u)
    G = uh.grad
    eps = 0.5 * (G + np.einsum("ij...->ji...", G))
    tr_eps = np.trace(eps, axis1=0, axis2=1)
    sig = lmbda * tr_eps * np.eye(3)[:, :, None, None] + 2.0 * mu * eps
    tr_sig = np.trace(sig, axis1=0, axis2=1)                 # = (3*lmbda+2*mu)*tr_eps
    dev = sig - (tr_sig / 3.0) * np.eye(3)[:, :, None, None]
    vm = np.sqrt(1.5 * np.sum(dev * dev, axis=(0, 1)))       # (nelems,nqp)
    vm_max = float(vm.max())
    vm_elem = vm.max(axis=-1)

    interp = basis.interpolator(u)
    probe = None
    if a.probe:
        p = np.array([float(v) for v in a.probe.split(",")], dtype=float).reshape(3, 1)
        d = interp(p)[:, 0]
        probe = dict(point=a.probe, ux=float(d[0]), uy=float(d[1]), uz=float(d[2]),
                     magnitude=float(np.linalg.norm(d)))

    report = dict(
        model=a.step_file,
        material=dict(E_MPa=E, nu=nu),
        load_N=F.tolist(),
        mesh=dict(elements=int(tets.shape[1]), dofs=int(len(u)), elem_size_mm=a.elem_size),
        max_von_mises_MPa=vm_max,
        p95_von_mises_MPa=float(np.percentile(vm_elem, 95)),
        probe=probe,
        max_abs_displacement_mm=float(np.abs(interp(m.p)).max()),
        equilibrium_relative_error=equilibrium_error,
        safety_factor_vs_yield=(a.yield_mpa / vm_max) if a.yield_mpa else None,
        runtime_s=round(time.time() - t0, 1),
    )
    with open(a.out, "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report, indent=2))

    if a.png:
        plot_results(m, u, vm_elem, a.png, a.scale)


def plot_results(m, u, vm_elem, path, scale=1.0):
    """Deformed boundary surface colored by per-element VM."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    faces = np.vstack([m.t[[0, 1, 2]].T, m.t[[0, 1, 3]].T,
                       m.t[[0, 2, 3]].T, m.t[[1, 2, 3]].T])   # (4*nelems, 3)
    key = np.sort(faces, axis=1)
    uniq, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    bmask = cnt[inv] == 1
    bf = faces[bmask]
    b_elem = np.tile(np.arange(m.t.shape[1]), 4)[bmask]

    basis = Basis(m, ElementVector(ElementTetP2()))
    ip = basis.interpolator(u)
    un = ip(m.p).T                                   # (nvert, 3)
    P = m.p.T + scale * un

    cols = vm_elem[b_elem]
    vmax = max(float(cols.max()), 1e-9)

    fig = plt.figure(figsize=(10, 8), dpi=110)
    ax = fig.add_subplot(111, projection="3d")
    pc = Poly3DCollection(P[bf], cmap="inferno", norm=plt.Normalize(0, vmax))
    pc.set_array(cols)
    ax.add_collection3d(pc)
    for axis, arr in zip("xyz", P.T):
        getattr(ax, f"set_{axis}lim")(arr.min(), arr.max())
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.set_zlabel("z [mm]")
    ax.set_title(f"deformed (x{scale:g}) — VM max {vmax:.2f} MPa")
    fig.colorbar(pc, ax=ax, shrink=0.6, label="VM [MPa]")
    fig.savefig(path, bbox_inches="tight")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
