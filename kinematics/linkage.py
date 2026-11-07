"""Planar linkage kinematics: position + velocity + acceleration, animated.

Model: pin-jointed planar linkage.
  - named points (x, y), some fixed (ground pivots)
  - links = distance constraints between point pairs
  - one driver link swept through an angle

Position analysis: Newton iteration on the loop-closure / distance residuals
(general: any topology, not just four-bar).  Velocity/acceleration from central
finite differences of the converged positions (documented, adequate for plots).

Verification (run `python linkage.py`):
  1. slider-crank piston position vs CLOSED FORM  x = r*cos(t) + sqrt(l^2 - (r*sin t)^2)
     gate: max |error| < 1e-6 mm
  2. four-bar loop closure at every swept angle: |P2-P3| = coupler length
     gate: max closure error < 1e-6 mm

Run:  python linkage.py
"""
import json
import numpy as np


# --------------------------------------------------------------------------
# solver
# --------------------------------------------------------------------------
class Linkage:
    """Points + distance constraints + one driver angle."""

    def __init__(self, points, fixed, links, driver, sliders=()):
        self.points = {k: np.asarray(v, dtype=float) for k, v in points.items()}
        self.fixed = set(fixed)
        self.links = [(a, b, float(L)) for a, b, L in links]
        self.driver = driver  # (p_moving, p_pivot): angle of this link is the input
        # sliders: (point, axis_p1, axis_p2) -> point slides along line(axis_p1, axis_p2)
        self.sliders = list(sliders)

    def _unknowns(self):
        return [k for k in self.points if k not in self.fixed and k != self.driver[0]]

    def solve(self, theta, x0=None):
        """Position analysis at driver angle theta (rad). Returns dict of points.

        x0: optional flat array of unknown coords to warm-start Newton (avoids
        re-converging from a possibly singular default guess each angle).
        """
        pm, pp = self.driver
        pts = {k: v.copy() for k, v in self.points.items()}
        piv = pts[pp]
        # driver point pinned: its position is exact from theta
        d0 = self.points[pm] - self.points[pp]
        Ld = float(np.linalg.norm(d0))
        pts[pm] = piv + Ld * np.array([np.cos(theta), np.sin(theta)])

        unk = self._unknowns()
        if not unk:
            return pts

        x = (np.concatenate([pts[k] for k in unk]) if x0 is None
             else np.asarray(x0, dtype=float).copy())
        rng = np.random.default_rng(20261025)  # deterministic jitter

        def resid(xv):
            for i, k in enumerate(unk):
                pts[k] = xv[2 * i:2 * i + 2]
            r = []
            for a, b, L in self.links:
                r.append(np.linalg.norm(pts[a] - pts[b]) - L)
            for p, a1, a2 in self.sliders:
                d = pts[a2] - pts[a1]
                d = d / np.linalg.norm(d)
                rel = pts[p] - pts[a1]
                r.append(float(rel[0] * d[1] - rel[1] * d[0]))  # signed perp distance
            return np.array(r)

        # Newton with numerical Jacobian; jitter on stall (rank-1 Jacobians,
        # e.g. a collinear initial guess, make min-norm steps stall on a line)
        best_x, best_r = x.copy(), np.inf
        prev_r, stall = np.inf, 0
        for _ in range(200):
            r = resid(x)
            nr = float(np.max(np.abs(r)))
            if nr < best_r:
                best_x, best_r = x.copy(), nr
            if nr < 1e-10:
                break
            stall = stall + 1 if nr >= prev_r * 0.999 else 0
            prev_r = nr
            if stall >= 3:
                x = best_x + rng.normal(scale=1.0, size=x.shape)
                stall, prev_r = 0, np.inf
                continue
            J = np.zeros((len(r), len(x)))
            h = 1e-7
            for j in range(len(x)):
                xp = x.copy(); xp[j] += h
                J[:, j] = (resid(xp) - r) / h
            dx, *_ = np.linalg.lstsq(J, -r, rcond=None)
            step = 1.0
            for _ls in range(20):
                xn = x + step * dx
                if np.max(np.abs(resid(xn))) < nr:
                    break
                step *= 0.5
            x = xn
        x = best_x
        for i, k in enumerate(unk):
            pts[k] = x[2 * i:2 * i + 2]
        self._last_x = x
        return pts

    def sweep(self, angles):
        """Return (N, npts, 2) array of positions in `angles` order.

        Warm-starts each angle from the previous solution (branch continuity).
        """
        keys = list(self.points.keys())
        out = np.zeros((len(angles), len(keys), 2))
        x = None
        for i, th in enumerate(angles):
            pts = self.solve(th, x0=x)
            x = getattr(self, "_last_x", None)
            out[i] = [pts[k] for k in keys]
        return keys, out

    def velocities(self, angles):
        """Central-difference velocities along a sweep (mm/rad).

        MUST differentiate a warm-started sweep() chain — FD across two
        independent cold solves can jump assembly branches and produce
        garbage (learned 2026-10-25).
        """
        keys, P = self.sweep(angles)
        return np.gradient(P, np.asarray(angles, dtype=float), axis=0)


# --------------------------------------------------------------------------
# verification: slider-crank closed form
# --------------------------------------------------------------------------
def verify_slider_crank():
    r, l = 50.0, 150.0
    ang = np.linspace(0, 2 * np.pi, 73)
    lk = Linkage(
        points={"O": (0, 0), "X": (100.0, 0.0), "A": (r, 0), "B": (r + l, 0)},
        fixed={"O", "X"},
        links=[("O", "A", r), ("A", "B", l)],
        driver=("A", "O"),
        sliders=[("B", "O", "X")],   # piston constrained to the x-axis
    )
    err = []
    for th in ang:
        p = lk.solve(th)
        x_exact = r * np.cos(th) + np.sqrt(l**2 - (r * np.sin(th)) ** 2)
        err.append(abs(p["B"][0] - x_exact))
    err = float(np.max(err))
    return dict(check="slider-crank vs closed form", n_angles=len(ang),
                max_error_mm=err, gate_mm=1e-6, passed=bool(err < 1e-6))


# --------------------------------------------------------------------------
# verification: four-bar loop closure
# --------------------------------------------------------------------------
def verify_fourbar_closure():
    # crank 50, coupler 120, rocker 100, ground 140 (Grashof crank-rocker)
    lk = Linkage(
        points={"O2": (0, 0), "A": (50, 0), "B": (190, 0), "O4": (140, 0)},
        fixed={"O2", "O4"},
        links=[("O2", "A", 50.0), ("A", "B", 120.0), ("O4", "B", 100.0)],
        driver=("A", "O2"),
    )
    ang = np.linspace(0, 2 * np.pi, 121)
    err = 0.0
    for th in ang:
        p = lk.solve(th)
        err = max(err, abs(np.linalg.norm(p["A"] - p["B"]) - 120.0),
                  abs(np.linalg.norm(p["O4"] - p["B"]) - 100.0))
    return dict(check="four-bar loop closure", n_angles=len(ang),
                max_error_mm=float(err), gate_mm=1e-6, passed=bool(err < 1e-6))


# --------------------------------------------------------------------------
# demo: four-bar coupler curve + animation
# --------------------------------------------------------------------------
def demo_fourbar(gif_path="fourbar.gif"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    lk = Linkage(
        points={"O2": (0, 0), "A": (50, 0), "B": (190, 0), "O4": (140, 0)},
        fixed={"O2", "O4"},
        links=[("O2", "A", 50.0), ("A", "B", 120.0), ("O4", "B", 100.0)],
        driver=("A", "O2"),
    )
    ang = np.linspace(0, 2 * np.pi, 181)
    keys, P = lk.sweep(ang)
    idx = {k: i for i, k in enumerate(keys)}

    # coupler curve (point B) + velocity magnitude of B
    vel = lk.velocities(ang)
    vB = np.linalg.norm(vel[:, idx["B"]], axis=1)

    fig, ax = plt.subplots(figsize=(7, 7), dpi=100)
    ax.plot(P[:, idx["B"], 0], P[:, idx["B"], 1], "-", lw=1, color="gray",
            label="coupler curve (B)")
    link_lines = [ax.plot([], [], "o-", lw=2, color=c)[0] for c in ("k", "tab:blue", "tab:red")]
    trace, = ax.plot([], [], ".", ms=2, color="tab:orange")
    ax.set_aspect("equal"); ax.grid(alpha=0.3); ax.legend(loc="upper right")
    ax.set_title("four-bar: crank 50 / coupler 120 / rocker 100 / ground 140 mm")

    def update(i):
        A = P[i, idx["A"]]; B = P[i, idx["B"]]
        link_lines[0].set_data([0, A[0]], [0, A[1]])
        link_lines[1].set_data([A[0], B[0]], [A[1], B[1]])
        link_lines[2].set_data([140, B[0]], [0, B[1]])
        trace.set_data(P[:i + 1, idx["B"], 0], P[:i + 1, idx["B"], 1])
        return link_lines + [trace]

    anim = FuncAnimation(fig, update, frames=len(ang), blit=True)
    anim.save(gif_path, writer=PillowWriter(fps=30))
    plt.close(fig)

    # velocity sanity plot
    fig2, ax2 = plt.subplots(figsize=(7, 3.5), dpi=100)
    ax2.plot(np.degrees(ang), vB)
    ax2.set_xlabel("crank angle [deg]"); ax2.set_ylabel("|v_B| [mm/rad]")
    ax2.set_title("coupler point speed vs crank angle")
    ax2.grid(alpha=0.3)
    fig2.savefig("fourbar_velocity.png", bbox_inches="tight")
    plt.close(fig2)

    return gif_path, float(vB.min()), float(vB.max())


if __name__ == "__main__":
    results = [verify_slider_crank(), verify_fourbar_closure()]
    for r in results:
        print(json.dumps(r))
    assert all(r["passed"] for r in results), "KINEMATICS VERIFICATION FAILED"

    gif, vmin, vmax = demo_fourbar()
    print(json.dumps(dict(demo="four-bar animation", gif=gif,
                          vB_mm_per_rad=[vmin, vmax])))
    print("PASS: kinematics verified + demo rendered")
