"""Sheet-metal flat pattern -> laser DXF (kerf-compensated), with self-verification.

Geometry convention (90-deg-friendly, works for any bend angle):
  A_out, B_out : OUTSIDE flange lengths (mm), measured to the outside bend surface
  W            : part width (mm)
  t            : thickness (mm)
  R            : inside bend radius (mm)
  theta        : bend angle (deg)
  K            : K-factor (neutral axis offset / t, typical 0.33-0.50)
  kerf         : laser kerf width (mm); cut path is offset OUTWARD by kerf/2 so
                 the cooled part comes out at nominal size

Formulas:
  BA   = theta_rad * (R + K*t)               (bend allowance)
  OSSB = (R + t) * tan(theta/2)              (outside setback)
  flat = A_out + B_out - 2*OSSB + BA

DXF layers:
  CUT   : kerf-compensated closed outline (LWPOLYLINE) — the laser path
  BEND  : bend centerline + tangent lines (dashed)
  ETCH  : labels/dims (engrave or ignore)

Units: mm ($INSUNITS=4).  Verification gate (run `python flat_pattern.py`):
  ezdxf audit 0 errors AND re-measured DXF geometry reproduces the analytic
  flat length / bend-line position within 1e-6 mm.

Run:  python flat_pattern.py
"""
import json
import math

import ezdxf

# ---- parameters -----------------------------------------------------------
A_OUT, B_OUT = 100.0, 50.0   # outside flange lengths, mm
WIDTH = 60.0
T = 2.0                      # thickness
R = 3.0                      # inside bend radius
THETA = 90.0                 # bend angle, deg
K = 0.44                     # K-factor (soft steel ~0.44, PETG/alu varies)
KERF = 0.15                  # laser kerf, mm
OUT = "flat_pattern.dxf"


def bend_math(a_out, b_out, t, r, theta_deg, k):
    th = math.radians(theta_deg)
    ba = th * (r + k * t)                 # bend allowance
    ossb = (r + t) * math.tan(th / 2)     # outside setback
    flat = a_out + b_out - 2 * ossb + ba
    return dict(bend_allowance=ba, outside_setback=ossb, flat_length=flat,
                bend_line_from_edge=a_out - ossb)


def build_dxf(path):
    bm = bend_math(A_OUT, B_OUT, T, R, THETA, K)
    flat = bm["flat_length"]
    off = KERF / 2.0                      # kerf compensation, outward

    doc = ezdxf.new("R2010", setup=True)
    doc.header["$INSUNITS"] = 4           # mm
    msp = doc.modelspace()

    doc.layers.add("CUT", color=1)
    doc.layers.add("BEND", color=5, linetype="DASHED")
    doc.layers.add("ETCH", color=3)

    # CUT: kerf-compensated rectangle (the laser path)
    x0, y0 = -off, -off
    x1, y1 = flat + off, WIDTH + off
    msp.add_lwpolyline(
        [(x0, y0), (x1, y0), (x1, y1), (x0, y1)],
        close=True, dxfattribs={"layer": "CUT"},
    )

    # BEND: centerline (dashed) + tangent lines
    xb = bm["bend_line_from_edge"]
    msp.add_line((xb, y0), (xb, y1), dxfattribs={"layer": "BEND"})
    for xt in (xb - bm["outside_setback"], xb + bm["outside_setback"]):
        msp.add_line((xt, y0), (xt, y1), dxfattribs={"layer": "BEND"})

    # ETCH: labels
    msp.add_text(
        f"L-bracket flat  A={A_OUT:g} B={B_OUT:g} W={WIDTH:g} t={T:g} "
        f"R={R:g} {THETA:g}deg K={K:g}",
        height=3.0, dxfattribs={"layer": "ETCH"},
    ).set_placement((0, WIDTH + 4))
    msp.add_text(
        f"flat={flat:.3f}  BA={bm['bend_allowance']:.3f}  "
        f"OSSB={bm['outside_setback']:.3f}  kerf={KERF:g}",
        height=3.0, dxfattribs={"layer": "ETCH"},
    ).set_placement((0, WIDTH + 9))

    doc.saveas(path)
    return bm


def verify(path):
    doc = ezdxf.readfile(path)
    auditor = doc.audit()
    msp = doc.modelspace()

    cut = [e for e in msp if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == "CUT"]
    assert len(cut) == 1, f"expected 1 CUT polyline, got {len(cut)}"
    pl = cut[0]
    pts = [(p[0], p[1]) for p in pl.get_points()]
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    w_out = max(xs) - min(xs) - KERF       # de-compensated part width
    h_out = max(ys) - min(ys) - KERF

    bm = bend_math(A_OUT, B_OUT, T, R, THETA, K)

    bends = [e for e in msp if e.dxftype() == "LINE" and e.dxf.layer == "BEND"]
    bend_xs = sorted({round(e.dxf.start.x, 6) for e in bends})
    # centerline is the middle x of the 3 bend lines (center ± setback)
    centerline = bend_xs[len(bend_xs) // 2] if len(bend_xs) == 3 else None

    checks = dict(
        audit_errors=len(auditor.errors),
        polyline_closed=bool(pl.closed),
        flat_length_from_dxf=w_out,
        flat_length_formula=bm["flat_length"],
        flat_length_error=abs(w_out - bm["flat_length"]),
        width_from_dxf=h_out,
        width_error=abs(h_out - WIDTH),
        bend_centerline_from_dxf=centerline,
        bend_centerline_formula=bm["bend_line_from_edge"],
        bend_centerline_error=(abs(centerline - bm["bend_line_from_edge"])
                               if centerline is not None else None),
    )
    return checks


if __name__ == "__main__":
    bm = build_dxf(OUT)
    print(json.dumps(dict(formula=bm), indent=2))
    checks = verify(OUT)
    print(json.dumps(checks, indent=2))

    ok = (checks["audit_errors"] == 0 and checks["polyline_closed"]
          and checks["flat_length_error"] < 1e-6 and checks["width_error"] < 1e-6
          and checks["bend_centerline_error"] is not None
          and checks["bend_centerline_error"] < 1e-6)
    assert ok, "FLAT PATTERN VERIFICATION FAILED"

    # preview render (ezdxf drawing add-on -> matplotlib)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ezdxf.addons.drawing import RenderContext, Frontend
    from ezdxf.addons.drawing.matplotlib import MatplotlibBackend

    doc = ezdxf.readfile(OUT)
    fig = plt.figure(figsize=(10, 5), dpi=110)
    ax = fig.add_axes([0, 0, 1, 1])
    Frontend(RenderContext(doc), MatplotlibBackend(ax)).draw_layout(doc.modelspace())
    ax.set_aspect("equal")
    fig.savefig("flat_pattern_preview.png", bbox_inches="tight")
    print("wrote flat_pattern_preview.png")
    print("PASS: flat pattern verified + preview rendered")
