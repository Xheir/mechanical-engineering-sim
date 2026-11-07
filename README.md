# Mechanical Engineering Simulation Toolkit

Reproducible Python examples for three practical mechanical-engineering tasks:

- **Linear FEA:** STEP geometry → tetrahedral mesh → displacement and von Mises stress
- **Planar kinematics:** slider-crank and four-bar position/velocity analysis
- **Sheet-metal fabrication:** bend allowance → kerf-compensated laser DXF

The examples are intentionally small and inspectable. They are **engineering analysis examples, not certification software**. Results require engineering judgment, appropriate material data, boundary conditions, mesh studies, and experimental correlation before being used for a real part or safety decision.

## Repository layout

```text
fea/
  beam_check.py       # analytic cantilever regression test + FEA cross-check
  fea_step.py         # generic STEP → linear-elastic FEA report
  beam_model.py       # cadgen source for the example STEP beam
  beam.step           # example geometry
  beam_report.json    # recorded example report
  beam_def.png        # recorded deformed-shape visualization

dxf/
  flat_pattern.py             # DXF generator + self-verification
  flat_pattern.dxf            # example output
  flat_pattern_preview.png    # recorded preview

kinematics/
  linkage.py            # generic constrained planar solver + verification/demo
  fourbar.gif           # recorded animation
  fourbar_velocity.png  # recorded velocity plot
```

## Safety and trust boundary

This repository is designed to run **offline on local input files**. The analysis scripts:

- do not open a network listener;
- do not send data to a server or collect telemetry;
- do not read credentials, environment secrets, browser data, or personal files;
- do not execute shell commands or dynamically evaluate code;
- write only explicitly requested solver outputs and plots in the working directory;
- do not control machines, actuators, printers, or other external equipment.

`fea_step.py` reads a STEP path supplied by the operator and writes a mesh beside that input plus the requested JSON/PNG outputs. Run it in a disposable working directory when processing untrusted files, because CAD parsers and native numerical dependencies are third-party software.

These properties are reviewed against the published source, but no public repository can promise that every future dependency or contributor change is harmless. Review changes before running them, keep dependencies updated, and pin the environment for reproducibility.

## Requirements

- Python 3.12 recommended
- `numpy`, `scipy`, `scikit-fem`, `gmsh`, `meshio`, `ezdxf`, `matplotlib`
- For the optional STEP generator: `cadgen` (not required by the FEA solver itself)

Install the pinned Python dependencies in a virtual environment:

```bash
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
# Windows Git Bash
source .venv/Scripts/activate

python -m pip install -r requirements.txt
```

If the default package index is unavailable in your region, use a trusted mirror explicitly. Do not execute downloaded install scripts or commit a private package index credential.

## Verification

Run these checks from the repository root:

```bash
python kinematics/linkage.py
python fea/beam_check.py
python dxf/flat_pattern.py
```

Expected gates:

- slider-crank and four-bar closure error `< 1e-6 mm`;
- cantilever tip-deflection error `< 5%` against the analytic reference;
- midspan von Mises error `< 10%` against the beam reference;
- DXF audit errors `0`;
- DXF flat-length, width, and bend-line errors `< 1e-6 mm`.

The repository's recorded baseline was generated with:

```text
numpy==2.5.4
scipy==1.18.1
scikit-fem==12.0.2
gmsh==4.15.2
meshio==5.3.5
ezdxf==1.4.5
matplotlib==3.11.2
```

The current source includes assertions for the core numerical gates. A passing gate means the example agrees with its stated reference; it does **not** certify a new geometry, material, load case, or manufactured part.

## Generic STEP FEA example

```bash
python fea/fea_step.py path/to/part.step \
  --E 2100 --nu 0.38 \
  --fix=-100,-15,-3:-98,15,3 \
  --load-box=99.99,-15,-3:100.01,15,3 \
  --load 0,0,-15 \
  --probe 100,0,3 \
  --elem-size 4 \
  --yield-mpa 50 \
  --out report.json \
  --png deflection.png
```

Units are mm, N, and MPa. The fixed and loaded boxes must be chosen for the actual geometry. A result from the default beam is not evidence for another part.

## Limitations

- The FEA model is linear elastic and does not cover plasticity, contact, buckling, fatigue, fracture, thermal loads, or nonlinear material behavior.
- The default beam is a regression example, not a validated design.
- The linkage solver uses numerical Jacobians and finite-difference velocities; singular configurations and branch selection require review.
- The DXF generator assumes the stated bend convention, K-factor, and kerf. Fabrication must be checked against the actual machine, material, and process.
- Third-party dependency versions are pinned for the recorded baseline but should still be reviewed for security and compatibility before use.

## License

MIT. See [LICENSE](LICENSE).
