# Cook's membrane benchmark

This directory contains a small, reproducible geometry generator for the
published linear Cook's membrane benchmark.

Reference case:

- trapezoid vertices: `(0, 0)`, `(48, 44)`, `(48, 60)`, `(0, 44)` mm;
- thickness: `5 mm` in `z`;
- material: `E = 70 MPa`, `nu = 0.33`;
- fixed: left face `x = 0`;
- load: `6.25 N/mm^2` shear traction in `+y` on the right face;
- symmetry: `u_z = 0` on `z = 0` and `z = 5`;
- reported point: `A = (48, 60, 0)` in this solver's coordinate convention.

The published CalculiX very-fine quadratic-tetrahedral reference is
`u_y = 32.26662 mm`.

Source: [CoFEA Cook's membrane benchmark](https://github.com/spolanski/CoFEA/tree/master/benchmarks/02-cooks-membrane).

Generate the local STEP input from the repository root:

```bash
python benchmarks/cook_membrane.py
```

Then run `fea/fea_step.py` with the command in the root README.
