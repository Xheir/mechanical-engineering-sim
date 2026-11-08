"""Generate the 3-D Cook's membrane benchmark STEP geometry.

Geometry: trapezoid in x-y, extruded from z=0 to z=5 mm.
Vertices: (0,0), (48,44), (48,60), (0,44) mm.

"""

from pathlib import Path

import gmsh


OUT = Path(__file__).with_name("cook_membrane.step")

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("cook_membrane")

points = [
    gmsh.model.occ.addPoint(0, 0, 0),
    gmsh.model.occ.addPoint(48, 44, 0),
    gmsh.model.occ.addPoint(48, 60, 0),
    gmsh.model.occ.addPoint(0, 44, 0),
]
lines = [
    gmsh.model.occ.addLine(points[0], points[1]),
    gmsh.model.occ.addLine(points[1], points[2]),
    gmsh.model.occ.addLine(points[2], points[3]),
    gmsh.model.occ.addLine(points[3], points[0]),
]
loop = gmsh.model.occ.addCurveLoop(lines)
surface = gmsh.model.occ.addPlaneSurface([loop])
gmsh.model.occ.extrude([(2, surface)], 0, 0, 5)
gmsh.model.occ.synchronize()
gmsh.write(str(OUT))
gmsh.finalize()
print(f"wrote {OUT}")
