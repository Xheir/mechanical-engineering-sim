import json
import os
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]


def simulation_python() -> str:
    return os.environ.get("MECH_SIM_PYTHON", sys.executable)


def run_script(folder: str, script: str) -> str:
    with TemporaryDirectory(prefix="mechanical-engineering-sim-test-") as workdir:
        result = subprocess.run(
            [simulation_python(), str(ROOT / folder / script)],
            cwd=workdir,
            capture_output=True,
            text=True,
            check=True,
        )
    return result.stdout


def parse_json_objects(output: str) -> list[dict]:
    decoder = json.JSONDecoder()
    objects = []
    cursor = 0
    while (start := output.find("{", cursor)) >= 0:
        try:
            value, consumed = decoder.raw_decode(output[start:])
        except json.JSONDecodeError:
            cursor = start + 1
            continue
        objects.append(value)
        cursor = start + consumed
    return objects


def test_kinematics_verification():
    output = run_script("kinematics", "linkage.py")
    assert "PASS: kinematics verified + demo rendered" in output
    checks = parse_json_objects(output)
    assert all(item["passed"] for item in checks[:2])


def test_fea_regression_verification():
    output = run_script("fea", "beam_check.py")
    assert "PASS: FEA matches beam theory within tolerance" in output
    report = parse_json_objects(output)[0]
    assert report["tip_deflection_mm"]["error_pct"] < 5.0
    assert report["midspan_VM_MPa"]["error_pct"] < 10.0


def test_flat_pattern_verification():
    output = run_script("dxf", "flat_pattern.py")
    assert "PASS: flat pattern verified + preview rendered" in output
    checks = parse_json_objects(output)[1]
    assert checks["audit_errors"] == 0
    assert checks["flat_length_error"] < 1e-6
    assert checks["bend_centerline_error"] < 1e-6
