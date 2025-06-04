import re
from pathlib import Path
from typing import Optional, Tuple
import numpy as np


def _extract_dipole_components(lines):
    pattern = re.compile(r"Dipole moment .*?X=\s*([-\d\.E+e]+).*?Y=\s*([-\d\.E+e]+).*?Z=\s*([-\d\.E+e]+)")
    matches = pattern.findall("".join(lines))
    if not matches:
        return None
    dx, dy, dz = map(float, matches[-1])
    return np.array([dx, dy, dz])


def _extract_polar_matrix(lines):
    header_idx = None
    for i, line in enumerate(lines):
        if "Polarizability" in line and ("SCF" in line or "freq" in line):
            header_idx = i
    if header_idx is None:
        return None
    matrix = []
    for l in lines[header_idx + 1 : header_idx + 4]:
        nums = re.findall(r"[-\d\.E+e]+", l)
        if len(nums) >= 3:
            matrix.append([float(n) for n in nums[:3]])
    if len(matrix) == 3:
        return np.array(matrix)
    return None


def analyze_polarizability(log_file: str) -> Optional[Tuple[float, float]]:
    path = Path(log_file)
    if not path.exists():
        return None
    lines = path.read_text().splitlines()
    dipole = _extract_dipole_components(lines)
    polar = _extract_polar_matrix(lines)
    if dipole is None or polar is None:
        return None
    norm = np.linalg.norm(dipole)
    if norm == 0:
        return None
    p = dipole / norm
    axial = float(p @ polar @ p)
    isotropic = float(np.trace(polar))
    perpendicular = (isotropic - axial) / 2
    return axial, perpendicular


__all__ = ["analyze_polarizability"]
