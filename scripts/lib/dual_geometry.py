"""Streaming reader for the solver's polygon checkpoints (versions 1–6)."""
from contextlib import contextmanager
from dataclasses import dataclass
import math
from pathlib import Path
import struct

import numpy as np


@dataclass(frozen=True)
class GeometryHeader:
    nx: int
    ny: int
    ratio: int
    spacing: float
    depth: float
    epsilon: float
    version: int

    @property
    def fine_cells(self):
        return self.nx * self.ny * self.ratio ** 2


def _header(data):
    if len(data) != 44:
        raise ValueError('Truncated dual geometry header')
    magic, version, nx, ny, ratio, spacing, depth, epsilon = struct.unpack('<5I3d', data)
    if (magic != 0x44554752 or version not in range(1, 7)
            or nx < 1 or ny < 1 or not 1 <= ratio <= 32
            or nx * ny * ratio ** 2 > 10_000_000
            or not all(math.isfinite(v) for v in (spacing, depth, epsilon))
            or spacing <= 0 or depth <= 0 or not 0 <= epsilon < .01):
        raise ValueError('Invalid dual geometry header')
    return GeometryHeader(nx, ny, ratio, spacing, depth, epsilon, version)


def read_geometry_header(path):
    with Path(path).open('rb') as stream:
        return _header(stream.read(44))


@contextmanager
def geometry_records(path, *, include_water=True):
    """Yield a header and an iterator of (fine cell, phase, polygon rings).

    Coordinates remain in the cell's original unit square. Holes and vertex
    order are preserved; callers choose display units and drawing precision.
    The iterator checks the complete footer when exhausted.
    """
    path = Path(path)
    with path.open('rb') as stream:
        size = path.stat().st_size

        def take(count):
            if count > size - stream.tell():
                raise ValueError(f'Truncated fine geometry: {path}')
            data = stream.read(count)
            if len(data) != count:
                raise ValueError(f'Truncated fine geometry: {path}')
            return data

        def skip(count):
            if count > size - stream.tell():
                raise ValueError(f'Truncated fine geometry: {path}')
            stream.seek(count, 1)

        header = _header(take(44))

        def ring(count, keep):
            if count < 3:
                raise ValueError('A polygon ring needs at least three vertices')
            if not keep:
                skip(count * 16)
                return None
            points = np.frombuffer(take(count * 16), dtype='<f8').reshape(-1, 2)
            if not np.isfinite(points).all():
                raise ValueError('Nonfinite polygon coordinates')
            return points

        def records():
            for cell in range(header.fine_cells):
                count, = struct.unpack('<I', take(4))
                if not 1 <= count <= 4096:
                    raise ValueError('Invalid fine-cell piece count')
                for _ in range(count):
                    phase, vertices = struct.unpack('<iI', take(8))
                    if phase not in (-1, 0, 1, 2):
                        raise ValueError(f'Invalid mineral phase {phase}')
                    keep = include_water or phase >= 0
                    rings = [ring(vertices, keep)]
                    holes = struct.unpack('<I', take(4))[0] if header.version >= 2 else 0
                    if holes > (size - stream.tell()) // 4:
                        raise ValueError('Invalid polygon hole count')
                    for _ in range(holes):
                        vertices, = struct.unpack('<I', take(4))
                        points = ring(vertices, keep)
                        if keep:
                            rings.append(points)
                    if keep:
                        yield cell, phase, rings
            cells = header.nx * header.ny
            if header.version >= 3:
                skip(cells * 144)  # Nucleation hazards and event counters.
            if header.version >= 4:
                skip(cells * 72)   # Retained volume and eight-component inventory.
            if header.version >= 5:
                neck, coupling, initialized, count = struct.unpack('<ddII', take(24))
                if (not 0 <= neck <= 1 or not 0 < coupling <= 1
                        or initialized not in (0, 1) or count > header.fine_cells * 4096):
                    raise ValueError('Invalid local-water header')
                skip(count * 76)
            if header.version >= 6:
                skip(24)          # Deferred mineral-volume ledger.
            if stream.read(1):
                raise ValueError(f'Trailing fine geometry data: {path}')

        yield header, records()


def read_geometry(path):
    """Collect rings in physical micrometres for plotting and legacy viewers."""
    polygons, phases, owners = [], [], []
    with geometry_records(path) as (header, records):
        nx, ny, ratio = header.nx, header.ny, header.ratio
        width = nx * ratio
        for cell, phase, rings in records:
            offset = np.array([cell % width, cell // width])
            parent = (cell // width // ratio) * nx + cell % width // ratio
            for ring in rings:
                polygons.append((ring + offset) * header.spacing / ratio * 1e6)
                phases.append(phase)
                owners.append(parent)
    cfg = (nx, ny, ratio, header.spacing, header.depth, header.epsilon)
    return cfg, polygons, np.array(phases), np.array(owners)
