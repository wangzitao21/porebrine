"""Shared case definitions and small I/O helpers for the mineral studies."""
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import ROOT

STUDY = ROOT / 'cases/carnallite_coating/output/dual_resolution'
DUAL_MINERALS = ('Halite', 'Sylvite', 'Carnallite')
DEFAULT_MODES = ('carnallite_only', 'unblocked', 'coupled')
LABELS = {'carnallite_only': 'Carnallite only', 'frozen_seeds': 'Fixed seeds',
          'unblocked': 'Reactions without shielding', 'coupled': 'Reactions with shielding',
          'mixed': 'Three initial minerals', 'coated': 'Three seeded minerals'}
COMPARISON_SETTINGS = {'end': 2., 'dt': .01, 'inlet_flow': 8e-15, 'nucleation_prefactor': 1e10}


@dataclass(frozen=True)
class StudyCase:
    mode: str
    tag: str
    rate_scale: float
    title: str
    description: str

    @property
    def name(self):
        return self.mode + self.tag


COMPARISON_CASES = (
    StudyCase('coupled', '_poisson_rate_1', 1., '溶蚀与新生沉淀',
              '饱和石盐水流过光卤石，产生石盐和钾石盐沉淀；生长位置、接触面积和附近水体随时间更新。'),
    StudyCase('coupled', '_poisson_rate_10', 10., '盐类反应加快',
              '初始矿物、入口卤水和成核参数相同。石盐与钾石盐反应系数提高 10 倍，观察覆盖与溶液变化的竞争。'),
    StudyCase('carnallite_only', '_poisson_control', 1., '单矿物对照',
              '相同初始光卤石和入口卤水，仅启用光卤石反应，用于对比沉淀带来的影响。'),
)


def save_rows(path, rows):
    rows = iter(rows)
    first = next(rows, None)
    if first is None:
        return
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(first))
        writer.writeheader()
        writer.writerow(first)
        writer.writerows(rows)


def flatten(value, prefix=''):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from flatten(item, f'{prefix}.{key}' if prefix else key)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from flatten(item, f'{prefix}[{index}]')
    else:
        yield {'parameter': prefix, 'value': value}


def saturation(chemistry):
    lo, hi = 4., 8.
    for _ in range(56):
        mid = (lo + hi) / 2
        if chemistry.fractions([mid, 0, 0, 0, mid, 0, 0])[1].si[0] > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def reactive_water(run, stamp, mineral='Carnallite'):
    path = Path(run) / f'foam/output/dual_geometry/{stamp}.reactors.csv'
    if not path.is_file():
        return None
    values = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
    if values.size == 0:
        return {'mK': 0., 'affinity': 0.}
    area = values['A_' + mineral]
    weights = area / area.sum() if area.sum() else area
    return {'mK': float(weights @ values['mK']),
            'affinity': float(weights @ (1 - 10 ** values['SI' + mineral]))}
