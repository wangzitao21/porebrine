from pathlib import Path
import numpy as np
from .case_input import result_path

MINERALS = ['Halite', 'Sylvite', 'Carnallite', 'Polyhalite', 'Gypsum', 'Bischofite', 'Calcite']
COMPONENTS = ['Na', 'K', 'Ca', 'Mg', 'Cl', 'S', 'C']
COLUMNS = ['x_m', 'y_m', 'epsilon', 'Ux_m_s', 'Uy_m_s', 'p_Pa'] + MINERALS + ['inert', 'rho_kg_m3'] + [f'm{x}' for x in COMPONENTS] + [f'SI{x}' for x in MINERALS] + ['pH', 'waterActivity', 'qx_s', 'qy_s']
FINAL_COLUMNS = ['x_m', 'y_m', 'epsilon', 'Ux_m_s', 'Uy_m_s', 'p_Pa', 'rho_kg_m3', 'water_equivalent_kg_m3_bulk'] + [f'aqueous_{n}_mol_m3_bulk' for n in COMPONENTS] + [f'mineral_{n}_mol_m3_bulk' for n in MINERALS]


def spatial_order(data):
    scale=max(float(np.max(np.abs(data[:,:2]))),1e-30)
    coordinates=np.rint(data[:,:2]/scale*1e12).astype(np.int64)
    return np.lexsort((coordinates[:,0],coordinates[:,1]))


def partitions(case):
    case = Path(case)
    if (case/'input_data').is_dir():
        case = result_path(case)
    if (case/'data').is_dir():
        return [case/'data']
    return sorted((p for p in case.glob('processor[0-9]*') if p.is_dir() and p.name[9:].isdigit()), key=lambda p:int(p.name[9:])) or [case]


def output_directory(part):
    return part/'output' if (part/'output').is_dir() else part


def frames(case):
    return sorted([p.stem for p in (output_directory(partitions(case)[0]) / 'frames').glob('*.f32')], key=float)


def read_frame(case, time):
    arrays = []
    for part in partitions(case):
        with (output_directory(part) / 'frames' / f'{time}.f32').open('rb') as f:
            header = np.fromfile(f, dtype='<u4', count=4)
            if len(header) != 4 or tuple(header[:2]) != (0x53414C54, 1) or header[3] != len(COLUMNS):
                raise ValueError('Invalid snapshot header')
            values = np.fromfile(f, dtype='<f4')
            if len(values) != header[2]*header[3]:
                raise ValueError('Invalid snapshot length')
            arrays.append(values.reshape(-1, len(COLUMNS)))
    data = np.concatenate(arrays).astype(float)
    data = data[spatial_order(data)]
    if len(np.unique(data[:, :2], axis=0)) != len(data):
        raise ValueError('Duplicate cell coordinates')
    if not np.isfinite(data).all():
        raise ValueError('Nonfinite snapshot value')
    return data


def history(case):
    return np.atleast_1d(np.genfromtxt(output_directory(partitions(case)[0]) / 'history.csv', delimiter=',', names=True))


def final_state(case):
    arrays = []
    for part in partitions(case):
        with (output_directory(part) / 'final_state.f64').open('rb') as f:
            header = np.fromfile(f, dtype='<u4', count=4)
            if len(header) != 4 or tuple(header[:2]) != (0x53414C54, 2) or header[3] != 22:
                raise ValueError('Invalid final-state header')
            a = np.fromfile(f, dtype='<f8')
            if len(a) != header[2]*header[3] or not np.isfinite(a).all():
                raise ValueError('Invalid final-state values')
            a = a.reshape(-1, 22)
            arrays.append(a)
    data = np.concatenate(arrays)
    return data[spatial_order(data)]


def metrics(case):
    import csv
    with (output_directory(partitions(case)[0]) / 'summary.csv').open() as f:
        return {r['metric']: float(r['value']) for r in csv.DictReader(f)}
