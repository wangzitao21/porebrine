"""Export recorded mineral geometry and transport fields to one offline HTML."""
from pathlib import Path
import argparse
import base64
import csv
import json
import sys
import zlib

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from lib.fields import COLUMNS, frames, history, metrics, read_frame
from lib.dual_geometry import geometry_records
from lib.study import COMPARISON_CASES, DUAL_MINERALS as MINERALS, STUDY

DISPLAY_UNIT_UM = 1e-3


def table(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def compact(values):
    return [float(f'{x:.7g}') if np.isfinite(x) else None for x in values]


def display_ring(points):
    """Quantise only drawing coordinates, then drop duplicate/collinear pixels.

    The coordinate unit stays below 0.1 screen pixel at the viewer's maximum
    zoom. Physical volumes and all validation use the original binary rings.
    """
    q = np.rint(points / DISPLAY_UNIT_UM).astype(np.int64)
    q = q[np.any(q != np.roll(q, 1, axis=0), axis=1)]
    while len(q) >= 3:
        before = q - np.roll(q, 1, axis=0)
        after = np.roll(q, -1, axis=0) - q
        keep = ((before[:, 0] * after[:, 1] - before[:, 1] * after[:, 0]) != 0) | (np.sum(before * after, axis=1) < 0)
        if keep.all():
            return np.round(q * DISPLAY_UNIT_UM, 4).ravel().tolist()
        q = q[keep]
    return []


def geometry(path):
    """Keep compound polygon holes; validate its volume against field output."""
    polygons = [[], [], []]
    with geometry_records(path, include_water=False) as (header, records):
        nx, ny, ratio = header.nx, header.ny, header.ratio
        H, depth, epsilon = header.spacing, header.depth, header.epsilon
        fractions = np.zeros((nx * ny, 3))
        scale = H / ratio * 1e6
        width = nx * ratio
        for cell, phase, rings in records:
            offset = np.array([cell % width, cell // width])
            parent = (cell // width // ratio) * nx + cell % width // ratio
            areas, output = [], []
            for ring in rings:
                following = np.roll(ring, -1, axis=0)
                areas.append(abs(float(np.sum(ring[:, 0] * following[:, 1] - ring[:, 1] * following[:, 0]))) / 2)
                output.append(display_ring((ring + offset) * scale))
            area = areas[0] - sum(areas[1:])
            if area < -1e-12:
                raise ValueError('Polygon holes exceed exterior')
            fractions[parent, phase] += area * (1 - epsilon) / ratio ** 2
            if output[0]:
                polygons[phase].append([output[0], *(r for r in output[1:] if r)])
    return {'nx': nx, 'ny': ny, 'ratio': ratio, 'H': H * 1e6,
            'h': scale, 'depth': depth * 1e6, 'width': nx * H * 1e6,
            'height': ny * H * 1e6}, polygons, fractions


def build_run(study, name, title, description, until=None):
    run = study / name
    output = run / 'foam/output'
    hist = history(run / 'foam')
    prefix_path = output / 'history_prefix.csv'
    if prefix_path.exists():
        prefix = np.atleast_1d(np.genfromtxt(prefix_path, delimiter=',', names=True))
        if prefix.dtype.names != hist.dtype.names:
            raise ValueError('Restart history schema mismatch')
        overlap = np.flatnonzero(abs(prefix['time_s'] - hist['time_s'][0]) < 1e-10)
        if len(overlap) != 1:
            raise ValueError('Missing restart boundary sample')
        scale = max(float(prefix['min_Carnallite'][0]), 1e-30)
        for mineral in MINERALS:
            for key in ('min_', 'prec_', 'diss_'):
                if abs(prefix[key + mineral][overlap[0]] - hist[key + mineral][0]) > 1e-8 * scale:
                    raise ValueError('Restart changed the mineral history')
        hist = np.concatenate([prefix[prefix['time_s'] < hist['time_s'][0] - 1e-10], hist])
    parameters = {r['parameter']: r['value'] for r in table(run / 'parameters.csv')}
    stamps = frames(run / 'foam')
    if until is not None:
        stamps = [stamp for stamp in stamps if float(stamp) <= until + 1e-10]
        hist = hist[hist['time_s'] <= until + 1e-10]
    times = [float(t) for t in stamps]
    expected_end = float(parameters['simulation.end_time_s'])
    requested_end = expected_end if until is None else until
    if not times or abs(times[-1] - requested_end) > 1e-10:
        raise ValueError(f'Incomplete result sequence: {name}')
    complete = abs(times[-1] - expected_end) <= 1e-10 and (output / 'summary.csv').exists()
    if until is None:
        metric = metrics(run / 'foam')
    else:
        # A progress preview reports only the actually saved records.
        metric = {'physical_time_s': times[-1], 'max_component_balance_error':
                  max(float(np.max(np.abs(hist[key]))) for key in hist.dtype.names if key.startswith('balance_'))}
    metric['max_component_balance_error'] = max(metric['max_component_balance_error'],
        max(float(np.max(np.abs(hist[key]))) for key in hist.dtype.names if key.startswith('balance_')))
    if metric['max_component_balance_error'] > 1e-7:
        raise ValueError(f'Conservation check failed: {name}')
    if int(parameters.get('geometry.initial_seed_count', '-1')) != 0:
        raise ValueError(f'This viewer requires zero pre-existing precipitate seeds: {name}')
    if name.startswith('coupled_') and parameters.get('geometry.nucleation_process') != 'area-weighted hazard with reproducible exponential event thresholds':
        raise ValueError(f'This viewer requires the corrected nucleation process: {name}')
    if any(abs(hist['min_' + m][0]) > 1e-30 for m in ('Halite', 'Sylvite')):
        raise ValueError(f'This viewer starts with carnallite only: {name}')
    record = {'id': name, 'title': title, 'description': description,
              'metrics': metric, 'times': times, 'frames': [], 'curve': [],
              'complete': complete, 'target_end_time_s': expected_end,
              'parameters': {k: parameters[k] for k in (
                  'simulation.end_time_s', 'physics.temperature_C',
                  'initial_solution_mol_kg_water.Na', 'boundaries[0].value_m3_s')},
              'hashes': table(run / 'program_hashes.csv')}
    record['parameters']['initial_seed_count'] = 0
    record['parameters']['nucleation_prefactor_m2_s'] = float(parameters.get('reactions.nucleation[0].prefactor_m2_s', 0))
    record['parameters']['kinetic_coefficients_mol_m2_s'] = [float(parameters['reactions.rate_coefficients_mol_m2_s.' + m]) for m in MINERALS]
    if (run / 'restart.csv').exists():
        record['restart'] = {r['property']: r['value'] for r in table(run / 'restart.csv')}
        record['prefix_hashes'] = table(run / 'prefix_program_hashes.csv')
    max_mapping = 0.
    for stamp, time in zip(stamps, times):
        cfg, polygons, fractions = geometry(output / f'dual_geometry/{stamp}.bin')
        if 'mesh' in record and record['mesh'] != cfg:
            raise ValueError('Changing mesh dimensions')
        record['mesh'] = cfg
        frame = read_frame(run / 'foam', stamp)
        surface = np.atleast_1d(np.genfromtxt(output / f'dual_geometry/{stamp}.csv', delimiter=',', names=True))
        n = cfg['nx'] * cfg['ny']
        if len(frame) != n or len(surface) != n:
            raise ValueError('Cell count mismatch')
        surface = surface[np.argsort(surface['cell'])]
        if not np.array_equal(surface['cell'], np.arange(n)):
            raise ValueError('Invalid surface cell indices')
        expected_xy = np.array([[(i % cfg['nx'] + .5) * cfg['H'],
                                 (i // cfg['nx'] + .5) * cfg['H']] for i in range(n)])
        if np.max(np.abs(frame[:, :2] * 1e6 - expected_xy)) > 2e-5:
            raise ValueError('Frame ordering or geometry extent mismatch')
        phase = frame[:, [COLUMNS.index(m) for m in MINERALS]]
        if time == 0 and (polygons[0] or polygons[1] or np.max(abs(phase[:, :2])) > 1e-30):
            raise ValueError('Initial geometry contains pre-existing halite or sylvite')
        mapping = float(np.max(np.abs(fractions - phase)))
        max_mapping = max(max_mapping, mapping)
        if mapping > 6e-8:
            raise ValueError(f'Geometry and phase field mismatch: {name}/{stamp}: {mapping}')
        mobile = frame[:, COLUMNS.index('epsilon')]
        retained = surface['retained_water_fraction']
        water = mobile + retained
        if np.max(np.abs(water + phase.sum(axis=1) - 1)) > 2e-7:
            raise ValueError('Solid/water volume closure failed')
        fields = {c: compact(frame[:, COLUMNS.index('m' + c)]) for c in ('K', 'Mg', 'Na', 'Cl')}
        fields.update(speed=compact(np.hypot(frame[:, COLUMNS.index('Ux_m_s')], frame[:, COLUMNS.index('Uy_m_s')]) * 1e6),
                      water=compact(water), mobile=compact(mobile), retained=compact(retained))
        for mineral in MINERALS:
            covered = surface[mineral + '_covered_m2']
            total_area = covered + surface[mineral + '_area_m2']
            fields['cover_' + mineral] = compact(np.divide(100 * covered, total_area,
                out=np.full(n, np.nan), where=total_area > 0))
        reactive = table(output / f'dual_geometry/{stamp}.reactors.csv')
        local = np.zeros((n, 2))
        mineral_K = np.zeros((3, 2))
        for row in reactive:
            cell = int(row['cell'])
            a = np.array([float(row['A_' + m]) for m in MINERALS])
            k = float(row['mK'])
            local[cell] += [a.sum() * k, a.sum()]
            mineral_K[:, 0] += a * k
            mineral_K[:, 1] += a
        fields['localK'] = compact(np.divide(local[:, 0], local[:, 1], out=np.full(n, np.nan), where=local[:, 1] > 0))
        reaction = np.column_stack([surface[m + '_reaction_area_m2'] for m in MINERALS])
        covered = np.column_stack([surface[m + '_covered_m2'] for m in MINERALS])
        wet = np.column_stack([surface[m + '_area_m2'] for m in MINERALS])
        fields['area'] = compact(reaction.sum(axis=1) * 1e12)
        index = int(np.argmin(abs(hist['time_s'] - time)))
        if abs(float(hist['time_s'][index]) - time) > 1e-9:
            raise ValueError('Missing matching history sample')
        amounts = [float(hist['min_' + m][index] * 1e12) for m in MINERALS]
        if min(amounts) < -1e-8:
            raise ValueError('Negative mineral inventory')
        for k, mineral in enumerate(MINERALS):
            expected = hist['min_' + mineral][0] + hist['prec_' + mineral][index] - hist['diss_' + mineral][index]
            if abs(amounts[k] * 1e-12 - expected) > 1e-8 * hist['min_Carnallite'][0]:
                raise ValueError('Mineral inventory disagrees with precipitation/dissolution history')
        potential = covered.sum(axis=0) + wet.sum(axis=0)
        record['curve'].append({
            'time': time, 'amount': amounts,
            'dissolved': [float(hist['diss_' + m][index] * 1e12) for m in MINERALS],
            'precipitated': [float(hist['prec_' + m][index] * 1e12) for m in MINERALS],
            'area': (reaction.sum(axis=0) * 1e9).tolist(),
            'contact': compact(np.divide(covered.sum(axis=0), potential, out=np.zeros(3), where=potential > 0) * 100),
            'localK': compact(np.divide(mineral_K[:, 0], mineral_K[:, 1], out=np.full(3, np.nan), where=mineral_K[:, 1] > 0)),
            'nuclei': [int(surface['nuclei_' + m].sum()) for m in MINERALS],
        })
        coverage_focus = {}
        for key, areas in [('all', covered.sum(axis=1)), *[(m, covered[:, k]) for k, m in enumerate(MINERALS)]]:
            cell = int(np.argmax(areas))
            coverage_focus[key] = expected_xy[cell].tolist() if areas[cell] > 0 else None
        record['frames'].append({'polygons': polygons, 'fields': fields, 'coverageFocus': coverage_focus})
    record['metrics']['viewer_max_geometry_mapping_error'] = max_mapping
    if np.any(np.diff(np.array([r['nuclei'] for r in record['curve']]), axis=0) < 0):
        raise ValueError('Nucleation history is not cumulative')
    if name.startswith('coupled_') and until is None:
        for phase, mineral in enumerate(('Halite', 'Sylvite')):
            if record['curve'][-1]['nuclei'][phase] <= 0 or hist['prec_' + mineral][-1] <= 0:
                raise ValueError(f'The selected result does not demonstrate new {mineral} precipitation: {name}')
    return record


def export(study, output, until=None, tag=''):
    runs = [build_run(study, case.name + tag, case.title, case.description, until=until)
            for case in COMPARISON_CASES]
    if any(r['times'] != runs[0]['times'] or r['mesh'] != runs[0]['mesh'] for r in runs):
        raise ValueError('Cases have different stored times or extents')
    ranges = {}
    for key in runs[0]['frames'][0]['fields']:
        values = [v for run in runs for frame in run['frames'] for v in frame['fields'][key] if v is not None]
        lo, hi = min(values, default=0.), max(values, default=0.)
        if key in ('water', 'mobile', 'retained'):
            lo, hi = (0., 1.) if key != 'retained' else (0., hi)
        elif key in ('K', 'Mg', 'speed', 'localK', 'area'):
            lo = 0.
        elif key.startswith('cover_'):
            lo, hi = 0., 100.
        ranges[key] = [lo, hi if hi > lo else lo + 1]
    packed = []
    for run in runs:
        raw = json.dumps(run, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()
        packed.append(base64.b64encode(zlib.compress(raw, 9)).decode())
        print(f'{run["id"]}: {len(run["times"])} frames, {len(raw)/1e6:.2f} MB expanded, '
              f'mapping error {run["metrics"]["viewer_max_geometry_mapping_error"]:.3g}', flush=True)
    metadata = {'cases': [{'id': r['id'], 'title': r['title'], 'description': r['description']} for r in runs],
                'ranges': ranges, 'times': runs[0]['times'], 'mesh': runs[0]['mesh'],
                'display_coordinate_unit_um': DISPLAY_UNIT_UM,
                'complete': all(r['complete'] for r in runs),
                'target_end_time_s': runs[0]['target_end_time_s']}
    template = Path(__file__).with_name('reactive_results_viewer.html').read_text(encoding='utf-8')
    for placeholder in ('__METADATA__', '__PAYLOADS__'):
        if template.count(placeholder) != 1:
            raise ValueError(f'Template placeholder mismatch: {placeholder}')
    html = template.replace('__METADATA__', json.dumps(metadata, ensure_ascii=False, separators=(',', ':')))
    html = html.replace('__PAYLOADS__', json.dumps(packed, separators=(',', ':')))
    output.write_text(html, encoding='utf-8')
    print(f'{output}\n{output.stat().st_size/1e6:.2f} MB; {sum(len(r["times"]) for r in runs)} verified frames')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, default=STUDY)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--tag', default='', help='Suffix used with dual_resolution.py --comparison')
    parser.add_argument('--until', type=float, help='Explicitly labelled progress preview through this saved time')
    args = parser.parse_args()
    if '/' in args.tag or '\\' in args.tag:
        parser.error('tag must not contain a path separator')
    export(args.study, args.output or args.study / f'三矿物反应运移_交互结果{args.tag}.html', args.until, args.tag)
