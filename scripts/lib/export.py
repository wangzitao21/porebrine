import csv
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from . import ROOT
from .fields import COLUMNS, FINAL_COLUMNS, frames, read_frame, final_state, partitions, spatial_order


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def export(work, run_directory):
    work, target = Path(work).resolve(), Path(run_directory).resolve()
    if work.parent != target or work.name != '.work':
        raise ValueError('Export only from this run\'s own .work directory')
    data, record = target/'data', target/'record'
    data.mkdir(exist_ok=True)
    if any(data.iterdir()):
        raise ValueError('Result data directory must be empty')
    parts = partitions(work)
    output_parts = [part/'output' for part in parts]
    (data/'frames').mkdir()
    times = frames(work)
    if not times:
        raise ValueError('Solver produced no snapshots')
    if len(parts) == 1:
        # Keep native float32/float64 bytes; renaming avoids a second large copy.
        for stamp in times:
            shutil.move(str(output_parts[0]/'frames'/f'{stamp}.f32'), data/'frames'/f'{stamp}.f32')
        shutil.move(str(output_parts[0]/'final_state.f64'), data/'final_state.f64')
    else:
        for stamp in times:
            values = read_frame(work, stamp).astype('<f4')
            with (data/'frames'/f'{stamp}.f32').open('wb') as stream:
                np.array([0x53414C54, 1, len(values), len(COLUMNS)], dtype='<u4').tofile(stream)
                values.tofile(stream)
        values = final_state(work).astype('<f8')
        with (data/'final_state.f64').open('wb') as stream:
            np.array([0x53414C54, 2, len(values), len(FINAL_COLUMNS)], dtype='<u4').tofile(stream)
            values.tofile(stream)
    for name in ('summary.csv', 'history.csv', 'boundary_history.csv', 'time_steps.csv'):
        shutil.move(str(output_parts[0]/name), data/name)
    interface_parts = [folder/'interfaces' for folder in output_parts]
    if (output_parts[0]/'dual_geometry').is_dir():
        shutil.move(str(output_parts[0]/'dual_geometry'), data/'dual_geometry')
    if any(folder.exists() for folder in interface_parts):
        if not all(folder.exists() for folder in interface_parts):
            raise ValueError('Missing interface diagnostics on an MPI rank')
        if len(parts) == 1:
            shutil.move(str(interface_parts[0]), data/'interfaces')
        else:
            (data/'interfaces').mkdir()
            for stamp in times:
                rows = []
                for folder in interface_parts:
                    with (folder/f'{stamp}.csv').open() as stream:
                        rows.extend(csv.DictReader(stream))
                coordinates = np.array([[float(r['x_m']), float(r['y_m'])] for r in rows])
                rows = [rows[i] for i in spatial_order(coordinates)]
                with (data/'interfaces'/f'{stamp}.csv').open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                    writer.writeheader(); writer.writerows(rows)
    profiles = []
    for folder in output_parts:
        with (folder/'chemistry_cells.csv').open() as stream:
            profiles.extend(csv.DictReader(stream))
    profiles.sort(key=lambda row: int(row['global_cell']))
    if not profiles or [int(row['global_cell']) for row in profiles] != list(range(len(profiles))):
        raise ValueError('Incomplete or duplicate chemistry cell profile')
    with (data/'chemistry_cells.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=profiles[0].keys())
        writer.writeheader(); writer.writerows(profiles)
    for name, columns in (('frame_columns.csv', COLUMNS), ('final_state_columns.csv', FINAL_COLUMNS)):
        with (data/name).open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(['column', 'field'])
            writer.writerows(enumerate(columns))
    paths = sorted((ROOT/'src').rglob('*'))
    paths += sorted(p for p in (ROOT/'scripts').rglob('*') if 'trash' not in p.relative_to(ROOT/'scripts').parts)
    paths += [ROOT/'bin/saltLakeFoam', ROOT/'bin/libSaltChemistry.so']
    with (record/'program_hashes.csv').open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['path', 'bytes', 'sha256'])
        for path in paths:
            if path.is_file() and '__pycache__' not in path.parts:
                writer.writerow([path.relative_to(ROOT), path.stat().st_size, sha256(path)])
    settings = {str(path.relative_to(work)): path.read_text()
                for folder in ('constant', 'system') for path in (work/folder).glob('*') if path.is_file() and path.suffix != '.bin'}
    (record/'generated_settings.json').write_text(json.dumps(settings, indent=2)+'\n', encoding='utf-8')
    (record/'data_format.json').write_text(json.dumps({
        'byte_order': 'little-endian', 'header': 'four uint32 values: magic, version, cells, columns',
        'snapshots': {'dtype': 'float32', 'columns': COLUMNS},
        'final_state': {'dtype': 'float64', 'columns': FINAL_COLUMNS},
    }, indent=2)+'\n', encoding='utf-8')
