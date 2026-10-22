import argparse
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np

from lib.case_input import case_path, default_config, input_hashes, read_inputs, read_json


def create(reference, template=None, config_file=None, geometry_file=None, initial_file=None,
           nx=16, ny=8, spacing=1e-5, dual_geometry_file=None):
    destination = case_path(reference, create=True)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.prepare-', dir=destination.parent) as temporary:
        source = Path(temporary)/'input_data'
        if template:
            if any((config_file, geometry_file, initial_file, dual_geometry_file)):
                raise ValueError('--template cannot be combined with separate input files')
            template = Path(template)
            if not template.exists():
                template = case_path(template)
            template = template if template.name == 'input_data' else template/'input_data'
            read_inputs(template)
            before = input_hashes(template)
            shutil.copytree(template, source)
            if input_hashes(source) != before or input_hashes(template) != before:
                raise ValueError('Template inputs changed while they were being copied')
        else:
            source.mkdir()
            config = read_json(config_file) if config_file else default_config(nx, ny, spacing)
            (source/'config.json').write_text(json.dumps(config, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            if geometry_file:
                shutil.copy2(geometry_file, source/'geometry.npz')
            else:
                m = config['mesh']
                np.savez_compressed(source/'geometry.npz', solid_fractions=np.zeros((8, m['ny'], m['nx'])))
            if initial_file:
                shutil.copy2(initial_file, source/'initial_conditions.npz')
            if dual_geometry_file:
                shutil.copy2(dual_geometry_file, source/'dual_geometry.bin')
        config, _, _ = read_inputs(source)
        destination.mkdir()  # No replacement of existing cases or their run history.
        shutil.move(str(source), str(destination/'input_data'))
        (destination/'output').mkdir()
    return destination, config


def main():
    parser = argparse.ArgumentParser(description='Create a case with input_data/ and output/; never overwrite an existing case.')
    parser.add_argument('case', help='New name or a path under cases/')
    parser.add_argument('--template', type=Path, help='Copy only the scientific inputs from any existing case')
    parser.add_argument('--config', type=Path, help='Scientific config.json')
    parser.add_argument('--geometry', type=Path, help='geometry.npz: solid_fractions[Halite,Sylvite,Carnallite,Polyhalite,Gypsum,Bischofite,Calcite,inert,y,x]')
    parser.add_argument('--initial-conditions', type=Path, help='Optional component molalities [y,x], pressure_Pa [y,x], velocity_m_s [y,x,2]')
    parser.add_argument('--dual-geometry', type=Path, help='Fine mineral polygon geometry for dualMineralGeometry')
    parser.add_argument('--nx', type=int, default=16)
    parser.add_argument('--ny', type=int, default=8)
    parser.add_argument('--spacing-m', type=float, default=1e-5)
    args = parser.parse_args()
    try:
        path, config = create(args.case, args.template, args.config, args.geometry,
                              args.initial_conditions, args.nx, args.ny, args.spacing_m, args.dual_geometry)
    except (ValueError, OSError, KeyError) as error:
        parser.error(str(error))
    size = sum(p.stat().st_size for p in (path/'input_data').iterdir())
    print(f'Created {path}: {config["mesh"]["nx"]} x {config["mesh"]["ny"]} cells, {size} input bytes')


if __name__ == '__main__':
    main()
