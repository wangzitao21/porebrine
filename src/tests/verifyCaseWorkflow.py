import copy
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import time
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from prepare import create
from run import case_lock, clean_work, run_case
from lib.case_input import default_config, input_hashes, read_inputs, read_json, result_path, validate_config
from lib.foam_input import generate
from lib.native_run import command


class CaseWorkflow(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='.workflow-test-', dir=ROOT/'cases')
        self.addCleanup(self.remove_temporary)
        self.root = Path(self.temporary.name)

    def remove_temporary(self):
        path = Path(self.temporary.name).resolve()
        if path.parent != (ROOT/'cases').resolve() or not path.name.startswith('.workflow-test-'):
            raise ValueError('Unsafe test cleanup')
        self.temporary.cleanup()

    def native_input(self, config=None, solid=None, fields=None):
        directory = self.root/'source/input_data'
        directory.mkdir(parents=True)
        config = config or default_config()
        (directory/'config.json').write_text(json.dumps(config)+'\n')
        if solid is None:
            solid = np.zeros((8, config['mesh']['ny'], config['mesh']['nx']))
        np.savez_compressed(directory/'geometry.npz', solid_fractions=solid)
        if fields is not None:
            np.savez_compressed(directory/'initial_conditions.npz', **fields)
        return directory

    def test_unknown_keys_units_and_overlapping_ports_rejected(self):
        for change in ('unknown', 'unit', 'overlap', 'alignment', 'water', 'rates'):
            config = default_config()
            if change == 'unknown':
                config['physics']['viscosity'] = config['physics'].pop('viscosity_Pa_s')
            elif change == 'unit':
                config['boundaries'][0]['value_m'] = config['boundaries'][0].pop('value_Pa')
            elif change == 'overlap':
                config['boundaries'].append(dict(config['boundaries'][0], name='overlap'))
            elif change == 'alignment':
                config['boundaries'][0]['range_m'][0] = 2e-6
            elif change == 'water':
                config['initial_solution_mol_kg_water']['Na'] = float('nan')
            else:
                del config['reactions']['rate_coefficients_mol_m2_s']['Halite']
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate_config(config)

    def test_geometry_and_spatial_fields_not_silently_reshaped(self):
        directory = self.native_input()
        for values in (np.zeros((8, 16, 8)), np.ones((8, 8, 16)), np.full((8, 8, 16), np.nan), np.zeros((8, 8, 16), dtype=complex)):
            np.savez_compressed(directory/'geometry.npz', solid_fractions=values)
            with self.assertRaises(ValueError):
                read_inputs(directory)
        np.savez_compressed(directory/'geometry.npz', solid_fractions=np.zeros((8, 8, 16)))
        np.savez_compressed(directory/'initial_conditions.npz', Na=np.zeros((8, 16)))
        with self.assertRaises(ValueError):
            read_inputs(directory)

    def test_json_duplicates_rejected_and_windows_bom_accepted(self):
        path = self.root/'values.json'
        path.write_text('{"Na": 1, "Na": 2}')
        with self.assertRaises(ValueError):
            read_json(path)
        path.write_text('{"Na": 1}', encoding='utf-8-sig')
        self.assertEqual(read_json(path), {'Na': 1})

    def test_template_copy_does_not_copy_results_or_replace_case(self):
        source = self.native_input()
        (source.parent/'output/run_0001').mkdir(parents=True)
        (source.parent/'output/run_0001/large_result.bin').write_bytes(b'old data')
        destination, _ = create(self.root/'clone', template=source.parent)
        self.assertEqual({p.name for p in destination.iterdir()}, {'input_data', 'output'})
        self.assertEqual(input_hashes(source), input_hashes(destination/'input_data'))
        self.assertEqual(list((destination/'output').iterdir()), [])
        with self.assertRaises(ValueError):
            create(destination, template=source.parent)

    def test_generated_inputs_preserve_values_and_support_spatial_states(self):
        config = default_config()
        config['mesh']['thickness_m'] = 2e-5
        config['boundaries'].append({'name': 'top_port', 'side': 'north', 'range_m': [0., 2e-5],
                                    'type': 'head', 'value_m': .001,
                                    'solution_mol_kg_water': config['initial_solution_mol_kg_water'].copy()})
        solid = np.zeros((8, 8, 16)); solid[0, 3:5, 7:9] = .25
        initial = {name: np.full((8, 16), value) for name, value in config['initial_solution_mol_kg_water'].items()}
        initial['velocity_m_s'] = np.zeros((8, 16, 2))
        initial['pressure_Pa'] = np.zeros((8, 16))
        source = self.native_input(config, solid, initial)
        before = input_hashes(source)
        work = self.root/'generated'
        generate(source, work, ranks=2)
        self.assertEqual(input_hashes(source), before)
        self.assertEqual({p.name for p in work.iterdir()}, {'0', 'constant', 'system'})
        values = (work/'0/epsilon').read_text().split('internalField nonuniform List<scalar>\n128\n(\n')[1].split('\n)')[0]
        np.testing.assert_array_equal(np.fromstring(values, sep=' '), (1-solid.sum(axis=0)).ravel())
        self.assertIn('nonuniform List<vector>\n128', (work/'0/U').read_text())
        self.assertIn('top_port', (work/'system/createPatchDict').read_text())
        self.assertIn('numberOfSubdomains 2;', (work/'system/decomposeParDict').read_text())

    def test_failed_runs_are_retained_and_never_replace_history(self):
        source = self.native_input()
        case, _ = create(self.root/'case', template=source.parent)
        before = input_hashes(case/'input_data')
        with patch('run.execute', side_effect=RuntimeError('injected solver failure')):
            for _ in range(2):
                with self.assertRaisesRegex(RuntimeError, 'injected solver failure'):
                    run_case(case)
        self.assertEqual(input_hashes(case/'input_data'), before)
        self.assertEqual({p.name for p in case.iterdir()}, {'input_data', 'output'})
        for name in ('run_0001', 'run_0002'):
            output = case/'output'/name
            self.assertTrue((output/'.work/0').is_dir())
            self.assertEqual(json.loads((output/'logs/status.json').read_text())['state'], 'failed')
        with self.assertRaises(ValueError):
            result_path(case)
        with self.assertRaises(ValueError):
            clean_work(case/'input_data', case/'output/run_0001')

    def test_concurrent_run_is_rejected_before_creating_output(self):
        source = self.native_input()
        case, _ = create(self.root/'case', template=source.parent)
        with case_lock(case/'output'):
            with self.assertRaisesRegex(RuntimeError, 'already'):
                run_case(case)
        self.assertEqual(list((case/'output').glob('run_*')), [])

    def test_timeout_terminates_child_process_group(self):
        late_file = self.root/'late.bin'
        child_code = f'from pathlib import Path; import time; time.sleep(.6); Path({str(late_file)!r}).write_bytes(b"late")'
        parent_code = f'import subprocess,sys,time; child=subprocess.Popen([sys.executable,"-c",{child_code!r}]); print(child.pid,flush=True); time.sleep(20)'
        with self.assertRaises(subprocess.TimeoutExpired):
            command([sys.executable, '-c', parent_code], self.root, self.root/'timeout.log', .2)
        time.sleep(.7)
        self.assertFalse(late_file.exists())


if __name__ == '__main__':
    unittest.main()
