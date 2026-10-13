import argparse
import json
import subprocess
import tempfile
from pathlib import Path

from lib import ROOT
from lib.case_input import case_path, input_hashes, read_inputs, result_path
from lib.native_run import command
from lib.validation import validate_result


def core_tests():
    tests = {
        'verifyContracts': True, 'verifyDerivatives': True, 'verifyMineralInterface': True,
        'verifyProductTransport': False, 'verifySurfaceReaction': True, 'verifyNucleation': True,
        'verifyCapReactions': True, 'verifySurfaceFilm': True, 'verifySurfacePartition': False,
    }
    for name in tests:
        if not (ROOT/'bin/tests'/name).is_file():
            raise ValueError('Verification binaries are missing; build with --checks first')
    with tempfile.TemporaryDirectory(prefix='.core-check-', dir=ROOT/'bin') as temporary:
        work = Path(temporary)
        for name, database in tests.items():
            args = [str(ROOT/'bin/tests'/name)]
            if database:
                args.append(str(ROOT/'src/chemistry/pitzer.dat'))
            command(args, ROOT, work/f'{name}.log', 60)
            print(f'PASS {name}')
    return len(tests)


def main():
    parser = argparse.ArgumentParser(description='Validate any case inputs and its completed numerical result.')
    parser.add_argument('case', nargs='?', help='Case name or path')
    parser.add_argument('--run', help='Specific run_XXXX; defaults to latest completed run')
    parser.add_argument('--inputs-only', action='store_true')
    parser.add_argument('--core-tests', action='store_true', help='Also run the existing chemistry/interface regression binaries')
    parser.add_argument('--time-refinement', action='store_true', help='Compare the selected result with one dt/2 run; no grid convergence claim')
    args = parser.parse_args()
    if not args.case and not args.core_tests:
        parser.error('Supply a case or --core-tests')
    if args.time_refinement and (not args.case or args.inputs_only):
        parser.error('--time-refinement requires a completed case result')
    try:
        if args.case:
            case = case_path(args.case)
            config, _, _ = read_inputs(case/'input_data')
            print(f'Valid inputs: {case}; {config["mesh"]["nx"]} x {config["mesh"]["ny"]} cells')
            if not args.inputs_only:
                result = result_path(case, args.run)
                report = validate_result(result)
                report['current_inputs_match_run'] = input_hashes(case/'input_data') == input_hashes(result/'record/input_data')
                print(json.dumps(report, indent=2))
                if args.time_refinement:
                    from lib.accuracy import time_refinement
                    accuracy = time_refinement(result)
                    print(f'Passed dt/2 comparison; maximum error / limit = {accuracy["maximum_limit_ratio"]:.6g}')
        if args.core_tests:
            print(f'Passed {core_tests()} core regression checks')
    except (ValueError, OSError, RuntimeError, subprocess.TimeoutExpired) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
