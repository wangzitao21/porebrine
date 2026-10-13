import argparse

from lib.case_input import result_path
from lib.plotting import plot_result


def main():
    parser = argparse.ArgumentParser(description='Make a numerical overview for any completed case run.')
    parser.add_argument('case', help='Case name, case path, or run directory')
    parser.add_argument('--run', help='Specific run_XXXX; defaults to latest completed run')
    args = parser.parse_args()
    print(plot_result(result_path(args.case, args.run)))


if __name__ == '__main__':
    main()
