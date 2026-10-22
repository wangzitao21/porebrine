"""Verify bounded liquid-only routing across a thin four-cell pore strip."""
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.mesh import make_case, head


def execute(args, work, log, success=True):
    with log.open("w") as stream:
        result = subprocess.run(args, cwd=work, stdout=stream, stderr=subprocess.STDOUT, timeout=120)
    if success and result.returncode:
        raise RuntimeError(log.read_text()[-4000:])
    if not success and (not result.returncode or "configured open-link budget" not in log.read_text()):
        raise AssertionError(f"Expected physically inaccessible routing: {log}")


def make(work, hops, blocked=False, liquid_fraction=0, isolated=False):
    h, ep = 5e-6, 1e-5
    mineral = np.full((2, 5), 1-ep)
    mineral[1, :4] = .9 * (1-ep)
    mineral[1, 4] = 0
    if blocked:
        mineral[1, 2] = 1-ep
    field, put = make_case(work, 1-mineral, spacing=h, ports=2)
    field("inert", 0)
    field("normalX", 0)
    field("normalY", 1)
    for j in range(7):
        field(f"mineral{j}", mineral if j == 0 else 0)
    put("constant/reactiveProperties", head("reactiveProperties") + f"""
database "{ROOT}/src/chemistry/pitzer.dat";
dx {h}; dy {h}; epsilonMin {ep};
diffusivity [0 2 -1 0 0 0 0] 1e-9;
mineralInterfaceModel heterogeneous;
solidOwnerSurfacePolicy fluidSide;
surfaceRoutingHops {hops};
liquidControlVolumeFraction {liquid_fraction};
isolatedLiquidControlVolumes {'true' if isolated else 'false'};
""")
    put("system/decomposeParDict", head("decomposeParDict") + "numberOfSubdomains 2;method scotch;\n")


def rows(work, ranks):
    folders = [work] if ranks == 1 else [work / f"processor{i}" for i in range(ranks)]
    result = []
    for folder in folders:
        with (folder / "surface_routes.csv").open() as stream:
            result.extend({key: float(value) for key, value in row.items()} for row in csv.DictReader(stream))
    return sorted(result, key=lambda row: tuple(row[k] for k in ("sx", "sy", "rx", "ry", "area")))


def main():
    output = ROOT / "cases/carnallite_coating/output/routing_checks"
    output.mkdir(parents=True, exist_ok=True)
    binary = ROOT / "bin/tests/verifySurfaceRouting"
    foam_src = Path(os.environ["FOAM_SRC"])
    foam_lib = Path(os.environ["FOAM_LIBBIN"])
    mpi_lib = foam_lib / os.environ["FOAM_MPI"]
    command = ["g++", "-std=c++17", "-O2", "-fopenmp", "-DWM_DP", "-DWM_LABEL_SIZE=32", "-DNoRepository"]
    command += [f"-I{ROOT / folder}" for folder in ("src/chemistry", "src/geometry", "src/flow")]
    command += [f"-I{foam_src / folder / 'lnInclude'}" for folder in ("OpenFOAM", "OSspecific/POSIX", "finiteVolume", "meshTools")]
    command += [str(ROOT / "src/tests/verifySurfaceRouting.C"), f"-L{foam_lib}", f"-L{mpi_lib}",
                f"-Wl,-rpath,{foam_lib}", f"-Wl,-rpath,{mpi_lib}",
                "-lfiniteVolume", "-lmeshTools", "-lOpenFOAM", "-lPstream", "-o", str(binary)]
    execute(command, ROOT, output / "compile.log")
    collected = {}
    with tempfile.TemporaryDirectory(prefix="porebrine-route-check-") as tmp:
        for name, hops, blocked, ranks, success, liquid_fraction, isolated in (
            ("two_links_reject", 2, False, 1, False, 0, False),
            ("four_links_serial", 4, False, 1, True, 0, False),
            ("four_links_parallel", 4, False, 2, True, 0, False),
            ("solid_barrier_reject", 4, True, 1, False, 0, False),
            ("finite_liquid_serial", 2, False, 1, True, .05, False),
            ("finite_liquid_parallel", 2, False, 2, True, .05, False),
            ("isolated_pocket_serial", 2, True, 1, True, .05, False),
            ("pocket_fallback_serial", 4, True, 1, True, .35, True),
            ("pocket_fallback_parallel", 4, True, 2, True, .35, True),
        ):
            work = Path(tmp) / name
            make(work, hops, blocked, liquid_fraction, isolated)
            execute(["blockMesh"], work, output / f"{name}_mesh.log")
            args = [str(binary)]
            if ranks > 1:
                execute(["decomposePar"], work, output / f"{name}_decompose.log")
                args = ["mpirun", "--bind-to", "none", "-np", str(ranks), str(binary), "-parallel"]
            execute(args, work, output / f"{name}.log", success)
            if success:
                values = rows(work, ranks)
                area = sum(r["area"] for r in values)
                assert abs(area / (5e-6)**2 - (5.1 if blocked else 5.9)) < 1e-9, area
                far = [r for r in values if r["sx"] == 0 and r["sy"] == 1]
                if liquid_fraction:
                    assert far and all(r["rx"] == 0 and r["ry"] == 1 and r["velocity"] == -1 for r in far)
                else:
                    assert far and all(r["rx"] == 4 and r["ry"] == 1 and 0 < r["velocity"] < 1e-3 for r in far)
                assert all(r["sy"] == 1 or r["sx"] == 4 for r in values)
                assert all(r["molar_flux"] > 0 for r in values)
                collected[name] = np.array([[r[k] for k in ("sx", "sy", "rx", "ry", "area", "velocity", "molar_flux")] for r in values])
            print(f"PASS {name}", flush=True)
        for prefix in ("four_links", "finite_liquid", "pocket_fallback"):
            np.testing.assert_allclose(collected[prefix+"_serial"], collected[prefix+"_parallel"],
                                       rtol=1e-11, atol=1e-25)
        print("PASS serial_parallel_route_area_and_flux", flush=True)
        link_binary = ROOT / "bin/tests/verifySurfaceLinkExchange"
        link_command = command.copy()
        link_command[link_command.index(str(ROOT / "src/tests/verifySurfaceRouting.C"))] = str(ROOT / "src/tests/verifySurfaceLinkExchange.C")
        link_command[-1] = str(link_binary)
        execute(link_command, ROOT, output / "compile_links.log")
        for distance in (2, 4):
            for ranks in (1, 2):
                name = f"all_links_distance_{distance}_ranks_{ranks}"
                work = Path(tmp) / name
                _, put = make_case(work, np.ones((7, 9)), spacing=5e-6, ports=7)
                put("constant/linkTest", head("linkTest") +
                    f"nx 9;ny 7;spacing 5e-6;maximumDistance {distance};\n")
                put("system/decomposeParDict", head("decomposeParDict") + "numberOfSubdomains 2;method scotch;\n")
                execute(["blockMesh"], work, output / f"{name}_mesh.log")
                args = [str(link_binary)]
                if ranks > 1:
                    execute(["decomposePar"], work, output / f"{name}_decompose.log")
                    args = ["mpirun", "--bind-to", "none", "-np", str(ranks), str(link_binary), "-parallel"]
                execute(args, work, output / f"{name}.log")
                print(f"PASS {name}", flush=True)


if __name__ == "__main__":
    main()
