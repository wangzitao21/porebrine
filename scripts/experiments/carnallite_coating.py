"""Reproducible exploratory carnallite / saturated-NaCl flow experiment."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.case_input import COMPONENTS, MINERALS, default_config
from lib.chemistry import Chemistry
from lib.foam_input import generate
from lib.native_run import execute
from lib.fields import COLUMNS, frames, read_frame, history, metrics, final_state, partitions
from lib.study import flatten, saturation, save_rows


def probe(destination):
    destination.mkdir(parents=True, exist_ok=True)
    rows = []
    with Chemistry() as chemistry:
        molality = saturation(chemistry)
        for extent in [0, .001, .01, .03, .1, .3, .5, 1, 1.5, 2, 3, 4]:
            water = 1 + 6 * .018016 * extent
            totals = np.array([molality, extent, 0, extent, molality + 3 * extent, 0, 0]) / water
            _, state = chemistry.fractions(totals)
            row = {"carnallite_added_mol_per_initial_kg_water": extent,
                   "free_water_kg": water, "density_kg_m3": state.rho,
                   "water_activity": state.aw}
            row.update({f"molality_{name}": value for name, value in zip(COMPONENTS, totals)})
            row.update({f"SI_{name}": state.si[j] for j, name in enumerate(MINERALS)})
            rows.append(row)
    save_rows(destination / "chemical_path.csv", rows)
    print(f"NaCl saturation molality = {molality:.14g}", flush=True)
    for row in rows:
        print("extent={:.3g} SI(H,S,C)=({:.5g},{:.5g},{:.5g})".format(
            row["carnallite_added_mol_per_initial_kg_water"],
            row["SI_Halite"], row["SI_Sylvite"], row["SI_Carnallite"]), flush=True)
    return molality


def geometry(nx, ny, spacing):
    # A fixed physical grain is integrated by subpixel quadrature at every mesh.
    fractions = np.zeros((8, ny, nx))
    offsets = (np.arange(24) + .5) / 24
    occupied = np.zeros((ny, nx))
    for oy in offsets:
        y = (np.arange(ny)[:, None] + oy) * spacing
        for ox in offsets:
            x = (np.arange(nx)[None, :] + ox) * spacing
            occupied += (x - 120e-6) ** 2 + (y - 100e-6) ** 2 <= (40e-6) ** 2
    fractions[2] = np.minimum(occupied / len(offsets) ** 2, 1 - 1e-5)
    return fractions


def integrated_geometry(nx, ny, spacing):
    """Integrate one physical circle; close only quadrature error on cut cells."""
    radius, cx, cy = 40e-6, 120e-6, 100e-6
    nodes, weights = np.polynomial.legendre.leggauss(48)
    occupied = np.zeros((ny, nx))
    for iy in range(ny):
        yl, yh = iy*spacing-cy, (iy+1)*spacing-cy
        for ix in range(nx):
            xl, xh = ix*spacing-cx, (ix+1)*spacing-cx
            if max(abs(xl), abs(xh))**2+max(abs(yl), abs(yh))**2 <= radius**2:
                occupied[iy, ix] = 1
                continue
            a, b = max(xl, -radius), min(xh, radius)
            if a >= b:
                continue
            cuts = [a, b]
            for y in (yl, yh):
                if abs(y) < radius:
                    q = np.sqrt(radius**2-y*y)
                    cuts += [x for x in (-q, q) if a < x < b]
            cuts = sorted(set(cuts))
            for a, b in zip(cuts[:-1], cuts[1:]):
                x = (a+b)/2+(b-a)/2*nodes
                h = np.sqrt(np.maximum(0, radius**2-x*x))
                height = np.maximum(0, np.minimum(yh, h)-np.maximum(yl, -h))
                occupied[iy, ix] += (b-a)/2*np.dot(weights, height)/spacing**2
    cut = (occupied > 0) & (occupied < 1)
    occupied[cut] *= (np.pi*radius**2/spacing**2-np.sum(occupied == 1))/occupied[cut].sum()
    if occupied.min() < 0 or occupied.max() > 1+1e-12:
        raise RuntimeError('Circle quadrature failed')
    fractions = np.zeros((8, ny, nx))
    fractions[2] = occupied*(1-1e-5)
    return fractions


def settings(mode, molality, refinement, dt, end, ranks, threads, nucleation, seed,
             carnallite_rate, liquid_volume_threshold, routing_hops):
    nx, ny, spacing = 64 * refinement, 40 * refinement, 5e-6 / refinement
    config = default_config(nx, ny, spacing)
    config["title"] = f"Carnallite in saturated NaCl: {mode}"
    config["mesh"]["thickness_m"] = 20e-6
    solution = dict.fromkeys(COMPONENTS, 0.0)
    solution.update(Na=molality, Cl=molality)
    config["initial_solution_mol_kg_water"] = solution
    for port in config["boundaries"]:
        port["solution_mol_kg_water"] = solution.copy()
    config["boundaries"][0]["value_Pa"] = .02
    reaction = config["reactions"]
    reaction["interface_model"] = "heterogeneous"
    reaction["rate_coefficients_mol_m2_s"]["Carnallite"] = carnallite_rate
    reaction["nucleation_seed"] = seed
    reaction["surface_blocking"] = mode not in ("dissolution_only", "precipitation_unblocked")
    reaction["product_layer"].update(model="porous", porosity=.35, tortuosity=2.0)
    if mode == "impermeable":
        reaction["product_layer"]["model"] = "impermeable"
    if mode != "dissolution_only":
        reaction["rate_coefficients_mol_m2_s"].update(Halite=.02, Sylvite=.01)
        for product in ("Halite", "Sylvite"):
            reaction["nucleation"].append({
                "product": product, "substrate": "Carnallite",
                "prefactor_m2_s": nucleation, "contact_angle_deg": 60.0,
                "embryo_radius_m": 1e-7, "effective_barrier": .02,
            })
    config["simulation"].update(end_time_s=end, time_step_s=dt,
                                minimum_time_step_s=1e-7,
                                maximum_time_step_s=dt,
                                adaptive_time_step=False)
    config["solver"].update(chemistry_threads=threads, mpi_ranks=ranks,
                            timeout_s=3600.0, coupling_tolerance=1e-8)
    config["output"].update(frame_interval_s=min(.25, end), make_figures=False)
    config["advanced"] = {"solidOwnerSurfacePolicy": "fluidSide", "surfaceRoutingHops": routing_hops,
                          "liquidControlVolumeFraction": liquid_volume_threshold}
    return config, geometry(nx, ny, spacing)


def run(args):
    study = ROOT / "cases" / "carnallite_coating" / "output"
    study.mkdir(parents=True, exist_ok=True)
    if args.coating_morphology is None:
        molality = probe(study)
    else:
        with Chemistry() as chemistry:
            molality = saturation(chemistry)
    target = study / args.name
    target.mkdir()
    logs = target / "logs"
    logs.mkdir()
    config, solids = settings(args.mode, molality, args.refinement, args.dt,
                             args.end, args.ranks, args.threads, args.nucleation, args.seed,
                             args.carnallite_rate, args.liquid_volume_threshold, args.routing_hops)
    config["output"]["frame_interval_s"] = min(args.frame_interval, args.end)
    config["advanced"]["nucleusResolutionModel"] = args.nucleus_representation
    if args.coating_morphology is not None:
        reaction = config['reactions']
        reaction['nucleation'] = []
        reaction['surface_blocking'] = True
        reaction['product_layer'].update(model='porous', porosity=.35 if args.coating_morphology == 'porous' else 0., tortuosity=2.)
        if args.coating_morphology != 'baseline':
            reaction['rate_coefficients_mol_m2_s'].update(Halite=.02, Sylvite=.01)
            shares = dict.fromkeys(MINERALS, 0.)
            shares.update(Halite=.95, Sylvite=.05)
            reaction['initial_coating'] = {
                'substrate': 'Carnallite', 'morphology': 'particles' if args.coating_morphology == 'particles' else 'film',
                'volume_m3': 2*np.pi*40e-6*20e-6*args.coating_solid_thickness,
                'volume_shares': shares, 'particle_radius_m': 1.5e-6, 'contact_angle_deg': 90.,
            }
        else:
            reaction['rate_coefficients_mol_m2_s'].update(Halite=0., Sylvite=0.)
        # Tiny cut-cell pockets communicate through the existing finite route
        # resistance; a coating is stored only in a sufficiently large liquid
        # receiver, not in epsilonMin numerical water inside the primary.
        config['advanced']['liquidControlVolumeFraction'] = .35
        config['advanced']['isolatedLiquidControlVolumes'] = True
        config['advanced']['surfaceRoutingHops'] = 4
        config['advanced']['interfaceNormalModel'] = 'heightFit'
        # Fine-grid initial Stokes solves can need more SIMPLE iterations;
        # retain the same residual tolerances and increase only the ceiling.
        config['advanced']['maxFlowIterations'] = 8000
        solids = integrated_geometry(config['mesh']['nx'], config['mesh']['ny'], config['mesh']['spacing_m'])
    save_rows(target / "parameters.csv", flatten(config))
    np.savez_compressed(target / "initial_geometry.npz", solid_fractions=solids)
    save_rows(target / "program_hashes.csv", (
        {"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for path in [ROOT / "bin/saltLakeFoam", Path(__file__),
                     *sorted((ROOT / "src").rglob("*.hpp")), *sorted((ROOT / "src/flow").glob("*.H")),
                     ROOT / "src/flow/saltLakeFoam.C", ROOT / "src/chemistry/pitzer.dat"]))
    # The standard generator requires JSON input; keep its adapter files temporary.
    with tempfile.TemporaryDirectory(prefix="porebrine-carnallite-") as temporary:
        inputs = Path(temporary) / "input_data"
        inputs.mkdir()
        (inputs / "config.json").write_text(json.dumps(config), encoding="utf-8")
        np.savez_compressed(inputs / "geometry.npz", solid_fractions=solids)
        generate(inputs, target / "foam", args.ranks)
    print(f"RUN {target} mode={args.coating_morphology or args.mode} cells={solids.shape[1:]} dt={args.dt} end={args.end}", flush=True)
    started = time.monotonic()
    # Concurrent independent runs must not all bind to the same first cores.
    result = execute(target / "foam", logs, args.ranks, args.threads, 3600, bind_cores=False)
    result["total_wall_s"] = time.monotonic() - started
    save_rows(target / "runtime.csv", ({"metric": key, "value": value} for key, value in result.items()))
    analyse(target)


def interface_rows(work, stamp):
    rows = []
    for part in partitions(work):
        with (part / "output/interfaces" / f"{stamp}.csv").open() as stream:
            rows.extend({key: float(value) for key, value in row.items()} for row in csv.DictReader(stream))
    return rows


def analyse(target):
    target = Path(target)
    work = target / "foam"
    fields = {name: i for i, name in enumerate(COLUMNS)}
    rows = []
    for stamp in frames(work):
        frame = read_frame(work, stamp)
        interfaces = interface_rows(work, stamp)
        area = sum(r["bare_m2_Carnallite"] for r in interfaces)
        def weighted(name):
            return sum(r[name] * r["bare_m2_Carnallite"] for r in interfaces) / area if area else 0
        mineral = frame[:, [fields[name] for name in MINERALS]]
        row = {"time_s": float(stamp), "carnallite_support_m2": area,
               "coverage_Carnallite_area_weighted": weighted("coverage_Carnallite"),
               "layer_coverage_Carnallite_area_weighted": weighted("layer_coverage_Carnallite"),
               "layer_transmission_Carnallite_area_weighted": weighted("layer_transmission_Carnallite"),
               "layer_thickness_max_m": max(r["layer_thickness_m_Carnallite"] for r in interfaces),
               "mixed_mineral_cells_above_1e-6": int(np.sum(np.sum(mineral > 1e-6, axis=1) > 1)),
               "three_mineral_cells_above_1e-6": int(np.sum(np.sum(mineral > 1e-6, axis=1) == 3)),
               "minimum_porosity": frame[:, fields["epsilon"]].min(),
               "maximum_solid_sum": mineral.sum(axis=1).max()}
        for name in MINERALS[:3]:
            row[f"solid_fraction_sum_{name}"] = frame[:, fields[name]].sum()
            for prefix in ("births", "nuclei", "promoted", "subgrid_fraction", "attached_fraction", "residue_fraction"):
                row[f"{prefix}_{name}"] = sum(r[f"{prefix}_{name}"] for r in interfaces)
            wet = frame[:, fields["epsilon"]] > .01
            row[f"maximum_wet_SI_{name}"] = frame[wet, fields[f"SI{name}"]].max()
        rows.append(row)
    save_rows(target / "interface_history.csv", rows)
    state = final_state(work)
    if not np.isfinite(state).all():
        raise RuntimeError("Nonfinite final state")
    summary = metrics(work)
    hist = history(work)
    if abs(summary["physical_time_s"] - hist["time_s"][-1]) > 1e-8:
        raise RuntimeError("Missing final history")
    if summary["max_component_balance_error"] > 1e-7:
        raise RuntimeError("Conservation criterion failed")
    print("RESULT", target.name, "wall_s", summary["solver_wall_s"],
          "balance", summary["max_component_balance_error"], flush=True)
    print("FINAL", {key: value for key, value in rows[-1].items()
                    if key in ("time_s", "coverage_Carnallite_area_weighted",
                               "mixed_mineral_cells_above_1e-6", "three_mineral_cells_above_1e-6")
                    or key.startswith(("births_", "solid_fraction_sum_"))}, flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--analyse", type=Path)
    parser.add_argument("--name")
    parser.add_argument("--mode", choices=("dissolution_only", "precipitation_unblocked", "porous", "impermeable"), default="porous")
    parser.add_argument("--end", type=float, default=10.0)
    parser.add_argument("--dt", type=float, default=.02)
    parser.add_argument("--refinement", type=int, default=1)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--nucleation", type=float, default=1e11)
    parser.add_argument("--carnallite-rate", type=float, default=.05)
    # Below the PLIC volume-collapse tolerance: every geometrically retained
    # liquid pocket has its own finite inventory, including isolated pockets.
    parser.add_argument("--liquid-volume-threshold", type=float, default=1e-14)
    parser.add_argument("--routing-hops", type=int, default=2)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--nucleus-representation", choices=("population", "plic"), default="population",
                        help="Retain local precipitate populations, or use the legacy radius-to-PLIC conversion")
    parser.add_argument("--frame-interval", type=float, default=.05)
    parser.add_argument('--coating-morphology', choices=('baseline', 'particles', 'dense', 'porous'))
    parser.add_argument('--coating-solid-thickness', type=float, default=2e-7,
                        help='Fixed total coating volume divided by the physical initial grain area, in m')
    args = parser.parse_args()
    if args.probe:
        probe(ROOT / "cases/carnallite_coating/output")
    elif args.analyse:
        analyse(args.analyse)
    else:
        if not args.name or Path(args.name).name != args.name:
            parser.error("--name must be a new single-directory run name")
        run(args)


if __name__ == "__main__":
    main()
