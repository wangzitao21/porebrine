"""Inspect and compare native experiment output; all panels use stored data."""
from pathlib import Path
import argparse
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.fields import COLUMNS, FINAL_COLUMNS, frames, history, metrics, read_frame, final_state
from carnallite_coating import interface_rows, save_rows

VM = np.array([27.1, 37.5, 173.7, 0, 0, 0, 0]) * 1e-6
LABEL = {"dissolution_only": "Dissolution only",
         "precipitation_unblocked": "Precipitation, no surface barrier",
         "porous": "Precipitation, porous coating",
         "impermeable": "Precipitation, impermeable coating"}
COLOR = {"dissolution_only": "#555555", "precipitation_unblocked": "#0072B2",
         "porous": "#D55E00", "impermeable": "#009E73"}
INDEX = {name: i for i, name in enumerate(COLUMNS)}


def parameters(run):
    with (run / "parameters.csv").open() as stream:
        return {r["parameter"]: r["value"] for r in csv.DictReader(stream)}


def inspect(run):
    work = run / "foam"
    hist, diagnostic, settings = history(work), metrics(work), parameters(run)
    interface = np.atleast_1d(np.genfromtxt(run / "interface_history.csv", delimiter=",", names=True, deletechars=""))
    state = final_state(work)
    minerals = state[:, 15:22]
    solid_fraction = minerals[:, :3] @ VM[:3]
    volume = float(settings["mesh.spacing_m"]) ** 2 * float(settings["mesh.thickness_m"])
    initial_moles = hist["min_Carnallite"][0]
    removed = initial_moles - hist["min_Carnallite"][-1]
    product_volume = sum(hist[f"min_{name}"][-1] * VM[j] for j, name in enumerate(("Halite", "Sylvite")))
    birth_sylvite = interface["time_s"][interface["births_Sylvite"] > 0]
    birth_halite = interface["time_s"][interface["births_Halite"] > 0]
    row = {
        "run": run.name, "mode": settings["title"].split(": ")[-1],
        "nucleus_representation": settings.get("advanced.nucleusResolutionModel", "plic"),
        "nx": int(settings["mesh.nx"]), "ny": int(settings["mesh.ny"]),
        "time_step_s": float(settings["simulation.time_step_s"]),
        "end_time_s": diagnostic["physical_time_s"],
        "initial_Carnallite_mol": initial_moles,
        "final_Carnallite_mol": hist["min_Carnallite"][-1],
        "Carnallite_dissolved_fraction": removed / initial_moles,
        "net_Carnallite_removed_volume_m3": removed * VM[2],
        "final_Halite_mol": hist["min_Halite"][-1],
        "final_Sylvite_mol": hist["min_Sylvite"][-1],
        "net_precipitate_volume_m3": product_volume,
        "precipitate_to_removed_volume_ratio": product_volume / (removed * VM[2]) if removed else 0,
        "initial_porosity": hist["porosity"][0], "final_porosity": hist["porosity"][-1],
        "Qin_initial_m3_s": hist["Qin_m3_s"][0], "Qin_final_m3_s": hist["Qin_m3_s"][-1],
        "coverage_final": interface["coverage_Carnallite_area_weighted"][-1],
        "coverage_peak": max(interface["coverage_Carnallite_area_weighted"]),
        "thickness_max_m": max(interface["layer_thickness_max_m"]),
        "first_Halite_birth_frame_s": birth_halite[0] if len(birth_halite) else "",
        "first_Sylvite_birth_frame_s": birth_sylvite[0] if len(birth_sylvite) else "",
        "mixed_cells_final": int(interface["mixed_mineral_cells_above_1e-6"][-1]),
        "three_mineral_cells_final": int(interface["three_mineral_cells_above_1e-6"][-1]),
        "three_mineral_cells_peak": int(max(interface["three_mineral_cells_above_1e-6"])),
        "births_Halite": interface["births_Halite"][-1],
        "births_Sylvite": interface["births_Sylvite"][-1],
        "promoted_Halite": interface["promoted_Halite"][-1],
        "promoted_Sylvite": interface["promoted_Sylvite"][-1],
        "mass_balance_error_max": diagnostic["max_component_balance_error"],
        "phase_volume_closure_error": np.max(np.abs(state[:, 2] + solid_fraction - 1)),
        "minimum_final_mineral_mol_m3": minerals.min(),
        "solver_wall_s": diagnostic["solver_wall_s"],
        "kinetic_relative_tolerance": float(settings["solver.kinetic_relative_tolerance"]),
    }
    if row["phase_volume_closure_error"] > 1e-9 or minerals.min() < -1e-9:
        raise ValueError(f"Invalid final mineral volume ledger: {run}")
    if abs(np.sum(minerals[:, 2]) * volume - hist["min_Carnallite"][-1]) > initial_moles * 1e-10:
        raise ValueError(f"Final inventory differs from integrated history: {run}")
    return row, hist, interface


def map_axes(ax, values, nx, ny, *, cmap, vmin=0, vmax=None, title=""):
    im = ax.imshow(values.reshape(ny, nx), origin="lower",
                   extent=[0, 320, 0, 200], interpolation="nearest",
                   cmap=cmap, vmin=vmin, vmax=vmax)
    ax.add_patch(Circle((120, 100), 40, fill=False, ls="--", lw=.8, ec="#555555"))
    ax.set(xlabel="x (µm)", ylabel="y (µm)", title=title, xlim=(50, 200), ylim=(30, 170))
    return im


def plot(study, names):
    runs = [study / name for name in names]
    inspected = [inspect(run) for run in runs]
    save_rows(study / "comparison.csv", [item[0] for item in inspected])
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for run, (row, hist, interface) in zip(runs, inspected):
        mode = row["mode"]
        label = LABEL.get(mode, run.name)
        if run.name != mode:
            label += f" ({run.name})"
        color = COLOR[mode]
        time = hist["time_s"]
        axes[0, 0].plot(time, hist["min_Carnallite"] / hist["min_Carnallite"][0],
                        label=label, color=color)
        axes[0, 1].plot(time, (hist["min_Halite"] * VM[0] + hist["min_Sylvite"] * VM[1]) /
                       (hist["min_Carnallite"][0] * VM[2]) * 100, color=color)
        axes[1, 0].plot(interface["time_s"], interface["coverage_Carnallite_area_weighted"] * 100,
                       color=color)
        axes[1, 1].plot(time, (hist["porosity"] - hist["porosity"][0]) * 100, color=color)
    for ax in axes.flat:
        ax.set_xlabel("Time (s)")
        ax.grid(alpha=.17)
    axes[0, 0].set_ylabel("Carnallite remaining / initial")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set_ylabel("Product volume / initial rock volume (%)")
    axes[1, 0].set_ylabel("Area-weighted attached coverage (%)")
    axes[1, 1].set_ylabel("Porosity gain (percentage points)")
    fig.suptitle("Saturated-NaCl injection into a 2D carnallite-grain domain\n"
                 "Exploratory kinetics and nucleation parameters; constant pressure drop", fontsize=12)
    fig.savefig(study / "comparison.png", dpi=170)
    plt.close(fig)

    selected = next((r for r, item in zip(runs, inspected) if item[0]["mode"] == "porous"), runs[-1])
    work = selected / "foam"
    stamp = frames(work)[-1]
    frame = read_frame(work, stamp)
    nx, ny = len(np.unique(frame[:, 0])), len(np.unique(frame[:, 1]))
    rows = sorted(interface_rows(work, stamp), key=lambda row: (round(row["y_m"], 12), round(row["x_m"], 12)))
    maps = [
        (frame[:, INDEX["Carnallite"]], "Greys", 1, "Carnallite volume fraction"),
        (frame[:, INDEX["Halite"]], "Blues", None, "Halite volume fraction"),
        (frame[:, INDEX["Sylvite"]], "YlOrBr", None, "Sylvite volume fraction"),
        (np.array([r["coverage_Carnallite"] if r["bare_m2_Carnallite"] > 0 else np.nan for r in rows]),
         "magma", 1, "Attached coverage of carnallite"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for ax, (values, cmap, vmax, title) in zip(axes.flat, maps):
        im = map_axes(ax, values, nx, ny, cmap=cmap, vmax=vmax, title=title)
        fig.colorbar(im, ax=ax, shrink=.8)
    fig.suptitle(f"Porous-coating case at {float(stamp):g} s · dashed circle: initial rock\n"
                 "Stored cell averages; individual crystal positions are not resolved", fontsize=12)
    fig.savefig(study / "mineral_maps.png", dpi=170)
    plt.close(fig)

    snapshots = frames(work)
    wanted = [snapshots[0], min(snapshots, key=lambda t: abs(float(t) - .5)),
              min(snapshots, key=lambda t: abs(float(t) - 1.5)), snapshots[-1]]
    colors = np.array([[.27, .29, .31], [.10, .47, .76], [.94, .62, .13]])
    fig, axes = plt.subplots(1, 4, figsize=(13, 4), constrained_layout=True)
    for ax, t in zip(axes, wanted):
        f = read_frame(work, t)
        phase = f[:, [INDEX["Carnallite"], INDEX["Halite"], INDEX["Sylvite"]]]
        rgb = np.clip(1 - phase.sum(axis=1)[:, None] + phase @ colors, 0, 1)
        ax.imshow(rgb.reshape(ny, nx, 3), origin="lower", extent=[0, 320, 0, 200], interpolation="nearest")
        ax.add_patch(Circle((120, 100), 40, fill=False, ls="--", lw=.7, ec="#777777"))
        ax.set(xlim=(50, 200), ylim=(30, 170), xlabel="x (µm)", ylabel="y (µm)", title=f"{float(t):g} s")
    fig.legend(handles=[Patch(color=c, label=n) for c, n in zip(colors, ["Carnallite", "Halite", "Sylvite"])],
               loc="outside lower center", ncol=3)
    fig.suptitle("Mineral-volume evolution · porous coating\nCell-average mixtures; white = pore liquid", fontsize=12)
    fig.savefig(study / "evolution.png", dpi=170)
    plt.close(fig)

    for row, _, _ in inspected:
        print(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("runs", nargs="+")
    args = parser.parse_args()
    plot(ROOT / "cases/carnallite_coating/output", args.runs)
