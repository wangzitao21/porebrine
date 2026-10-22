"""Audit the raw precipitate fields and compare representation closures."""
from pathlib import Path
import argparse
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.fields import COLUMNS, frames, read_frame
from carnallite_coating import interface_rows, save_rows
from plot_carnallite_coating import inspect, parameters


def audit(run):
    settings = parameters(run)
    retained = settings.get("advanced.nucleusResolutionModel", "plic") == "population"
    work = run / "foam"
    records = []
    for stamp in frames(work):
        field = read_frame(work, stamp)
        surface = sorted(interface_rows(work, stamp),
                         key=lambda row: (round(row["y_m"], 12), round(row["x_m"], 12)))
        radius = np.hypot(field[:, 0] * 1e6 - 120, field[:, 1] * 1e6 - 100)
        row = {"run": run.name, "representation": "population" if retained else "plic",
               "time_s": float(stamp)}
        for mineral in ("Halite", "Sylvite"):
            total = field[:, COLUMNS.index(mineral)]
            subgrid = np.array([r[f"subgrid_fraction_{mineral}"] for r in surface])
            births = np.array([r[f"births_{mineral}"] for r in surface])
            active = total > 1e-6
            row[f"{mineral}_radius_max_um_above_1e-6"] = float(radius[active].max()) if active.any() else 0
            row[f"{mineral}_cells_beyond_55um"] = int(np.sum(active & (radius > 55)))
            row[f"{mineral}_stock_without_local_birth_cells"] = int(np.sum(active & (births == 0)))
            row[f"{mineral}_promotions"] = sum(r[f"promoted_{mineral}"] for r in surface)
            row[f"{mineral}_total_minus_population_fraction_max"] = float(np.max(total - subgrid))
            row[f"{mineral}_volume_fraction_max"] = float(total.max())
            if retained:
                # These are invariants of the immobile-population closure,
                # not an assertion that real crystals cannot cross cells.
                np.testing.assert_allclose(total, subgrid, atol=3e-8, rtol=1e-6)
                assert row[f"{mineral}_stock_without_local_birth_cells"] == 0
                assert row[f"{mineral}_promotions"] == 0
        records.append(row)
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy", default="porous_final")
    parser.add_argument("--retained", default="porous_population")
    parser.add_argument("--fine", default="porous_population_fine")
    parser.add_argument("--validation", nargs="*", default=[])
    args = parser.parse_args()
    study = ROOT / "cases/carnallite_coating/output"
    names = [args.legacy, args.retained, args.fine]
    records = [row for name in names for row in audit(study / name)]
    save_rows(study / "morphology_audit.csv", records)
    old, new = (read_frame(study / name / "foam", "1") for name in names[:2])
    early_error = np.max(np.abs(old[:, COLUMNS.index("Halite")] - new[:, COLUMNS.index("Halite")]))
    print(f"Pre-promotion halite fraction difference at 1 s: {early_error:.9g}")
    for name in names:
        print(next(row for row in reversed(records) if row["run"] == name))
    save_rows(study / "numerical_validation.csv",
              [inspect(study / name)[0] for name in [*names, *args.validation]])
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.1), constrained_layout=True)
    fields = [read_frame(study / name / "foam", frames(study / name / "foam")[-1]) for name in names]
    maximum = max(float(f[:, COLUMNS.index("Halite")].max()) for f in fields)
    for ax, name, field, title in zip(axes, names, fields,
            ("Previous PLIC conversion / 5 µm",
             "Retained populations / 5 µm",
             "Retained populations / 2.5 µm")):
        settings = parameters(study / name)
        image = ax.imshow(field[:, COLUMNS.index("Halite")].reshape(
            int(settings["mesh.ny"]), int(settings["mesh.nx"])), origin="lower",
            extent=(0, 320, 0, 200), interpolation="nearest", cmap="Blues", vmin=0, vmax=maximum)
        ax.add_patch(Circle((120, 100), 40, fill=False, ls="--", lw=.9, ec="#555"))
        ax.set(title=title, xlabel="x (µm)", ylabel="y (µm)", xlim=(30, 220), ylim=(5, 195))
    fig.colorbar(image, ax=axes, shrink=.8, label="Halite volume fraction (linear, same scale)")
    fig.suptitle("Raw cell-average fields at 3 s · no smoothing or spatial filtering")
    fig.savefig(study / "morphology_comparison.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
