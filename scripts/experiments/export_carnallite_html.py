"""Export native simulation frames into one offline interactive comparison."""
from pathlib import Path
import argparse
import base64
import csv
import gzip
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.fields import COLUMNS, frames, history, read_frame
from carnallite_coating import interface_rows, save_rows
from plot_carnallite_coating import inspect, parameters

CHANNELS = [
    "Carnallite", "Halite", "Sylvite", "epsilon",
    "coverage_Carnallite", "subgrid_fraction_Halite", "subgrid_fraction_Sylvite",
    "SIHalite", "SISylvite", "SICarnallite", "Ux_m_s", "Uy_m_s",
    "bare_m2_Carnallite",
    "births_Halite", "births_Sylvite",
]


def dataset(run, label):
    settings = parameters(run)
    summary, hist, _ = inspect(run)
    work = run / "foam"
    stamps = frames(work)
    nx, ny = int(settings["mesh.nx"]), int(settings["mesh.ny"])
    data = []
    for stamp in stamps:
        f = read_frame(work, stamp)
        surface = sorted(interface_rows(work, stamp),
                         key=lambda r: (round(r["y_m"], 12), round(r["x_m"], 12)))
        coordinates = np.array([[r["x_m"], r["y_m"]] for r in surface])
        np.testing.assert_allclose(coordinates, f[:, :2], rtol=1e-6, atol=1e-12)
        assert len(f) == nx * ny
        packed = np.column_stack([
            f[:, COLUMNS.index(name)] if name in COLUMNS else
            np.array([r[name] for r in surface]) for name in CHANNELS
        ])
        assert np.isfinite(packed).all()
        np.testing.assert_allclose(packed[:, :4].sum(axis=1), 1, atol=2e-7, rtol=0)
        assert packed[:, :4].min() >= -1e-10
        assert np.all(packed[:, 5] <= packed[:, 1] + 2e-7)
        assert np.all(packed[:, 6] <= packed[:, 2] + 2e-7)
        data.append(packed)
    data = np.array(data, dtype="<f4")
    # Birth counters and native phase fields, rather than an artist's shell,
    # establish whether products can form in the newly opened rock interior.
    initial_solid = data[0, :, 0] > .9999
    initial_inner_volume = data[0, initial_solid, 0].sum()
    volume_um3 = float(settings["mesh.spacing_m"])**2 * float(settings["mesh.thickness_m"]) * 1e18
    distance = np.hypot(f[:, 0] - 120e-6, f[:, 1] - 100e-6) * 1e6
    front_rows = []
    previous_births = np.zeros((nx*ny, 2))
    for stamp, field in zip(stamps, data):
        new_births = field[:, 13:15] - previous_births
        if new_births.min() < -1e-6:
            raise ValueError("Cumulative birth counters went backwards")
        # Initially solid cells need actual dissolution before a birth;
        # epsilonMin numerical pore water alone is not physical access.
        early = initial_solid & (new_births.sum(axis=1) > 0) & (field[:, 0] >= data[0, :, 0] - 1e-7)
        if np.any(early):
            raise ValueError("New nucleus appeared in an unwetted initial solid cell")
        row = {"time_s": float(stamp), "initially_solid_cells": int(initial_solid.sum()),
               "opened_initially_solid_cells": int(np.sum(initial_solid & (field[:, 0] < .9999)))}
        for j, name in enumerate(("Halite", "Sylvite")):
            material = field[:, j+1]
            row[f"volume_um3_{name}"] = float(material.sum())*volume_um3
            row[f"inner_volume_um3_{name}"] = float(material[initial_solid].sum())*volume_um3
            row[f"births_{name}"] = float(field[:, 13+j].sum())
            row[f"inner_births_{name}"] = float(field[initial_solid, 13+j].sum())
            row[f"new_inner_births_{name}"] = float(new_births[initial_solid, j].sum())
            row[f"maximum_inward_cell_centre_um_{name}"] = float(np.maximum(0, 40-distance[material > 1e-6]).max(initial=0))
        front_rows.append(row)
        previous_births = field[:, 13:15].copy()
    save_rows(run / "front_mechanism.csv", front_rows)
    initial = hist["min_Carnallite"][0]
    curves = np.column_stack([
        hist["time_s"], hist["min_Carnallite"] / initial * 100,
        hist["min_Halite"] * 27.1 / (initial * 173.7) * 100,
        hist["min_Sylvite"] * 37.5 / (initial * 173.7) * 100,
        hist["porosity"] * 100,
    ])
    return {
        "id": run.name, "label": label, "nx": nx, "ny": ny,
        "dx": float(settings["mesh.spacing_m"]) * 1e6,
        "times": [float(t) for t in stamps], "channels": CHANNELS,
        "data": base64.b64encode(gzip.compress(data.tobytes(), mtime=0)).decode("ascii"),
        "values": int(data.size), "curves": curves.tolist(), "summary": summary,
        "settings": settings,
        "maxima": data.max(axis=(0, 1)).tolist(),
        "front": front_rows,
        "initial_inner_volume_um3": float(initial_inner_volume)*volume_um3,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", default="baseline_population")
    parser.add_argument("--coupled", default="porous_population")
    parser.add_argument("--morphology", action="store_true")
    parser.add_argument("--validation", nargs="*", default=[])
    parser.add_argument("--mechanism-audit", action="store_true")
    args = parser.parse_args()
    study = ROOT / "cases/carnallite_coating/output"
    main_cases = [dataset(study / args.coupled, "允许石盐、钾石盐沉淀"),
                  dataset(study / args.baseline, "只有光卤石溶解")]
    np.testing.assert_allclose(main_cases[0]["times"], main_cases[1]["times"], atol=1e-8, rtol=0)
    assert main_cases[0]["nx"] == main_cases[1]["nx"]
    validation = [inspect(study / name)[0] for name in args.validation]
    payload = {"cases": main_cases, "validation": validation}
    if args.mechanism_audit:
        for case in main_cases:
            if float(case["settings"].get("reactions.initial_coating.volume_m3", 0)) != 0:
                raise ValueError("A fresh-front audit cannot start with an initial coating")
            if case["front"][0]["volume_um3_Halite"] != 0 or case["front"][0]["volume_um3_Sylvite"] != 0:
                raise ValueError("A fresh-front audit must start without secondary minerals")
        np.testing.assert_equal(main_cases[0]["settings"]["mesh.spacing_m"], main_cases[1]["settings"]["mesh.spacing_m"])
        np.testing.assert_array_equal(
            np.load(study / args.coupled / "initial_geometry.npz")["solid_fractions"],
            np.load(study / args.baseline / "initial_geometry.npz")["solid_fractions"])
        for name in ("Halite", "Sylvite"):
            if main_cases[0]["front"][-1][f"inner_births_{name}"] <= 0:
                raise ValueError(f"No demonstrated fresh-front birth of {name}")
            if main_cases[1]["front"][-1][f"volume_um3_{name}"] != 0:
                raise ValueError("The dissolution-only control created a product")
        payload["mechanismAudit"] = True
    if args.morphology:
        with (study / "morphology_audit.csv").open() as stream:
            audit = list(csv.DictReader(stream))
        selected = [row for row in audit if float(row["time_s"]) == main_cases[0]["times"][-1]]
        assert any(row["run"] == args.coupled for row in selected)
        payload["morphology"] = {
            "final": selected,
            "legacySummary": inspect(study / next(
                row["run"] for row in selected if row["representation"] == "plic"))[0],
            "image": "data:image/png;base64," + base64.b64encode(
                (study / "morphology_comparison.png").read_bytes()).decode("ascii"),
        }
    template = Path(__file__).with_name("carnallite_animation.html").read_text(encoding="utf-8")
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    html = template.replace("__SIMULATION_DATA__", encoded.replace("</", "<\\/"))
    output = study / "carnallite_experiment.html"
    output.write_text(html, encoding="utf-8")
    print(f"{output} ({output.stat().st_size:,} bytes); {len(main_cases[0]['times'])} frames per case")


if __name__ == "__main__":
    main()
