"""Plot a completed SPH run using only its own history and ASCII VTP states.

Requires NumPy and Matplotlib. Particle output is optional: history-only runs
still produce response.png. No external reference files are read.
"""
from pathlib import Path
import argparse
import json
import re
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    return np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))


def read_time(path):
    with path.open("rb") as stream:
        for _, node in ET.iterparse(stream, events=["end"]):
            if node.tag == "DataArray":
                if node.get("Name") == "TimeValue":
                    return float((node.text or "").split()[0])
                node.clear()
    return None


def read_state(path):
    fields = {}
    required = {"Position", "InitialPosition", "JCPlasticStrain", "MaterialActive"}
    for _, node in ET.iterparse(path, events=["end"]):
        if node.tag == "DataArray":
            if node.get("Name") in required:
                if node.get("format", "ascii") != "ascii":
                    raise ValueError("Use the default ASCII state writer for this script")
                values = np.fromstring(node.text or "", sep=" ")
                components = int(node.get("NumberOfComponents", "1"))
                fields[node.get("Name")] = values.reshape(-1, components) if components > 1 else values
            node.clear()
    return fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, nargs="?", default=Path.cwd())
    args = parser.parse_args()
    run = args.run.resolve()
    config = json.loads((run / "case.json").read_text())
    data = read_csv(run / "history.csv")
    output = run / "figures"
    output.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    axes[0, 0].plot(data["extension_m"] * 1e3, data["reaction_N"] / 1e3)
    axes[0, 0].set(xlabel="Grip extension (mm)", ylabel="Right grip reaction (kN)")
    for field, axis, label in [("max_temperature_K", axes[0, 1], "Maximum recorded temperature (K)"),
                               ("max_peeq", axes[1, 0], "Maximum recorded plastic strain")]:
        axis.plot(data["time_s"] * 1e6, data[field])
        axis.set(xlabel="Time (us)", ylabel=label)
    axes[1, 1].plot(data["time_s"] * 1e6, data["failed_particles"] / config["particles"] * 100)
    axes[1, 1].set(xlabel="Time (us)", ylabel="Failed initial volume (%)")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
    fig.suptitle("Johnson-Cook tensile necking | raw SPH responses")
    fig.savefig(output / "response.png", dpi=180)
    plt.close(fig)
    paths = sorted((run / "output").glob("NeckingRod*.vtp"),
                   key=lambda path: int(re.findall(r"\d+", path.stem)[-1]))
    available = [(read_time(path), path) for path in paths]
    available = [(time, path) for time, path in available if time is not None]
    selected = []
    for microseconds in [35, 80]:
        if available:
            time, path = min(available, key=lambda entry: abs(entry[0] - microseconds * 1e-6))
            if abs(time - microseconds * 1e-6) <= 1e-9:
                selected.append((microseconds, path))
    if selected:
        width, height = 2.2 * len(selected) + 0.8, 8.3
        fig = plt.figure(figsize=(width, height))
        for column, (microseconds, path) in enumerate(selected):
            state = read_state(path)
            positions = state["Position"] * 1e3
            initial = state["InitialPosition"]
            active = (initial[:, 1] > 0) & (initial[:, 1] < config["spacing_m"]) & (state["MaterialActive"] > 0)
            axis = fig.add_axes([(0.55 + column * 2.2) / width, 0.6 / height, 1.5 / width, 6.375 / height])
            points = axis.scatter(positions[active, 0], positions[active, 2],
                                  c=state["JCPlasticStrain"][active], s=3, cmap="turbo", vmin=0, vmax=1.8, linewidths=0)
            axis.set(xlim=(-1, 1), ylim=(-4.25, 4.25), xlabel="x (mm)", ylabel="z (mm)", title="SPHinXsys\n%d us" % microseconds)
            axis.set_aspect("equal")
        colorbar = fig.add_axes([(width - 0.77) / width, 1.2 / height, 0.14 / width, 5.2 / height])
        fig.colorbar(points, cax=colorbar, label="Equivalent plastic strain")
        fig.suptitle("JC necking | active material\nInitial y=0 section\nTrue deformation x1", fontsize=10)
        fig.savefig(output / "necking.png", dpi=180)
        plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
