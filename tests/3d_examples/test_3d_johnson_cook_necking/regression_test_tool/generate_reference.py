"""Build standard SPHinXsys regression XML from completed SPH runs.

Only known 0.2 mm / 80 us / four-interval cases are accepted. The variance is
the larger of the observed squared deviation and the framework's 1% scale floor.
"""
from pathlib import Path
import argparse
import json
import xml.etree.ElementTree as ET

import numpy as np


def samples(run):
    config = json.loads((run / "case.json").read_text())
    if not (np.isclose(config["spacing_m"], 2e-4, rtol=1e-6) and
            np.isclose(config["end_time_s"], 8e-5, rtol=1e-6) and
            config["frames"] == 4 and config["particles"] == 1648):
        raise ValueError("Use the fixed 1648-particle regression case")
    history = np.genfromtxt(run / "history.csv", delimiter=",", names=True)
    energies, damages = [], []
    for frame in range(5):
        target = frame * 20e-6
        index = np.argmin(np.abs(history["time_s"] - target))
        if abs(history["time_s"][index] - target) > 1e-9:
            raise ValueError("Missing fixed-time history sample")
        paths = list((run / "output").glob("NeckingRod*%010d.vtp" % frame))
        if len(paths) != 1:
            raise ValueError("Record all five particle snapshots to generate references")
        root = ET.parse(paths[0])
        field = next(node for node in root.iter("DataArray") if node.get("Name") == "JCDamage")
        values = np.fromstring(field.text or "", sep=" ")
        if len(values) != 1648 or not np.all(np.isfinite(values)):
            raise ValueError("Invalid damage snapshot")
        energies.append(history["kinetic_J"][index])
        damages.append(values.mean())
    return energies, damages


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if len(args.run) < 2:
        parser.error("Supply at least two completed runs")
    results = [samples(run) for run in args.run]
    args.output.mkdir(parents=True, exist_ok=True)
    for column, quantity in enumerate(["TotalKineticEnergy", "TotalJCDamage"]):
        data = np.asarray([result[column] for result in results])
        mean = data.mean(axis=0)
        variance = np.maximum(np.max((data - mean) ** 2, axis=0), (mean * 0.01) ** 2)
        root = ET.Element("meanvariance")
        for tag, values in [("Mean_Element", mean), ("Variance_Element", variance)]:
            parent = ET.SubElement(root, tag)
            for index, value in enumerate(values):
                ET.SubElement(parent, "Snapshot_%d" % index, {quantity + "_0": "%.17g" % value})
        ET.indent(root, space="    ")
        name = "NeckingRod_" + quantity
        ET.ElementTree(root).write(args.output / (name + "_ensemble_averaged_mean_variance.xml"), encoding="utf-8", xml_declaration=True)
        (args.output / (name + "_runtimes.dat")).write_text("true\n%d\n0\n" % len(results))
    print(args.output)


if __name__ == "__main__":
    main()
