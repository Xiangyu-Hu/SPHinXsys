#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Compare CPU/device histories and render the device particle fields.

Requires numpy, matplotlib and vtk. On Linux, use xvfb-run if VTK needs a display.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cpu', type=Path)
    parser.add_argument('device', type=Path)
    parser.add_argument('--output', type=Path, default=Path('.'))
    parser.add_argument('--labels', nargs=2, default=['CPU', 'GPU'])
    args = parser.parse_args()
    cpu_metadata = json.loads((args.cpu / 'case.json').read_text())
    device_metadata = json.loads((args.device / 'case.json').read_text())
    for key in ['spacing', 'speed', 'cfl', 'end_time', 'particles']:
        if not np.isclose(cpu_metadata[key], device_metadata[key], rtol=1e-6, atol=0):
            raise ValueError(f'CPU and device cases differ in {key}')
    args.output.mkdir(parents=True, exist_ok=True)
    end_us = float(f'{cpu_metadata["end_time"] * 1e6:.6g}')
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.1))
    for case, label, color, style in [(args.cpu, args.labels[0], '#3564a5', '-'),
                                    (args.device, args.labels[1], '#ba413e', '--')]:
        data = np.genfromtxt(case / 'history.csv', delimiter=',', names=True)
        for ax, name, scale, ylabel in zip(axes,
                ['force_z', 'kinetic_energy', 'mean_damage'], [1e-3, 1, 1],
                ['Contact force / kN', 'Kinetic energy / J', 'Mean damage']):
            ax.plot(data['time'] * 1e6, data[name] * scale, style,
                    color=color, lw=1.5, label=label)
            ax.set(xlabel='Time / µs', ylabel=ylabel, xlim=(0, end_us))
            ax.spines[['top', 'right']].set_visible(False)
            ax.grid(alpha=.15)
    axes[0].legend(frameon=False, fontsize=9)
    fig.tight_layout(w_pad=2)
    fig.savefig(args.output / 'response.png', dpi=220, facecolor='white',
                metadata={'Software': None})
    plt.close(fig)

    # Reuse the CPU example's camera, particle rendering and fixed color scales.
    source = Path(__file__).resolve().parents[3] / '3d_examples/test_3d_taylor_bar_hjc/plot.py'
    spec = importlib.util.spec_from_file_location('hjc_plot', source)
    renderer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(renderer)
    renderer.contours(args.device, args.output, device_metadata['spacing'], [20, 40, 60])


if __name__ == '__main__':
    main()
