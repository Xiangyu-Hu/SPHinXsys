#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Plot impact histories and particle contours. Requires numpy, matplotlib, vtk.

Example: python plot.py impact --compare impact_half_dt
Use xvfb-run on Linux if VTK requires a display for offscreen rendering.
"""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read_vtp(path):
    tree = ET.parse(path)
    arrays = {}
    for node in tree.findall('.//DataArray'):
        if node.get('format') != 'ascii':
            raise ValueError('This reader expects the default ASCII VTP output')
        values = np.fromstring(node.text or '', sep=' ')
        components = int(node.get('NumberOfComponents', '1'))
        arrays[node.get('Name')] = values.reshape(-1, components) if components > 1 else values
    return arrays


def histories(case, comparison, output):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.1))
    for folder, color, style in [(case, '#3564a5', '-'), (comparison, '#ba413e', '--')]:
        if folder is None:
            continue
        metadata = json.loads((folder/'case.json').read_text())
        label = f'CFL {metadata["cfl"]:g}'
        data = np.genfromtxt(folder/'history.csv', delimiter=',', names=True)
        for ax, name, scale, ylabel in zip(axes, ['force_z', 'mean_velocity_z', 'mean_damage'],
                [1e-3, 1, 1], ['Contact force / kN', 'Axial velocity / m s$^{-1}$', 'Mean damage']):
            ax.plot(data['time']*1e6, data[name]*scale, style, color=color, lw=1.5, label=label)
            ax.set(xlabel='Time / µs', ylabel=ylabel, xlim=(0, data['time'][-1]*1e6))
            ax.spines[['top', 'right']].set_visible(False)
            ax.grid(alpha=.15)
    axes[0].legend(frameon=False, fontsize=9)
    fig.tight_layout(w_pad=2)
    fig.savefig(output/'response.png', dpi=220, facecolor='white')
    plt.close(fig)


def contours(case, output, spacing, times):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    frames = [(path, read_vtp(path)) for path in sorted((case/'output').glob('Concrete_ite_*.vtp'))]
    if not frames:
        raise ValueError('No Concrete VTP files; enable state_recording')
    snapshots = []
    for t in times:
        _, data = min(frames, key=lambda item: abs(item[1]['TimeValue'][0]*1e6-t))
        if abs(data['TimeValue'][0]*1e6-t) > 1e-3:
            raise ValueError(f'{t} µs is not a saved frame')
        snapshots.append(data)
    # Fixed scales across time. Particle values are never interpolated or smoothed.
    qmax = max(frame['VonMisesStress'].max() for frame in snapshots)*1e-6
    qlimit = np.ceil(qmax/50)*50
    width, height = 500, 570
    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(width*len(times)+110, height*2)
    window.SetMultiSamples(4)
    keep = []
    for row, (name, scale, upper, label) in enumerate([
            ('HJCDamage', 1, 1, 'Damage'),
            ('VonMisesStress', 1e-6, qlimit, 'q / MPa')]):
        lut = vtk.vtkColorTransferFunction()
        lut.SetColorSpaceToDiverging()
        lut.AddRGBPoint(0, .230, .299, .754)
        lut.AddRGBPoint(upper/2, .865, .865, .865)
        lut.AddRGBPoint(upper, .706, .016, .150)
        for col, data in enumerate(snapshots):
            renderer = vtk.vtkRenderer()
            renderer.SetBackground(.32, .355, .43)
            renderer.SetViewport(col*width/(width*len(times)+110), (1-row)*.5,
                                 (col+1)*width/(width*len(times)+110), (2-row)*.5)
            window.AddRenderer(renderer)
            pos = data['Position']*1000
            mask = pos[:, 1] >= 0
            points = vtk.vtkPoints()
            points.SetData(numpy_to_vtk(np.ascontiguousarray(pos[mask]), deep=True))
            cloud = vtk.vtkPolyData()
            cloud.SetPoints(points)
            values = numpy_to_vtk(np.ascontiguousarray(data[name][mask]*scale), deep=True)
            values.SetName(name)
            cloud.GetPointData().SetScalars(values)
            sphere = vtk.vtkSphereSource()
            sphere.SetRadius(spacing*1000*.46)
            sphere.SetThetaResolution(14)
            sphere.SetPhiResolution(10)
            mapper = vtk.vtkGlyph3DMapper()
            mapper.SetInputData(cloud)
            mapper.SetSourceConnection(sphere.GetOutputPort())
            mapper.ScalingOff()
            mapper.SetLookupTable(lut)
            mapper.SetScalarRange(0, upper)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetAmbient(.25)
            actor.GetProperty().SetDiffuse(.75)
            actor.GetProperty().SetSpecular(.25)
            actor.GetProperty().SetSpecularPower(22)
            renderer.AddActor(actor)
            wall = vtk.vtkCubeSource()
            wall.SetBounds(-7, 7, 0, 5, -1.5, 0)
            wm = vtk.vtkPolyDataMapper()
            wm.SetInputConnection(wall.GetOutputPort())
            wa = vtk.vtkActor()
            wa.SetMapper(wm)
            wa.GetProperty().SetColor(.78, .79, .81)
            renderer.AddActor(wa)
            camera = renderer.GetActiveCamera()
            camera.SetPosition(30, -60, 28)
            camera.SetFocalPoint(0, 0, 9.5)
            camera.SetViewUp(0, 0, 1)
            camera.ParallelProjectionOn()
            camera.SetParallelScale(12.6)
            renderer.ResetCameraClippingRange()
            title = vtk.vtkTextActor()
            title.SetInput(f'{times[col]:g} µs')
            title.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
            title.SetPosition(.42, .94)
            title.GetTextProperty().SetFontSize(22)
            title.GetTextProperty().SetColor(1, 1, 1)
            if row == 0:
                renderer.AddActor2D(title)
            keep.extend([cloud, points, values, sphere, mapper, actor, wall, wm, wa, title])
        bar_renderer = vtk.vtkRenderer()
        bar_renderer.SetBackground(.32, .355, .43)
        bar_renderer.SetViewport(width*len(times)/(width*len(times)+110), (1-row)*.5, 1, (2-row)*.5)
        window.AddRenderer(bar_renderer)
        bar = vtk.vtkScalarBarActor()
        bar.SetLookupTable(lut)
        bar.SetTitle(label)
        bar.SetNumberOfLabels(3)
        bar.SetLabelFormat('%.3g')
        bar.SetPosition(.02, .25)
        bar.SetWidth(.85)
        bar.SetHeight(.52)
        bar.UnconstrainedFontSizeOn()
        bar.GetLabelTextProperty().SetFontSize(16)
        bar.GetTitleTextProperty().SetFontSize(16)
        bar.GetTitleTextProperty().SetBold(False)
        bar.GetTitleTextProperty().SetItalic(False)
        bar.GetLabelTextProperty().SetItalic(False)
        bar_renderer.AddActor2D(bar)
        keep.extend([bar, lut])
    window.Render()
    capture = vtk.vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.SetScale(1)
    capture.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(str(output/'impact.png'))
    writer.SetInputConnection(capture.GetOutputPort())
    writer.Write()
    window.Finalize()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    parser.add_argument('--compare', type=Path)
    parser.add_argument('--output', type=Path, default=Path('.'))
    parser.add_argument('--times', nargs='+', type=float, default=[20, 40, 60])
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    histories(args.case, args.compare, args.output)
    metadata = json.loads((args.case/'case.json').read_text())
    contours(args.case, args.output, metadata['spacing'], args.times)
