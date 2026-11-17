# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# # README
#
# This script generates the field visualizations of the superconducting coaxial example for
# documentation.
#
# It reads the refined superconducting run written by `generate_kinetic_inductance_fields`
# in `coaxial_kinetic_inductance.jl`, which also calls this script. To run it on its own,
# assuming that the simulation output exists and ParaView is installed, run from the
# examples/coaxial/ directory:
# ```bash
# pvpython plot_magnetostatic_docs.py
# ```
#
# This creates `coaxial-5.png`, the sheet current on the conductors near the driven end,
# and `coaxial-6.png`, the magnetic flux density on the mid-length cross-section.

import paraview.simple as pv

paraview_dir = 'postpro/magnetostatic_superconductor_fields/paraview'


def make_view(size):
    view = pv.CreateView('RenderView')
    view.ViewSize = size
    view.ViewTime = 0.0
    view.OrientationAxesVisibility = 1
    view.UseColorPaletteForBackground = 0
    view.Background = [1, 1, 1]
    return view


def show_scalar_bar(display, view, array, title):
    bar = pv.GetScalarBar(pv.GetColorTransferFunction(array), view)
    bar.Title = title
    bar.ComponentTitle = ''
    bar.TitleColor = [0, 0, 0]
    bar.LabelColor = [0, 0, 0]
    bar.TitleFontSize = 18
    bar.LabelFontSize = 16
    bar.RangeLabelFormat = '%.3g'
    display.SetScalarBarVisibility(view, True)


# Sheet current on the driven end cap (attribute 3) and the first 12 μm of the inner and
# outer conductors (attribute 2), with the front half (y < 0) cut away.
boundary = pv.PVDReader(
    FileName=f'{paraview_dir}/magnetostatic_boundary/magnetostatic_boundary.pvd'
)
conductors = pv.Threshold(Input=boundary)
conductors.Scalars = ['CELLS', 'attribute']
conductors.LowerThreshold = 2
conductors.UpperThreshold = 3
near_end = pv.Clip(Input=conductors)
near_end.ClipType = 'Plane'
near_end.ClipType.Origin = [0, 0, 12]
near_end.ClipType.Normal = [0, 0, 1]
near_end.Invert = 1
cutaway = pv.Clip(Input=near_end)
cutaway.ClipType = 'Plane'
cutaway.ClipType.Origin = [0, 0, 0]
cutaway.ClipType.Normal = [0, 1, 0]
cutaway.Invert = 0
cutaway.UpdatePipeline(0.0)

view = make_view([1600, 1000])
display = pv.Show(cutaway, view)
pv.ColorBy(display, ('POINTS', 'J_s', 'Magnitude'))
lut = pv.GetColorTransferFunction('J_s')
lut.ApplyPreset('Viridis', True)
lut.AutomaticRescaleRangeMode = 'Never'
lut.RescaleTransferFunction(*cutaway.PointData['J_s'].GetRange(-1))
show_scalar_bar(display, view, 'J_s', 'Sheet current magnitude (A/m)')
view.CameraPosition = [7, -16, 11]
view.CameraFocalPoint = [0, 0, 5]
view.CameraViewUp = [0, 0, 1]
view.ResetCamera()
pv.Render(view)
pv.SaveScreenshot('coaxial-5.png', view, ImageResolution=[1600, 1000])
pv.Delete(view)

# Magnetic flux density on the mid-length cross-section, with arrows showing its azimuthal
# direction.
volume = pv.PVDReader(FileName=f'{paraview_dir}/magnetostatic/magnetostatic.pvd')
cross_section = pv.Slice(Input=volume)
cross_section.SliceType = 'Plane'
cross_section.SliceType.Origin = [0, 0, 20]
cross_section.SliceType.Normal = [0, 0, 1]
cross_section.UpdatePipeline(0.0)
arrows = pv.Glyph(Input=cross_section, GlyphType='Arrow')
arrows.OrientationArray = ['POINTS', 'B']
arrows.ScaleArray = ['POINTS', 'No scale array']
arrows.ScaleFactor = 0.22
arrows.GlyphMode = 'Every Nth Point'
arrows.Stride = 9
arrows.UpdatePipeline(0.0)

view = make_view([1300, 1100])
display = pv.Show(cross_section, view)
pv.ColorBy(display, ('POINTS', 'B', 'Magnitude'))
pv.GetColorTransferFunction('B').ApplyPreset('Inferno', True)
show_scalar_bar(display, view, 'B', 'Magnetic flux density magnitude (T)')
arrow_display = pv.Show(arrows, view)
pv.ColorBy(arrow_display, None)
arrow_display.DiffuseColor = [1, 1, 1]
view.CameraPosition = [0, 0, 60]
view.CameraFocalPoint = [0, 0, 20]
view.CameraViewUp = [0, 1, 0]
view.ResetCamera()
view.CameraParallelProjection = 1
pv.Render(view)
pv.SaveScreenshot('coaxial-6.png', view, ImageResolution=[1300, 1100])
pv.Delete(view)
