# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# # README
#
# This script generates a visualization of cpw2d boundary mode results for
# documentation.
#
# Assuming that you have already run the thick_impedance simulation and have
# ParaView installed, run from the examples/cpw2d/ directory:
# ```bash
# pvpython plot_cpw2d_docs.py
# ```
#
# This creates `cpw2d-1.png` showing the Ex component of the second propagating
# mode.

import paraview.simple as pv

# Read mode 2 data directly (Cycle000002 = second eigenmode).
reader = pv.OpenDataFile(
    'postpro/thick_impedance/paraview/boundarymode/Cycle000002/data.pvtu'
)

calc = pv.Calculator(Input=reader)
calc.ResultArrayName = 'Ex'
calc.Function = 'E_real_X'

view = pv.CreateRenderView()
view.ViewSize = [1260, 600]
view.Background = [0.82, 0.82, 0.76]

display = pv.Show(calc, view)
display.SetRepresentationType('Surface')
pv.ColorBy(display, ('POINTS', 'Ex'))

ctf = pv.GetColorTransferFunction('Ex')
ctf.ApplyPreset('Cool to Warm', True)
ctf.RescaleTransferFunction(-4e6, 4e6)

colorbar = pv.GetScalarBar(ctf, view)
colorbar.Title = 'Ex (V/m)'
colorbar.ComponentTitle = ''
colorbar.Visibility = 1
colorbar.Orientation = 'Horizontal'
colorbar.WindowLocation = 'Lower Center'

# 2D top-down view centered on the CPW trace and gaps.
pv.Render()
view.ResetCamera()
view.CameraParallelProjection = 1
view.CameraPosition = [511, 0, 100]
view.CameraFocalPoint = [511, 0, 0]
view.CameraViewUp = [0, 1, 0]
view.CameraParallelScale = 7
view.OrientationAxesVisibility = 0

pv.Render()
pv.SaveScreenshot('cpw2d-1.png', view, ImageResolution=[1260, 600])
