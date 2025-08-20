import paraview.simple as pv

reader = pv.OpenDataFile('/Users/gbozzola/repos/palace/examples/transmon/postpro/paraview/eigenmode/eigenmode.pvd')

slice_filter = pv.SliceWithPlane(Input=reader)
# For some reason, we have to slightly offset from the xy plane to get a good rendering.
slice_filter.PlaneType.Origin = [0, 0, 1e-3]
slice_filter.PlaneType.Normal = [0, 0, 1]

calc = pv.Calculator(Input=slice_filter)
calc.Function = 'log10(U_m)'
calc.ResultArrayName = 'log_U_m'

layout = pv.CreateLayout('Comparative View')
layout.SplitHorizontal(0, 0.5)

# Frame 0
view1 = pv.CreateView('RenderView')
view1.ViewTime = 0
view1.ViewSize = [500, 600]
layout.AssignView(1, view1)
display1 = pv.Show(calc, view1)
pv.ColorBy(display1, ('POINTS', 'log_U_m'))

# Frame 1
view2 = pv.CreateView('RenderView')
view2.ViewTime = 1
view2.ViewSize = [500, 600]
layout.AssignView(2, view2)
display2 = pv.Show(calc, view2)
pv.ColorBy(display2, ('POINTS', 'log_U_m'))

ctf = pv.GetColorTransferFunction('log_U_m')
colorbar = pv.GetScalarBar(ctf, view1)
colorbar.Title = 'U_m (J/m³)'
colorbar.ComponentTitle = ''

# Add titles
text1 = pv.Text()
text1.Text = 'First Mode'
text_display1 = pv.Show(text1, view1)
text_display1.WindowLocation = 'Upper Center'

text2 = pv.Text()
text2.Text = 'Second Mode'
text_display2 = pv.Show(text2, view2)
text_display2.WindowLocation = 'Upper Center'

# Synchronize cameras
view1.ResetCamera()
view2.CameraPosition = view1.CameraPosition
view2.CameraFocalPoint = view1.CameraFocalPoint
view2.CameraViewUp = view1.CameraViewUp
view2.CameraParallelScale = view1.CameraParallelScale

pv.Render()
pv.SaveScreenshot('slice_comparison.png', layout, ImageResolution=[1000, 600])
