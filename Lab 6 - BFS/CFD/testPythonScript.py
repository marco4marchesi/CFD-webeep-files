# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
bfs_fdsvtu = FindSource('bfs_fds.vtu')

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=bfs_fdsvtu)

# Properties modified on plotOverLine1
plotOverLine1.Point1 = [-0.8, 0.0, 0.0]
plotOverLine1.Point2 = [-0.8, 0.04, 0.0]

# show data in view
plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1, 'XYChartRepresentation')

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView1, layout=layout1, hint=0)

# Properties modified on plotOverLine1Display_1
plotOverLine1Display_1.SeriesPlotCorner = ['Density', '0', 'Eddy_Viscosity', '0', 'Heat_Flux', '0', 'Laminar_Viscosity', '0', 'Omega', '0', 'Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Pressure', '0', 'Pressure_Coefficient', '0', 'Residual_Omega', '0', 'Residual_Pressure', '0', 'Residual_TKE', '0', 'Residual_Velocity_Magnitude', '0', 'Residual_Velocity_X', '0', 'Residual_Velocity_Y', '0', 'Residual_Velocity_Z', '0', 'Skin_Friction_Coefficient_Magnitude', '0', 'Skin_Friction_Coefficient_X', '0', 'Skin_Friction_Coefficient_Y', '0', 'Skin_Friction_Coefficient_Z', '0', 'Turb_Kin_Energy', '0', 'Velocity_Magnitude', '0', 'Velocity_X', '0', 'Velocity_Y', '0', 'Velocity_Z', '0', 'Y_Plus', '0', 'arc_length', '0', 'vtkValidPointMask', '0']
plotOverLine1Display_1.SeriesLineStyle = ['Density', '1', 'Eddy_Viscosity', '1', 'Heat_Flux', '1', 'Laminar_Viscosity', '1', 'Omega', '1', 'Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Pressure', '1', 'Pressure_Coefficient', '1', 'Residual_Omega', '1', 'Residual_Pressure', '1', 'Residual_TKE', '1', 'Residual_Velocity_Magnitude', '1', 'Residual_Velocity_X', '1', 'Residual_Velocity_Y', '1', 'Residual_Velocity_Z', '1', 'Skin_Friction_Coefficient_Magnitude', '1', 'Skin_Friction_Coefficient_X', '1', 'Skin_Friction_Coefficient_Y', '1', 'Skin_Friction_Coefficient_Z', '1', 'Turb_Kin_Energy', '1', 'Velocity_Magnitude', '1', 'Velocity_X', '1', 'Velocity_Y', '1', 'Velocity_Z', '1', 'Y_Plus', '1', 'arc_length', '1', 'vtkValidPointMask', '1']
plotOverLine1Display_1.SeriesLineThickness = ['Density', '2', 'Eddy_Viscosity', '2', 'Heat_Flux', '2', 'Laminar_Viscosity', '2', 'Omega', '2', 'Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Pressure', '2', 'Pressure_Coefficient', '2', 'Residual_Omega', '2', 'Residual_Pressure', '2', 'Residual_TKE', '2', 'Residual_Velocity_Magnitude', '2', 'Residual_Velocity_X', '2', 'Residual_Velocity_Y', '2', 'Residual_Velocity_Z', '2', 'Skin_Friction_Coefficient_Magnitude', '2', 'Skin_Friction_Coefficient_X', '2', 'Skin_Friction_Coefficient_Y', '2', 'Skin_Friction_Coefficient_Z', '2', 'Turb_Kin_Energy', '2', 'Velocity_Magnitude', '2', 'Velocity_X', '2', 'Velocity_Y', '2', 'Velocity_Z', '2', 'Y_Plus', '2', 'arc_length', '2', 'vtkValidPointMask', '2']
plotOverLine1Display_1.SeriesMarkerStyle = ['Density', '0', 'Eddy_Viscosity', '0', 'Heat_Flux', '0', 'Laminar_Viscosity', '0', 'Omega', '0', 'Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Pressure', '0', 'Pressure_Coefficient', '0', 'Residual_Omega', '0', 'Residual_Pressure', '0', 'Residual_TKE', '0', 'Residual_Velocity_Magnitude', '0', 'Residual_Velocity_X', '0', 'Residual_Velocity_Y', '0', 'Residual_Velocity_Z', '0', 'Skin_Friction_Coefficient_Magnitude', '0', 'Skin_Friction_Coefficient_X', '0', 'Skin_Friction_Coefficient_Y', '0', 'Skin_Friction_Coefficient_Z', '0', 'Turb_Kin_Energy', '0', 'Velocity_Magnitude', '0', 'Velocity_X', '0', 'Velocity_Y', '0', 'Velocity_Z', '0', 'Y_Plus', '0', 'arc_length', '0', 'vtkValidPointMask', '0']
plotOverLine1Display_1.SeriesMarkerSize = ['Density', '4', 'Eddy_Viscosity', '4', 'Heat_Flux', '4', 'Laminar_Viscosity', '4', 'Omega', '4', 'Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Pressure', '4', 'Pressure_Coefficient', '4', 'Residual_Omega', '4', 'Residual_Pressure', '4', 'Residual_TKE', '4', 'Residual_Velocity_Magnitude', '4', 'Residual_Velocity_X', '4', 'Residual_Velocity_Y', '4', 'Residual_Velocity_Z', '4', 'Skin_Friction_Coefficient_Magnitude', '4', 'Skin_Friction_Coefficient_X', '4', 'Skin_Friction_Coefficient_Y', '4', 'Skin_Friction_Coefficient_Z', '4', 'Turb_Kin_Energy', '4', 'Velocity_Magnitude', '4', 'Velocity_X', '4', 'Velocity_Y', '4', 'Velocity_Z', '4', 'Y_Plus', '4', 'arc_length', '4', 'vtkValidPointMask', '4']

# update the view to ensure updated data information
renderView1.Update()

# save data
SaveData('C:/Users/marco/Desktop/UNI/2 MAGISTRALE/CFD/CFD-webeep-files/Lab 6 - BFS/CFD/plotOverLine_m0p8.csv', proxy=plotOverLine1, PointDataArrays=['Density', 'Eddy_Viscosity', 'Heat_Flux', 'Laminar_Viscosity', 'Omega', 'Pressure', 'Pressure_Coefficient', 'Residual_Omega', 'Residual_Pressure', 'Residual_TKE', 'Residual_Velocity', 'Skin_Friction_Coefficient', 'Turb_Kin_Energy', 'Velocity', 'Y_Plus', 'arc_length', 'vtkValidPointMask'])
