from paraview.simple import *
from paraview.servermanager import *
#from paraview import vtk
#import importlib.util

# get source
a1 = GetActiveSource() 

# plot over line
plotOverLine1 = PlotOverLine(Input=a1, Source='High Resolution Line Source')  

# select fields
passArrays1 = PassArrays(Input=plotOverLine1)
passArrays1.PointDataArrays = ['Velocity']

x_vec = [-0.8, -0.4, 0, 0.05,     0.1,    0.2,     0.3,      0.5,     1,       1.4   ] 
y_vec1 = [0,    0,   0, -0.0267, -0.0267, -0.0267, -0.0267, -0.0267, -0.0267, -0.0267]
for i in range(len(x_vec)) :
    plotOverLine1.Source.Point1 = [x_vec[i], y_vec1[i], 0 ]
    plotOverLine1.Source.Point2 = [x_vec[i], 0.04, 0 ]
    writer = CreateWriter('./bfsSimulationData_x{0}.csv',format(x_vec1[i]))
    writer.UpdatePipeline()
