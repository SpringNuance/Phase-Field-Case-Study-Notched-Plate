# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2023.HF4 replay file
# Internal Version: 2023_07_21-20.45.57 RELr425 183702
# Run by nguyenb5 on Tue Sep 17 21:07:58 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=90.2552032470703, 
    height=126.43286895752)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('damage_model_3D.cae')
#: The model database "C:\LocalUserData\User-data\nguyenb5\phase_field_isotropic_3d\benchmark_Emilio\damage_model_3D.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-2D'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p1 = mdb.models['Model-3D-mm'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p1 = mdb.models['Model-3D'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
#: 
#: Point 1: -500.E-03, 500.E-03, 20.E-03  Point 2: -500.E-03, 500.E-03, 0.
#:    Distance: 20.E-03  Components: 0., 0., -20.E-03
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.26309, 
    farPlane=3.55298, width=0.379137, height=0.178157, viewOffsetX=-0.199643, 
    viewOffsetY=0.410015)
p1 = mdb.models['Model-3D-mm'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
#: 
#: Point 1: -500.E-06, 500.E-06, 20.E-06  Point 2: -500.E-06, 500.E-06, 0.
#:    Distance: 20.E-06  Components: 0., 0., -20.E-06
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00203309, 
    farPlane=0.00378298, width=0.00180081, height=0.000846204, 
    viewOffsetX=0.000333392, viewOffsetY=0.000149075)
p1 = mdb.models['Model-2D'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
#: 
#: Point 1: -500.E-03, 500.E-03, 0.  Point 2: 500.E-03, 500.E-03, 0.
#:    Distance: 1.  Components: 1., 0., 0.
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.48565, 
    farPlane=3.1712, width=2.1148, height=0.993748, viewOffsetX=0.288104, 
    viewOffsetY=-0.00665538)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON, mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
p1 = mdb.models['Model-3D-mm'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
a = mdb.models['Model-3D-mm'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, optimizationTasks=OFF, 
    geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
#: 
#: Point 1: 500.E-06, -500.E-06, 20.E-06  Point 2: 500.E-06, 500.E-06, 20.E-06
#:    Distance: 1.E-03  Components: 0., 1.E-03, 0.
#: 
#: Point 1: 500.E-06, 500.E-06, 20.E-06  Point 2: 500.E-06, 500.E-06, 0.
#:    Distance: 20.E-06  Components: 0., 0., -20.E-06
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00205242, 
    farPlane=0.00376365, width=0.00218839, height=0.00102463, 
    viewOffsetX=0.000271352, viewOffsetY=-9.6791e-06)
a = mdb.models['Model-3D'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
#: 
#: Point 1: 500.E-03, 500.E-03, 20.E-03  Point 2: 500.E-03, -500.E-03, 20.E-03
#:    Distance: 1.  Components: 0., -1., 0.
#: 
#: Point 1: 500.E-03, -500.E-03, 20.E-03  Point 2: 500.E-03, -500.E-03, 0.
#:    Distance: 20.E-03  Components: 0., 0., -20.E-03
