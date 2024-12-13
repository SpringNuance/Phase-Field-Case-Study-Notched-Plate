# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023.HF4 replay file
# Internal Version: 2023_07_21-20.45.57 RELr425 183702
# Run by nguyenb5 on Sat Aug 24 12:59:54 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=118.47395324707, 
    height=144.303237915039)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='single_cube.odb')
#: Model: C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/reference/single_cube.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       5
#: Number of Node Sets:          9
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.39822, 
    farPlane=4.91885, width=2.44463, height=1.19723, viewOffsetX=0.186708, 
    viewOffsetY=-0.0311639)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.29895, 
    farPlane=4.96039, width=3.00153, height=1.46996, viewOffsetX=0.24574, 
    viewOffsetY=-0.096712)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.51305, 
    farPlane=4.80403, width=1.9043, height=0.932607, viewOffsetX=-0.0947971, 
    viewOffsetY=0.038596)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='COORD', outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, 'COORD1'), )
#: 
#: Node: CUBE-1.123
#:                                         1             2             3        Magnitude
#: Base coordinates:                  9.00000e-01,  9.00000e-01,  1.00000e+00,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   9.00000e-01,  9.00000e-01,  1.00000e+00,      -      
#: Deformed coordinates (scaled):     9.00000e-01,  9.00000e-01,  1.00000e+00,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.122
#:                                         1             2             3        Magnitude
#: Base coordinates:                  9.00000e-01,  1.00000e+00,  1.00000e+00,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   9.00000e-01,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (scaled):     9.00000e-01,  1.00000e+00,  1.00000e+00,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.133
#:                                         1             2             3        Magnitude
#: Base coordinates:                  9.00000e-01,  1.00000e+00,  9.00000e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   9.00000e-01,  1.00000e+00,  9.00000e-01,      -      
#: Deformed coordinates (scaled):     9.00000e-01,  1.00000e+00,  9.00000e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.1198
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00000e-01,  1.00000e-01,  1.00000e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   1.00000e-01,  1.00000e-01,  1.00000e-01,      -      
#: Deformed coordinates (scaled):     1.00000e-01,  1.00000e-01,  1.00000e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.1
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (scaled):     1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.2
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00000e+00,  9.00000e-01,  1.00000e+00,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   1.00000e+00,  9.00000e-01,  1.00000e+00,      -      
#: Deformed coordinates (scaled):     1.00000e+00,  9.00000e-01,  1.00000e+00,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.13
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00000e+00,  9.00000e-01,  9.00000e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   1.00000e+00,  9.00000e-01,  9.00000e-01,      -      
#: Deformed coordinates (scaled):     1.00000e+00,  9.00000e-01,  9.00000e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.145
#:                                         1             2             3        Magnitude
#: Base coordinates:                  9.00000e-01,  9.00000e-01,  8.00000e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   9.00000e-01,  9.00000e-01,  8.00000e-01,      -      
#: Deformed coordinates (scaled):     9.00000e-01,  9.00000e-01,  8.00000e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.12
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00000e+00,  1.00000e+00,  9.00000e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   1.00000e+00,  1.00000e+00,  9.00000e-01,      -      
#: Deformed coordinates (scaled):     1.00000e+00,  1.00000e+00,  9.00000e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: CUBE-1.12
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00000e+00,  1.00000e+00,  9.00000e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   1.00000e+00,  1.00000e+00,  9.00000e-01,      -      
#: Deformed coordinates (scaled):     1.00000e+00,  1.00000e+00,  9.00000e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Element: CUBE-1.1
#:   Type: C3D8T
#:   Material: STEEL
#:   Section: SET-1.Section-ASSEMBLY_CUBE-1_SET-1, Homogeneous Solid Section, Thickness = 0
#:   Connect: 122, 123, 134, 133, 1, 2, 13, 12
#:   COORD, COORD1 (Not averaged): 0.921132, 0.921132, 0.921132, 0.921132, 0.978868, 0.978868, 0.978868, 0.978868
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.67696, 
    farPlane=4.64012, width=0.711895, height=0.272018, viewOffsetX=0.0230391, 
    viewOffsetY=0.0617014)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='COORD', outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, 'COORD2'), )
#: 
#: Element: CUBE-1.1
#:   Type: C3D8T
#:   Material: STEEL
#:   Section: SET-1.Section-ASSEMBLY_CUBE-1_SET-1, Homogeneous Solid Section, Thickness = 0
#:   Connect: 122, 123, 134, 133, 1, 2, 13, 12
#:   COORD, COORD2 (Not averaged): 0.978868, 0.921132, 0.978868, 0.921132, 0.978868, 0.921132, 0.978868, 0.921132
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='COORD', outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, 'COORD3'), )
#: 
#: Element: CUBE-1.1
#:   Type: C3D8T
#:   Material: STEEL
#:   Section: SET-1.Section-ASSEMBLY_CUBE-1_SET-1, Homogeneous Solid Section, Thickness = 0
#:   Connect: 122, 123, 134, 133, 1, 2, 13, 12
#:   COORD, COORD3 (Not averaged): 0.978868, 0.978868, 0.921132, 0.921132, 0.978868, 0.978868, 0.921132, 0.921132
