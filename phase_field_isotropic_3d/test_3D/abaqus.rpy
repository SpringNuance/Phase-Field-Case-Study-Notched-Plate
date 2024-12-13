# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023.HF4 replay file
# Internal Version: 2023_07_21-20.45.57 RELr425 183702
# Run by nguyenb5 on Tue Sep 17 21:04:16 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=90.2552032470703, 
    height=132.6875)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='phase_field_3D_millimeter_UEL.odb')
#: Model: C:/LocalUserData/User-data/nguyenb5/phase_field_isotropic_3d/test_3D/phase_field_3D_millimeter_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       10
#: Number of Node Sets:          8
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.28501, 
    farPlane=3.5228, width=0.194515, height=0.0910742, viewOffsetX=-0.0294276, 
    viewOffsetY=0.0117957)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR15_CL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR14_SIG_H', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.2928, 
    farPlane=3.51501, width=0.146561, height=0.0686215, 
    viewOffsetX=-0.00832411, viewOffsetY=-0.0168381)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR13_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.28891, 
    farPlane=3.5189, width=0.170504, height=0.0798319, viewOffsetX=0.236568, 
    viewOffsetY=-0.176165)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR15_CL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.25395, 
    farPlane=3.55386, width=0.58446, height=0.273651, viewOffsetX=-0.00274149, 
    viewOffsetY=0.0283547)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.14712, 
    farPlane=3.66069, width=1.04642, height=0.489948, viewOffsetX=0.116248, 
    viewOffsetY=0.165081)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=500 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.22741, 
    farPlane=3.5804, width=0.549614, height=0.257335, viewOffsetX=0.0525312, 
    viewOffsetY=0.306016)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-1', frame=192)
