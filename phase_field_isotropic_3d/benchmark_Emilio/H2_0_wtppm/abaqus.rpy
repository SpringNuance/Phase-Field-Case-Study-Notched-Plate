# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023.HF4 replay file
# Internal Version: 2023_07_21-20.45.57 RELr425 183702
# Run by nguyenb5 on Tue Aug 13 17:11:57 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=114.666664123535, 
    height=112.359954833984)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='phase_field_plate_uel_0_wtppm.odb')
#: Model: C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/benchmark_Emilio/H2_0_wtppm/phase_field_plate_uel_0_wtppm.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.85364, 
    farPlane=5.14636, width=2.63703, height=0.968016, viewOffsetX=-0.0183845, 
    viewOffsetY=-0.0532835)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.81658, 
    farPlane=5.18342, width=2.94566, height=1.08131, viewOffsetX=0.037331, 
    viewOffsetY=-0.0158458)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_SH', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].animationController.setValues(
    animationType=TIME_HISTORY)
session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.00503, 
    farPlane=4.99497, width=1.94727, height=0.714813, viewOffsetX=0.067697, 
    viewOffsetY=0.0265539)
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.showLastFrame()
session.viewports['Viewport: 1'].animationController.showLastFrame()
session.viewports['Viewport: 1'].animationController.showLastFrame()
session.viewports['Viewport: 1'].animationController.showLastFrame()
session.viewports['Viewport: 1'].animationController.showLastFrame()
session.viewports['Viewport: 1'].animationController.showLastFrame()
session.viewports['Viewport: 1'].animationController.showLastFrame()
session.viewports['Viewport: 1'].animationController.setValues(
    animationType=NONE)
session.odbs['C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/benchmark_Emilio/H2_0_wtppm/phase_field_plate_uel_0_wtppm.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/benchmark_Emilio/H2_0_wtppm/phase_field_plate_uel_0_wtppm.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/benchmark_Emilio/H2_0_wtppm/phase_field_plate_uel_0_wtppm.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.83729, 
    farPlane=5.16271, width=2.78928, height=1.0239, viewOffsetX=0.167096, 
    viewOffsetY=-0.0411649)
session.viewports['Viewport: 1'].animationController.setValues(
    animationType=TIME_HISTORY)
session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_SH', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.90272, 
    farPlane=5.09728, width=2.68239, height=0.984665, viewOffsetX=0.170833, 
    viewOffsetY=-0.0453394)
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports['Viewport: 1'].animationController.setValues(
    animationType=TIME_HISTORY)
session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.81689, 
    farPlane=5.18311, width=2.60308, height=1.44874, viewOffsetX=0.25717, 
    viewOffsetY=-0.2539)
session.imageAnimationOptions.setValues(vpDecorations=ON, vpBackground=OFF, 
    compass=OFF, timeScale=1, frameRate=50)
session.writeImageAnimation(fileName='phasefield', format=AVI, canvasObjects=(
    session.viewports['Viewport: 1'], ))
