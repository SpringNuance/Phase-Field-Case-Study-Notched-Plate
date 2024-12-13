# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023.HF4 replay file
# Internal Version: 2023_07_21-20.45.57 RELr425 183702
# Run by nguyenb5 on Wed Aug 14 13:18:17 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=114.666664123535, 
    height=138.495361328125)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='phase_field_plate_uel_0p5_wtppm.odb')
#: Model: C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/benchmark_Emilio/H2_0p5_wtppm/phase_field_plate_uel_0p5_wtppm.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=NONE)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_SH', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.89115, 
    farPlane=5.10885, width=2.22079, height=1.0519, viewOffsetX=0.134455, 
    viewOffsetY=-0.0489944)
session.viewports['Viewport: 1'].animationController.setValues(
    animationType=TIME_HISTORY)
session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
session.imageAnimationOptions.setValues(vpDecorations=ON, vpBackground=OFF, 
    compass=OFF)
session.writeImageAnimation(
    fileName='C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/benchmark_Emilio/H2_0p5_wtppm/phase_field', 
    format=AVI, canvasObjects=(session.viewports['Viewport: 1'], ))
