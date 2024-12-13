# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2023.HF4 replay file
# Internal Version: 2023_07_21-20.45.57 RELr425 183702
# Run by nguyenb5 on Sun Sep  1 09:43:35 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=90.2552032470703, 
    height=154.80207824707)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('damage_model_3D.cae')
#: The model database "C:\LocalUserData\User-data\nguyenb5\phase_field_isotropic_3d\process_input_file_to_UEL_millimeter\damage_model_3D.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-2D-millimeter-ref'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p1 = mdb.models['Model-3D-millimeter'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.24353, 
    farPlane=3.57254, width=0.562754, height=0.314591, viewOffsetX=0.133613, 
    viewOffsetY=-0.104369)
p1 = mdb.models['Model-2D-millimeter-ref-Copy'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
del mdb.models['Model-2D-millimeter-ref-Copy']
p = mdb.models['Model-2D-millimeter-ref'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.23681, 
    farPlane=3.56977, width=0.546985, height=0.305776, viewOffsetX=0.171534, 
    viewOffsetY=-0.109897)
p1 = mdb.models['Model-3D-millimeter'].parts['plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON, mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
a = mdb.models['Model-3D-millimeter'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON, optimizationTasks=OFF, 
    geometricRestrictions=OFF, stopConditions=OFF)
mdb.models['Model-3D-millimeter'].fieldOutputRequests['F-Output-1'].setValues(
    variables=('RF', 'SDV', 'U'), frequency=1)
mdb.save()
#: The model database has been saved to "C:\LocalUserData\User-data\nguyenb5\phase_field_isotropic_3d\process_input_file_to_UEL_millimeter\damage_model_3D.cae".
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.95826, 
    farPlane=3.85781, width=2.20869, height=0.986586, viewOffsetX=0.0494762, 
    viewOffsetY=-0.02306)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
    bcs=OFF, predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
mdb.save()
#: The model database has been saved to "C:\LocalUserData\User-data\nguyenb5\phase_field_isotropic_3d\process_input_file_to_UEL_millimeter\damage_model_3D.cae".
