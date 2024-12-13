# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023.HF4 replay file
# Internal Version: 2023_07_21-20.45.57 RELr425 183702
# Run by nguyenb5 on Tue Sep 10 14:23:09 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=90.2552032470703, 
    height=126.87963104248)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='phase_field_plate_3D_UEL.odb')
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00190433, 
    farPlane=0.00390953, width=0.00259674, height=0.00102483, 
    viewOffsetX=8.17693e-05, viewOffsetY=-3.52058e-05)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM5', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM1', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM2', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM3', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM4', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00204483, 
    farPlane=0.00377123, width=0.00206314, height=0.000814241, 
    viewOffsetX=0.000145221, viewOffsetY=-5.72229e-06)
session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00212836, 
    farPlane=0.00368771, width=0.00121548, height=0.000479704, 
    viewOffsetX=8.2115e-05, viewOffsetY=0.000253814)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00220075, 
    farPlane=0.00361311, width=0.000723099, height=0.000285379, 
    viewOffsetX=-2.93607e-06, viewOffsetY=0.000309978)
odb = session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb']
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
    NODAL, ((COMPONENT, 'RF2'), )), ('U', NODAL, ((COMPONENT, 'U2'), )), ), 
    nodePick=(('PLATE-1', 135, (
    '[#5 #1ffe #0:12 #40000000 #80000001 #7ff #0:12', 
    ' #50000000 #e0000000 #1ff #0:12 #14000000 #f8000000 #7f', 
    ' #0:12 #5000000 #fe000000 #1f #0:429 #f0000000 #ffffffff', 
    ' #0:6 #3ffffff8 #0 #60 ]', )), ), )
xy1 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1']
xy2 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 3']
xy3 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 34']
xy4 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 35']
xy5 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 36']
xy6 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 37']
xy7 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 38']
xy8 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 39']
xy9 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 40']
xy10 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 41']
xy11 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 42']
xy12 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 43']
xy13 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 44']
xy14 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 45']
xy15 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 479']
xy16 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 481']
xy17 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 512']
xy18 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 513']
xy19 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 514']
xy20 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 515']
xy21 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 516']
xy22 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 517']
xy23 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 518']
xy24 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 519']
xy25 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 520']
xy26 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 521']
xy27 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 522']
xy28 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 523']
xy29 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 957']
xy30 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 959']
xy31 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 990']
xy32 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 991']
xy33 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 992']
xy34 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 993']
xy35 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 994']
xy36 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 995']
xy37 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 996']
xy38 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 997']
xy39 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 998']
xy40 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 999']
xy41 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1000']
xy42 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1001']
xy43 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1435']
xy44 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1437']
xy45 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1468']
xy46 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1469']
xy47 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1470']
xy48 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1471']
xy49 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1472']
xy50 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1473']
xy51 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1474']
xy52 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1475']
xy53 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1476']
xy54 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1477']
xy55 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1478']
xy56 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1479']
xy57 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1913']
xy58 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1915']
xy59 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1946']
xy60 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1947']
xy61 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1948']
xy62 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1949']
xy63 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1950']
xy64 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1951']
xy65 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1952']
xy66 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1953']
xy67 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1954']
xy68 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1955']
xy69 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1956']
xy70 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 1957']
xy71 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15741']
xy72 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15742']
xy73 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15743']
xy74 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15744']
xy75 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15745']
xy76 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15746']
xy77 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15747']
xy78 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15748']
xy79 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15749']
xy80 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15750']
xy81 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15751']
xy82 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15752']
xy83 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15753']
xy84 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15754']
xy85 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15755']
xy86 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15756']
xy87 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15757']
xy88 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15758']
xy89 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15759']
xy90 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15760']
xy91 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15761']
xy92 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15762']
xy93 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15763']
xy94 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15764']
xy95 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15765']
xy96 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15766']
xy97 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15767']
xy98 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15768']
xy99 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15769']
xy100 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15770']
xy101 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15771']
xy102 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15772']
xy103 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15773']
xy104 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15774']
xy105 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15775']
xy106 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15776']
xy107 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15972']
xy108 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15973']
xy109 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15974']
xy110 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15975']
xy111 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15976']
xy112 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15977']
xy113 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15978']
xy114 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15979']
xy115 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15980']
xy116 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15981']
xy117 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15982']
xy118 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15983']
xy119 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15984']
xy120 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15985']
xy121 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15986']
xy122 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15987']
xy123 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15988']
xy124 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15989']
xy125 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15990']
xy126 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15991']
xy127 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15992']
xy128 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15993']
xy129 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15994']
xy130 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15995']
xy131 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15996']
xy132 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15997']
xy133 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 15998']
xy134 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 16038']
xy135 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 16039']
xy136 = avg((xy1, xy2, xy3, xy4, xy5, xy6, xy7, xy8, xy9, xy10, xy11, xy12, 
    xy13, xy14, xy15, xy16, xy17, xy18, xy19, xy20, xy21, xy22, xy23, xy24, 
    xy25, xy26, xy27, xy28, xy29, xy30, xy31, xy32, xy33, xy34, xy35, xy36, 
    xy37, xy38, xy39, xy40, xy41, xy42, xy43, xy44, xy45, xy46, xy47, xy48, 
    xy49, xy50, xy51, xy52, xy53, xy54, xy55, xy56, xy57, xy58, xy59, xy60, 
    xy61, xy62, xy63, xy64, xy65, xy66, xy67, xy68, xy69, xy70, xy71, xy72, 
    xy73, xy74, xy75, xy76, xy77, xy78, xy79, xy80, xy81, xy82, xy83, xy84, 
    xy85, xy86, xy87, xy88, xy89, xy90, xy91, xy92, xy93, xy94, xy95, xy96, 
    xy97, xy98, xy99, xy100, xy101, xy102, xy103, xy104, xy105, xy106, xy107, 
    xy108, xy109, xy110, xy111, xy112, xy113, xy114, xy115, xy116, xy117, 
    xy118, xy119, xy120, xy121, xy122, xy123, xy124, xy125, xy126, xy127, 
    xy128, xy129, xy130, xy131, xy132, xy133, xy134, xy135))
xy136.setValues(
    sourceDescription='avg ( ( "U:U2 PI: PLATE-1 N: 1", "U:U2 PI: PLATE-1 N: 3", "U:U2 PI: PLATE-1 N: 34", "U:U2 PI: PLATE-1 N: 35", "U:U2 PI: PLATE-1 N: 36", "U:U2 PI: PLATE-1 N: 37", "U:U2 PI: PLATE-1 N: 38", "U:U2 PI: PLATE-1 N: 39", "U:U2 PI: PLATE-1 N: 40", "U:U2 PI: PLATE-1 N: 41", "U:U2 PI: PLATE-1 N: 42", "U:U2 PI: PLATE-1 N: 43", "U:U2 PI: PLATE-1 N: 44", "U:U2 PI: PLATE-1 N: 45", "U:U2 PI: PLATE-1 N: 479", "U:U2 PI: PLATE-1 N: 481", "U:U2 PI: PLATE-1 N: 512", "U:U2 PI: PLATE-1 N: 513", "U:U2 PI: PLATE-1 N: 514", "U:U2 PI: PLATE-1 N: 515", "U:U2 PI: PLATE-1 N: 516", "U:U2 PI: PLATE-1 N: 517", "U:U2 PI: PLATE-1 N: 518", "U:U2 PI: PLATE-1 N: 519", "U:U2 PI: PLATE-1 N: 520", "U:U2 PI: PLATE-1 N: 521", "U:U2 PI: PLATE-1 N: 522", "U:U2 PI: PLATE-1 N: 523", "U:U2 PI: PLATE-1 N: 957", "U:U2 PI: PLATE-1 N: 959", "U:U2 PI: PLATE-1 N: 990", "U:U2 PI: PLATE-1 N: 991", "U:U2 PI: PLATE-1 N: 992", "U:U2 PI: PLATE-1 N: 993", "U:U2 PI: PLATE-1 N: 994", "U:U2 PI: PLATE-1 N: 995", "U:U2 PI: PLATE-1 N: 996", "U:U2 PI: PLATE-1 N: 997", "U:U2 PI: PLATE-1 N: 998", "U:U2 PI: PLATE-1 N: 999", "U:U2 PI: PLATE-1 N: 1000", "U:U2 PI: PLATE-1 N: 1001", "U:U2 PI: PLATE-1 N: 1435", "U:U2 PI: PLATE-1 N: 1437", "U:U2 PI: PLATE-1 N: 1468", "U:U2 PI: PLATE-1 N: 1469", "U:U2 PI: PLATE-1 N: 1470", "U:U2 PI: PLATE-1 N: 1471", "U:U2 PI: PLATE-1 N: 1472", "U:U2 PI: PLATE-1 N: 1473", "U:U2 PI: PLATE-1 N: 1474", "U:U2 PI: PLATE-1 N: 1475", "U:U2 PI: PLATE-1 N: 1476", "U:U2 PI: PLATE-1 N: 1477", "U:U2 PI: PLATE-1 N: 1478", "U:U2 PI: PLATE-1 N: 1479", "U:U2 PI: PLATE-1 N: 1913", "U:U2 PI: PLATE-1 N: 1915", "U:U2 PI: PLATE-1 N: 1946", "U:U2 PI: PLATE-1 N: 1947", "U:U2 PI: PLATE-1 N: 1948", "U:U2 PI: PLATE-1 N: 1949", "U:U2 PI: PLATE-1 N: 1950", "U:U2 PI: PLATE-1 N: 1951", "U:U2 PI: PLATE-1 N: 1952", "U:U2 PI: PLATE-1 N: 1953", "U:U2 PI: PLATE-1 N: 1954", "U:U2 PI: PLATE-1 N: 1955", "U:U2 PI: PLATE-1 N: 1956", "U:U2 PI: PLATE-1 N: 1957", "U:U2 PI: PLATE-1 N: 15741", "U:U2 PI: PLATE-1 N: 15742", "U:U2 PI: PLATE-1 N: 15743", "U:U2 PI: PLATE-1 N: 15744", "U:U2 PI: PLATE-1 N: 15745", "U:U2 PI: PLATE-1 N: 15746", "U:U2 PI: PLATE-1 N: 15747", "U:U2 PI: PLATE-1 N: 15748", "U:U2 PI: PLATE-1 N: 15749", "U:U2 PI: PLATE-1 N: 15750", "U:U2 PI: PLATE-1 N: 15751", "U:U2 PI: PLATE-1 N: 15752", "U:U2 PI: PLATE-1 N: 15753", "U:U2 PI: PLATE-1 N: 15754", "U:U2 PI: PLATE-1 N: 15755", "U:U2 PI: PLATE-1 N: 15756", "U:U2 PI: PLATE-1 N: 15757", "U:U2 PI: PLATE-1 N: 15758", "U:U2 PI: PLATE-1 N: 15759", "U:U2 PI: PLATE-1 N: 15760", "U:U2 PI: PLATE-1 N: 15761", "U:U2 PI: PLATE-1 N: 15762", "U:U2 PI: PLATE-1 N: 15763", "U:U2 PI: PLATE-1 N: 15764", "U:U2 PI: PLATE-1 N: 15765", "U:U2 PI: PLATE-1 N: 15766", "U:U2 PI: PLATE-1 N: 15767", "U:U2 PI: PLATE-1 N: 15768", "U:U2 PI: PLATE-1 N: 15769", "U:U2 PI: PLATE-1 N: 15770", "U:U2 PI: PLATE-1 N: 15771", "U:U2 PI: PLATE-1 N: 15772", "U:U2 PI: PLATE-1 N: 15773", "U:U2 PI: PLATE-1 N: 15774", "U:U2 PI: PLATE-1 N: 15775", "U:U2 PI: PLATE-1 N: 15776", "U:U2 PI: PLATE-1 N: 15972", "U:U2 PI: PLATE-1 N: 15973", "U:U2 PI: PLATE-1 N: 15974", "U:U2 PI: PLATE-1 N: 15975", "U:U2 PI: PLATE-1 N: 15976", "U:U2 PI: PLATE-1 N: 15977", "U:U2 PI: PLATE-1 N: 15978", "U:U2 PI: PLATE-1 N: 15979", "U:U2 PI: PLATE-1 N: 15980", "U:U2 PI: PLATE-1 N: 15981", "U:U2 PI: PLATE-1 N: 15982", "U:U2 PI: PLATE-1 N: 15983", "U:U2 PI: PLATE-1 N: 15984", "U:U2 PI: PLATE-1 N: 15985", "U:U2 PI: PLATE-1 N: 15986", "U:U2 PI: PLATE-1 N: 15987", "U:U2 PI: PLATE-1 N: 15988", "U:U2 PI: PLATE-1 N: 15989", "U:U2 PI: PLATE-1 N: 15990", "U:U2 PI: PLATE-1 N: 15991", "U:U2 PI: PLATE-1 N: 15992", "U:U2 PI: PLATE-1 N: 15993", "U:U2 PI: PLATE-1 N: 15994", "U:U2 PI: PLATE-1 N: 15995", "U:U2 PI: PLATE-1 N: 15996", "U:U2 PI: PLATE-1 N: 15997", "U:U2 PI: PLATE-1 N: 15998", "U:U2 PI: PLATE-1 N: 16038", "U:U2 PI: PLATE-1 N: 16039" ) )')
tmpName = xy136.name
session.xyDataObjects.changeKey(tmpName, 'displacement')
xy1 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1']
xy2 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 3']
xy3 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 34']
xy4 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 35']
xy5 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 36']
xy6 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 37']
xy7 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 38']
xy8 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 39']
xy9 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 40']
xy10 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 41']
xy11 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 42']
xy12 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 43']
xy13 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 44']
xy14 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 45']
xy15 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 479']
xy16 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 481']
xy17 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 512']
xy18 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 513']
xy19 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 514']
xy20 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 515']
xy21 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 516']
xy22 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 517']
xy23 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 518']
xy24 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 519']
xy25 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 520']
xy26 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 521']
xy27 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 522']
xy28 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 523']
xy29 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 957']
xy30 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 959']
xy31 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 990']
xy32 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 991']
xy33 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 992']
xy34 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 993']
xy35 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 994']
xy36 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 995']
xy37 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 996']
xy38 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 997']
xy39 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 998']
xy40 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 999']
xy41 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1000']
xy42 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1001']
xy43 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1435']
xy44 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1437']
xy45 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1468']
xy46 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1469']
xy47 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1470']
xy48 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1471']
xy49 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1472']
xy50 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1473']
xy51 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1474']
xy52 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1475']
xy53 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1476']
xy54 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1477']
xy55 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1478']
xy56 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1479']
xy57 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1913']
xy58 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1915']
xy59 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1946']
xy60 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1947']
xy61 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1948']
xy62 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1949']
xy63 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1950']
xy64 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1951']
xy65 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1952']
xy66 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1953']
xy67 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1954']
xy68 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1955']
xy69 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1956']
xy70 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 1957']
xy71 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15741']
xy72 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15742']
xy73 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15743']
xy74 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15744']
xy75 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15745']
xy76 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15746']
xy77 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15747']
xy78 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15748']
xy79 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15749']
xy80 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15750']
xy81 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15751']
xy82 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15752']
xy83 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15753']
xy84 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15754']
xy85 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15755']
xy86 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15756']
xy87 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15757']
xy88 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15758']
xy89 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15759']
xy90 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15760']
xy91 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15761']
xy92 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15762']
xy93 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15763']
xy94 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15764']
xy95 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15765']
xy96 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15766']
xy97 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15767']
xy98 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15768']
xy99 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15769']
xy100 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15770']
xy101 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15771']
xy102 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15772']
xy103 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15773']
xy104 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15774']
xy105 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15775']
xy106 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15776']
xy107 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15972']
xy108 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15973']
xy109 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15974']
xy110 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15975']
xy111 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15976']
xy112 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15977']
xy113 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15978']
xy114 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15979']
xy115 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15980']
xy116 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15981']
xy117 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15982']
xy118 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15983']
xy119 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15984']
xy120 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15985']
xy121 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15986']
xy122 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15987']
xy123 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15988']
xy124 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15989']
xy125 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15990']
xy126 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15991']
xy127 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15992']
xy128 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15993']
xy129 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15994']
xy130 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15995']
xy131 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15996']
xy132 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15997']
xy133 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 15998']
xy134 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 16038']
xy135 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 16039']
xy136 = sum((xy1, xy2, xy3, xy4, xy5, xy6, xy7, xy8, xy9, xy10, xy11, xy12, 
    xy13, xy14, xy15, xy16, xy17, xy18, xy19, xy20, xy21, xy22, xy23, xy24, 
    xy25, xy26, xy27, xy28, xy29, xy30, xy31, xy32, xy33, xy34, xy35, xy36, 
    xy37, xy38, xy39, xy40, xy41, xy42, xy43, xy44, xy45, xy46, xy47, xy48, 
    xy49, xy50, xy51, xy52, xy53, xy54, xy55, xy56, xy57, xy58, xy59, xy60, 
    xy61, xy62, xy63, xy64, xy65, xy66, xy67, xy68, xy69, xy70, xy71, xy72, 
    xy73, xy74, xy75, xy76, xy77, xy78, xy79, xy80, xy81, xy82, xy83, xy84, 
    xy85, xy86, xy87, xy88, xy89, xy90, xy91, xy92, xy93, xy94, xy95, xy96, 
    xy97, xy98, xy99, xy100, xy101, xy102, xy103, xy104, xy105, xy106, xy107, 
    xy108, xy109, xy110, xy111, xy112, xy113, xy114, xy115, xy116, xy117, 
    xy118, xy119, xy120, xy121, xy122, xy123, xy124, xy125, xy126, xy127, 
    xy128, xy129, xy130, xy131, xy132, xy133, xy134, xy135))
xy136.setValues(
    sourceDescription='sum ( ( "RF:RF2 PI: PLATE-1 N: 1", "RF:RF2 PI: PLATE-1 N: 3", "RF:RF2 PI: PLATE-1 N: 34", "RF:RF2 PI: PLATE-1 N: 35", "RF:RF2 PI: PLATE-1 N: 36", "RF:RF2 PI: PLATE-1 N: 37", "RF:RF2 PI: PLATE-1 N: 38", "RF:RF2 PI: PLATE-1 N: 39", "RF:RF2 PI: PLATE-1 N: 40", "RF:RF2 PI: PLATE-1 N: 41", "RF:RF2 PI: PLATE-1 N: 42", "RF:RF2 PI: PLATE-1 N: 43", "RF:RF2 PI: PLATE-1 N: 44", "RF:RF2 PI: PLATE-1 N: 45", "RF:RF2 PI: PLATE-1 N: 479", "RF:RF2 PI: PLATE-1 N: 481", "RF:RF2 PI: PLATE-1 N: 512", "RF:RF2 PI: PLATE-1 N: 513", "RF:RF2 PI: PLATE-1 N: 514", "RF:RF2 PI: PLATE-1 N: 515", "RF:RF2 PI: PLATE-1 N: 516", "RF:RF2 PI: PLATE-1 N: 517", "RF:RF2 PI: PLATE-1 N: 518", "RF:RF2 PI: PLATE-1 N: 519", "RF:RF2 PI: PLATE-1 N: 520", "RF:RF2 PI: PLATE-1 N: 521", "RF:RF2 PI: PLATE-1 N: 522", "RF:RF2 PI: PLATE-1 N: 523", "RF:RF2 PI: PLATE-1 N: 957", "RF:RF2 PI: PLATE-1 N: 959", "RF:RF2 PI: PLATE-1 N: 990", "RF:RF2 PI: PLATE-1 N: 991", "RF:RF2 PI: PLATE-1 N: 992", "RF:RF2 PI: PLATE-1 N: 993", "RF:RF2 PI: PLATE-1 N: 994", "RF:RF2 PI: PLATE-1 N: 995", "RF:RF2 PI: PLATE-1 N: 996", "RF:RF2 PI: PLATE-1 N: 997", "RF:RF2 PI: PLATE-1 N: 998", "RF:RF2 PI: PLATE-1 N: 999", "RF:RF2 PI: PLATE-1 N: 1000", "RF:RF2 PI: PLATE-1 N: 1001", "RF:RF2 PI: PLATE-1 N: 1435", "RF:RF2 PI: PLATE-1 N: 1437", "RF:RF2 PI: PLATE-1 N: 1468", "RF:RF2 PI: PLATE-1 N: 1469", "RF:RF2 PI: PLATE-1 N: 1470", "RF:RF2 PI: PLATE-1 N: 1471", "RF:RF2 PI: PLATE-1 N: 1472", "RF:RF2 PI: PLATE-1 N: 1473", "RF:RF2 PI: PLATE-1 N: 1474", "RF:RF2 PI: PLATE-1 N: 1475", "RF:RF2 PI: PLATE-1 N: 1476", "RF:RF2 PI: PLATE-1 N: 1477", "RF:RF2 PI: PLATE-1 N: 1478", "RF:RF2 PI: PLATE-1 N: 1479", "RF:RF2 PI: PLATE-1 N: 1913", "RF:RF2 PI: PLATE-1 N: 1915", "RF:RF2 PI: PLATE-1 N: 1946", "RF:RF2 PI: PLATE-1 N: 1947", "RF:RF2 PI: PLATE-1 N: 1948", "RF:RF2 PI: PLATE-1 N: 1949", "RF:RF2 PI: PLATE-1 N: 1950", "RF:RF2 PI: PLATE-1 N: 1951", "RF:RF2 PI: PLATE-1 N: 1952", "RF:RF2 PI: PLATE-1 N: 1953", "RF:RF2 PI: PLATE-1 N: 1954", "RF:RF2 PI: PLATE-1 N: 1955", "RF:RF2 PI: PLATE-1 N: 1956", "RF:RF2 PI: PLATE-1 N: 1957", "RF:RF2 PI: PLATE-1 N: 15741", "RF:RF2 PI: PLATE-1 N: 15742", "RF:RF2 PI: PLATE-1 N: 15743", "RF:RF2 PI: PLATE-1 N: 15744", "RF:RF2 PI: PLATE-1 N: 15745", "RF:RF2 PI: PLATE-1 N: 15746", "RF:RF2 PI: PLATE-1 N: 15747", "RF:RF2 PI: PLATE-1 N: 15748", "RF:RF2 PI: PLATE-1 N: 15749", "RF:RF2 PI: PLATE-1 N: 15750", "RF:RF2 PI: PLATE-1 N: 15751", "RF:RF2 PI: PLATE-1 N: 15752", "RF:RF2 PI: PLATE-1 N: 15753", "RF:RF2 PI: PLATE-1 N: 15754", "RF:RF2 PI: PLATE-1 N: 15755", "RF:RF2 PI: PLATE-1 N: 15756", "RF:RF2 PI: PLATE-1 N: 15757", "RF:RF2 PI: PLATE-1 N: 15758", "RF:RF2 PI: PLATE-1 N: 15759", "RF:RF2 PI: PLATE-1 N: 15760", "RF:RF2 PI: PLATE-1 N: 15761", "RF:RF2 PI: PLATE-1 N: 15762", "RF:RF2 PI: PLATE-1 N: 15763", "RF:RF2 PI: PLATE-1 N: 15764", "RF:RF2 PI: PLATE-1 N: 15765", "RF:RF2 PI: PLATE-1 N: 15766", "RF:RF2 PI: PLATE-1 N: 15767", "RF:RF2 PI: PLATE-1 N: 15768", "RF:RF2 PI: PLATE-1 N: 15769", "RF:RF2 PI: PLATE-1 N: 15770", "RF:RF2 PI: PLATE-1 N: 15771", "RF:RF2 PI: PLATE-1 N: 15772", "RF:RF2 PI: PLATE-1 N: 15773", "RF:RF2 PI: PLATE-1 N: 15774", "RF:RF2 PI: PLATE-1 N: 15775", "RF:RF2 PI: PLATE-1 N: 15776", "RF:RF2 PI: PLATE-1 N: 15972", "RF:RF2 PI: PLATE-1 N: 15973", "RF:RF2 PI: PLATE-1 N: 15974", "RF:RF2 PI: PLATE-1 N: 15975", "RF:RF2 PI: PLATE-1 N: 15976", "RF:RF2 PI: PLATE-1 N: 15977", "RF:RF2 PI: PLATE-1 N: 15978", "RF:RF2 PI: PLATE-1 N: 15979", "RF:RF2 PI: PLATE-1 N: 15980", "RF:RF2 PI: PLATE-1 N: 15981", "RF:RF2 PI: PLATE-1 N: 15982", "RF:RF2 PI: PLATE-1 N: 15983", "RF:RF2 PI: PLATE-1 N: 15984", "RF:RF2 PI: PLATE-1 N: 15985", "RF:RF2 PI: PLATE-1 N: 15986", "RF:RF2 PI: PLATE-1 N: 15987", "RF:RF2 PI: PLATE-1 N: 15988", "RF:RF2 PI: PLATE-1 N: 15989", "RF:RF2 PI: PLATE-1 N: 15990", "RF:RF2 PI: PLATE-1 N: 15991", "RF:RF2 PI: PLATE-1 N: 15992", "RF:RF2 PI: PLATE-1 N: 15993", "RF:RF2 PI: PLATE-1 N: 15994", "RF:RF2 PI: PLATE-1 N: 15995", "RF:RF2 PI: PLATE-1 N: 15996", "RF:RF2 PI: PLATE-1 N: 15997", "RF:RF2 PI: PLATE-1 N: 15998", "RF:RF2 PI: PLATE-1 N: 16038", "RF:RF2 PI: PLATE-1 N: 16039" ) )')
tmpName = xy136.name
session.xyDataObjects.changeKey(tmpName, 'force')
xy1 = session.xyDataObjects['displacement']
xy2 = session.xyDataObjects['force']
xy3 = combine(xy1, xy2)
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
c1 = session.Curve(xyData=xy3)
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
#* There is no valid step data available on the database. If the analysis is 
#* running, the database must be closed and reopened once the results have been 
#* initialized. The requested operation has been cancelled.
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(width=0.00271254, 
    height=0.00107053, viewOffsetX=3.24919e-05, viewOffsetY=9.34282e-06)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR1_SIG11', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR2_SIG22', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR3_SIG33', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR4_SIG12', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR5_SIG13', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR6_SIG23', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR7_EPS11', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR8_EPS22', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR9_EPS33', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR10_EPS12', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR14_SIG_H', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00227702, 
    farPlane=0.00353903, width=0.000294536, height=0.000116242, 
    viewOffsetX=1.57775e-05, viewOffsetY=-2.66438e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR15_SIG_VONMISES', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR16_TRIAX', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.0022769, 
    farPlane=0.00353916, width=0.000333319, height=0.000131548, 
    viewOffsetX=1.93704e-05, viewOffsetY=-2.59464e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR19_HISTORY', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00228613, 
    farPlane=0.00352992, width=0.000253382, height=0.0001, 
    viewOffsetX=3.43937e-06, viewOffsetY=-2.14327e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR21_THETA_COVERAGE', outputPosition=INTEGRATION_POINT, 
    )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00221983, 
    farPlane=0.00359623, width=0.00064741, height=0.000255507, 
    viewOffsetX=4.18221e-05, viewOffsetY=-3.58851e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM1', outputPosition=INTEGRATION_POINT, )
session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=19 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00207896, 
    farPlane=0.00373715, width=0.00152068, height=0.000600154, 
    viewOffsetX=1.63765e-05, viewOffsetY=9.37564e-05)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00204707, 
    farPlane=0.00376901, width=0.00191785, height=0.000756899, 
    viewOffsetX=0.000116418, viewOffsetY=0.000100239)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR22_C_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00227405, 
    farPlane=0.00354202, width=0.000312927, height=0.0001235, 
    viewOffsetX=-3.78573e-05, viewOffsetY=-1.64173e-05)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM2', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=11 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=12 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=13 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=15 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=16 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=66 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=66 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM1', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM2', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00228608, 
    farPlane=0.00352993, width=0.000253382, height=0.0001, 
    viewOffsetX=-1.07908e-05, viewOffsetY=4.95798e-06)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UVARM3', outputPosition=INTEGRATION_POINT, )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR22_C_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=69 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR1_SIG11', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR13_EQPLAS', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR15_SIG_VONMISES', outputPosition=INTEGRATION_POINT, )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR4_SIG12', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR22_C_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR15_SIG_VONMISES', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00227534, 
    farPlane=0.00354065, width=0.000343634, height=0.000135619, 
    viewOffsetX=-1.23165e-06, viewOffsetY=-6.28008e-06)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR1_SIG11', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00228614, 
    farPlane=0.00352992, width=0.000253382, height=0.0001, 
    viewOffsetX=-2.61273e-05, viewOffsetY=6.35483e-06)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=11 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=12 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=13 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=15 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=16 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=17 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=18 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=19 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=20 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=21 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=21 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR15_SIG_VONMISES', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=11 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=12 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=13 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=15 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=16 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=17 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=18 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=19 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=20 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=21 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=22 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00225724, 
    farPlane=0.00355872, width=0.000415805, height=0.000164102, 
    viewOffsetX=4.97347e-06, viewOffsetY=1.63547e-05)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR16_TRIAX', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00219029, 
    farPlane=0.00362575, width=0.000829997, height=0.000327568, 
    viewOffsetX=2.11604e-05, viewOffsetY=-1.21129e-05)
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00229184, 
    farPlane=0.00351872, width=0.000868477, height=0.000342754, 
    cameraPosition=(0.00170849, 0.001106, 0.00208364), cameraUpVector=(
    -0.361468, 0.743515, -0.562607), cameraTarget=(2.70215e-05, 2.80906e-06, 
    -1.69113e-05), viewOffsetX=2.21414e-05, viewOffsetY=-1.26745e-05)
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.0023805, 
    farPlane=0.00343007, width=0.000253382, height=0.0001, 
    viewOffsetX=-1.05306e-05, viewOffsetY=-1.49324e-05)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00235005, 
    farPlane=0.0034605, width=0.000425964, height=0.000168111, 
    viewOffsetX=-3.21734e-05, viewOffsetY=-5.30754e-06)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR19_HISTORY', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00238048, 
    farPlane=0.00343007, width=0.000253382, height=0.0001, 
    viewOffsetX=-2.07933e-05, viewOffsetY=-7.7534e-06)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR20_GC', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00231565, 
    farPlane=0.0034949, width=0.000638255, height=0.000251894, 
    viewOffsetX=-4.68886e-05, viewOffsetY=5.49458e-06)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR21_THETA_COVERAGE', outputPosition=INTEGRATION_POINT, 
    )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00234985, 
    farPlane=0.0034607, width=0.000427168, height=0.000168587, 
    viewOffsetX=-5.49251e-05, viewOffsetY=-1.36697e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00233554, 
    farPlane=0.00347501, width=0.000580871, height=0.000229247, 
    viewOffsetX=-6.42119e-05, viewOffsetY=-9.63186e-06)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00227046, 
    farPlane=0.00354536, width=0.000376162, height=0.000148457, 
    viewOffsetX=-0.000258927, viewOffsetY=0.000137208)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.0019099, 
    farPlane=0.00390614, width=0.00243811, height=0.000962227, 
    viewOffsetX=7.73267e-05, viewOffsetY=-1.33443e-05)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00228608, 
    farPlane=0.00352993, width=0.000253382, height=0.0001, 
    viewOffsetX=-7.60121e-05, viewOffsetY=6.98539e-06)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00228614, 
    farPlane=0.00352992, width=0.000253382, height=0.0001, 
    viewOffsetX=-3.01091e-05, viewOffsetY=-2.19972e-05)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].view.setValues(width=0.000253382, 
    height=0.0001, viewOffsetX=-3.04857e-05, viewOffsetY=-9.33785e-06)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=11 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=12 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=13 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00217771, 
    farPlane=0.00363829, width=0.000907593, height=0.000358191, 
    viewOffsetX=2.11391e-05, viewOffsetY=6.05554e-06)
session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR14_SIG_H', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.002286, 
    farPlane=0.00352994, width=0.000253382, height=0.0001, 
    viewOffsetX=1.26317e-05, viewOffsetY=-1.41364e-05)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR16_TRIAX', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00228613, 
    farPlane=0.00352992, width=0.000253382, height=0.0001, 
    viewOffsetX=-6.88237e-05, viewOffsetY=1.29871e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(height=0.0001, 
    viewOffsetX=-8.13437e-05, viewOffsetY=-3.06374e-06)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR16_TRIAX', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00207341, 
    farPlane=0.00374264, width=0.00153376, height=0.000605315, 
    viewOffsetX=-5.57021e-05, viewOffsetY=-5.40045e-05)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00225119, 
    farPlane=0.00356487, width=0.000546487, height=0.000215677, 
    viewOffsetX=-0.000288668, viewOffsetY=5.88179e-05)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
#: The selected probe values were written to file "C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/abaqus.rpt".
#: 
#: ODB: phase_field_plate_3D_UEL.odb
#:    Number of nodes: 26230
#:    Number of elements: 20596
#:    Element types: C3D8 RNODE3D 
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.001993, 
    farPlane=0.00382306, width=0.00228806, height=0.000903009, 
    viewOffsetX=6.15423e-05, viewOffsetY=5.45453e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR13_EQPLAS', outputPosition=INTEGRATION_POINT, )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00206112, 
    farPlane=0.00375494, width=0.00163241, height=0.00064425, 
    viewOffsetX=4.25277e-05, viewOffsetY=3.46879e-05)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR16_TRIAX', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.0022358, 
    farPlane=0.00358024, width=0.000618038, height=0.000243915, 
    viewOffsetX=-9.22391e-05, viewOffsetY=1.33271e-05)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=11 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=12 )
session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00207816, 
    farPlane=0.00373791, width=0.00152631, height=0.000602376, 
    viewOffsetX=-0.000127038, viewOffsetY=0.000186819)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
#: Warning: The output database 'C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
leaf = dgo.LeafFromElementSets(elementSets=("PLATE-1.SOLID", ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.remove(leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.0018766, 
    farPlane=0.00393945, width=0.00272226, height=0.00107437, 
    viewOffsetX=0.000237378, viewOffsetY=-6.1051e-05)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=8 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=9 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=10 )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00209434, 
    farPlane=0.00372167, width=0.00142524, height=0.000562487, 
    viewOffsetX=-2.09591e-05, viewOffsetY=0.000269539)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=11 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=12 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=13 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=14 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=15 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=15 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=15 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=15 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=15 )
session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.odbs['C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb'].close(
    )
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       9
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR21_THETA_COVERAGE', outputPosition=INTEGRATION_POINT, 
    )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00192798, 
    farPlane=0.00388807, width=0.00247742, height=0.00097774, 
    viewOffsetX=0.000143128, viewOffsetY=5.53541e-05)
leaf = dgo.LeafFromElementSets(elementSets=("PLATE-1.SOLID", ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.remove(leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00200102, 
    farPlane=0.00381503, width=0.00223634, height=0.000882595, 
    viewOffsetX=-0.000186134, viewOffsetY=5.89706e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR6_SIG23', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR4_SIG12', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00201652, 
    farPlane=0.00379953, width=0.00212711, height=0.000839488, 
    viewOffsetX=-0.000189214, viewOffsetY=3.49317e-05)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR13_EQPLAS', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR14_SIG_H', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00228005, 
    farPlane=0.00353599, width=0.000275808, height=0.000108851, 
    viewOffsetX=-4.04298e-05, viewOffsetY=-1.99422e-06)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR15_SIG_VONMISES', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR16_TRIAX', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00221299, 
    farPlane=0.00360306, width=0.000775866, height=0.000306204, 
    viewOffsetX=-3.32113e-05, viewOffsetY=-2.90325e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR17_LODE', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00216476, 
    farPlane=0.00365129, width=0.00110917, height=0.000437745, 
    viewOffsetX=4.34119e-05, viewOffsetY=-7.68673e-05)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00223986, 
    farPlane=0.00357619, width=0.000592761, height=0.00023394, 
    viewOffsetX=-1.29864e-05, viewOffsetY=1.41931e-08)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00210371, 
    farPlane=0.00371234, width=0.00153009, height=0.000603865, 
    viewOffsetX=6.79351e-05, viewOffsetY=-3.21503e-06)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR18_PHI', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.0022852, 
    farPlane=0.00353085, width=0.000253382, height=0.0001, 
    viewOffsetX=-3.60961e-05, viewOffsetY=-9.66888e-06)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='SDV_AR23_CL_MOL', outputPosition=INTEGRATION_POINT, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=0.00197393, 
    farPlane=0.00384212, width=0.00241817, height=0.000954356, 
    viewOffsetX=-9.7984e-05, viewOffsetY=-5.48338e-06)
o1 = session.openOdb(
    name='C:/LocalUserData/User-data/nguyenb5/debug_phase_field_meter_many_SDV/phase_field_plate_3D_UEL.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
