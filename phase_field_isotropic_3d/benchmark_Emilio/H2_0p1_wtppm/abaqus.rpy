# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023.HF4 replay file
# Internal Version: 2023_07_21-20.45.57 RELr425 183702
# Run by nguyenb5 on Tue Aug 13 14:27:27 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=114.666664123535, 
    height=132.6875)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='phase_field_plate_uel_0p1_wtppm.odb')
#: Model: C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/benchmark_Emilio/H2_0p1_wtppm/phase_field_plate_uel_0p1_wtppm.odb
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
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(COMPONENT, 'RF2'), )
odb = session.odbs['C:/LocalUserData/User-data/nguyenb5/CP1000 (all combined)/phase_field_isotropic_3d/benchmark_Emilio/H2_0p1_wtppm/phase_field_plate_uel_0p1_wtppm.odb']
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
    NODAL, ((COMPONENT, 'RF2'), )), ('U', NODAL, ((COMPONENT, 'U2'), )), ), 
    nodePick=(('PLATE-1', 51, (
    '[#fffc0000 #1f #0:6 #80000000 #3f #0:163 #a000', 
    ' #0:6 #4 #20000010 #40004000 #44000004 #a600 #40000', 
    ' #6000000 #20000000 #0 #1 #8 #0:18 #3800', ' #0:132 #40000000 #0 #1000 ]', 
    )), ), )
xy1 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 19']
xy2 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 20']
xy3 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 21']
xy4 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 22']
xy5 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 23']
xy6 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 24']
xy7 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 25']
xy8 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 26']
xy9 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 27']
xy10 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 28']
xy11 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 29']
xy12 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 30']
xy13 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 31']
xy14 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 32']
xy15 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 33']
xy16 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 34']
xy17 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 35']
xy18 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 36']
xy19 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 37']
xy20 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 288']
xy21 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 289']
xy22 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 290']
xy23 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 291']
xy24 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 292']
xy25 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 293']
xy26 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 294']
xy27 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5550']
xy28 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5552']
xy29 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5763']
xy30 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5797']
xy31 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5822']
xy32 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5839']
xy33 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5855']
xy34 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5859']
xy35 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5883']
xy36 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5887']
xy37 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5898']
xy38 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5899']
xy39 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5902']
xy40 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5904']
xy41 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5939']
xy42 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5978']
xy43 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 5979']
xy44 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 6014']
xy45 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 6049']
xy46 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 6084']
xy47 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 6700']
xy48 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 6701']
xy49 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 6702']
xy50 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 10975']
xy51 = session.xyDataObjects['RF:RF2 PI: PLATE-1 N: 11021']
xy52 = sum((xy1, xy2, xy3, xy4, xy5, xy6, xy7, xy8, xy9, xy10, xy11, xy12, 
    xy13, xy14, xy15, xy16, xy17, xy18, xy19, xy20, xy21, xy22, xy23, xy24, 
    xy25, xy26, xy27, xy28, xy29, xy30, xy31, xy32, xy33, xy34, xy35, xy36, 
    xy37, xy38, xy39, xy40, xy41, xy42, xy43, xy44, xy45, xy46, xy47, xy48, 
    xy49, xy50, xy51))
xy52.setValues(
    sourceDescription='sum ( ( "RF:RF2 PI: PLATE-1 N: 19", "RF:RF2 PI: PLATE-1 N: 20", "RF:RF2 PI: PLATE-1 N: 21", "RF:RF2 PI: PLATE-1 N: 22", "RF:RF2 PI: PLATE-1 N: 23", "RF:RF2 PI: PLATE-1 N: 24", "RF:RF2 PI: PLATE-1 N: 25", "RF:RF2 PI: PLATE-1 N: 26", "RF:RF2 PI: PLATE-1 N: 27", "RF:RF2 PI: PLATE-1 N: 28", "RF:RF2 PI: PLATE-1 N: 29", "RF:RF2 PI: PLATE-1 N: 30", "RF:RF2 PI: PLATE-1 N: 31", "RF:RF2 PI: PLATE-1 N: 32", "RF:RF2 PI: PLATE-1 N: 33", "RF:RF2 PI: PLATE-1 N: 34", "RF:RF2 PI: PLATE-1 N: 35", "RF:RF2 PI: PLATE-1 N: 36", "RF:RF2 PI: PLATE-1 N: 37", "RF:RF2 PI: PLATE-1 N: 288", "RF:RF2 PI: PLATE-1 N: 289", "RF:RF2 PI: PLATE-1 N: 290", "RF:RF2 PI: PLATE-1 N: 291", "RF:RF2 PI: PLATE-1 N: 292", "RF:RF2 PI: PLATE-1 N: 293", "RF:RF2 PI: PLATE-1 N: 294", "RF:RF2 PI: PLATE-1 N: 5550", "RF:RF2 PI: PLATE-1 N: 5552", "RF:RF2 PI: PLATE-1 N: 5763", "RF:RF2 PI: PLATE-1 N: 5797", "RF:RF2 PI: PLATE-1 N: 5822", "RF:RF2 PI: PLATE-1 N: 5839", "RF:RF2 PI: PLATE-1 N: 5855", "RF:RF2 PI: PLATE-1 N: 5859", "RF:RF2 PI: PLATE-1 N: 5883", "RF:RF2 PI: PLATE-1 N: 5887", "RF:RF2 PI: PLATE-1 N: 5898", "RF:RF2 PI: PLATE-1 N: 5899", "RF:RF2 PI: PLATE-1 N: 5902", "RF:RF2 PI: PLATE-1 N: 5904", "RF:RF2 PI: PLATE-1 N: 5939", "RF:RF2 PI: PLATE-1 N: 5978", "RF:RF2 PI: PLATE-1 N: 5979", "RF:RF2 PI: PLATE-1 N: 6014", "RF:RF2 PI: PLATE-1 N: 6049", "RF:RF2 PI: PLATE-1 N: 6084", "RF:RF2 PI: PLATE-1 N: 6700", "RF:RF2 PI: PLATE-1 N: 6701", "RF:RF2 PI: PLATE-1 N: 6702", "RF:RF2 PI: PLATE-1 N: 10975", "RF:RF2 PI: PLATE-1 N: 11021" ) )')
tmpName = xy52.name
session.xyDataObjects.changeKey(tmpName, 'force')
xy1 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 19']
xy2 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 20']
xy3 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 21']
xy4 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 22']
xy5 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 23']
xy6 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 24']
xy7 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 25']
xy8 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 26']
xy9 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 27']
xy10 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 28']
xy11 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 29']
xy12 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 30']
xy13 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 31']
xy14 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 32']
xy15 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 33']
xy16 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 34']
xy17 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 35']
xy18 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 36']
xy19 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 37']
xy20 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 288']
xy21 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 289']
xy22 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 290']
xy23 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 291']
xy24 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 292']
xy25 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 293']
xy26 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 294']
xy27 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5550']
xy28 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5552']
xy29 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5763']
xy30 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5797']
xy31 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5822']
xy32 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5839']
xy33 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5855']
xy34 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5859']
xy35 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5883']
xy36 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5887']
xy37 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5898']
xy38 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5899']
xy39 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5902']
xy40 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5904']
xy41 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5939']
xy42 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5978']
xy43 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 5979']
xy44 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 6014']
xy45 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 6049']
xy46 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 6084']
xy47 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 6700']
xy48 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 6701']
xy49 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 6702']
xy50 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 10975']
xy51 = session.xyDataObjects['U:U2 PI: PLATE-1 N: 11021']
xy52 = avg((xy1, xy2, xy3, xy4, xy5, xy6, xy7, xy8, xy9, xy10, xy11, xy12, 
    xy13, xy14, xy15, xy16, xy17, xy18, xy19, xy20, xy21, xy22, xy23, xy24, 
    xy25, xy26, xy27, xy28, xy29, xy30, xy31, xy32, xy33, xy34, xy35, xy36, 
    xy37, xy38, xy39, xy40, xy41, xy42, xy43, xy44, xy45, xy46, xy47, xy48, 
    xy49, xy50, xy51))
xy52.setValues(
    sourceDescription='avg ( ( "U:U2 PI: PLATE-1 N: 19", "U:U2 PI: PLATE-1 N: 20", "U:U2 PI: PLATE-1 N: 21", "U:U2 PI: PLATE-1 N: 22", "U:U2 PI: PLATE-1 N: 23", "U:U2 PI: PLATE-1 N: 24", "U:U2 PI: PLATE-1 N: 25", "U:U2 PI: PLATE-1 N: 26", "U:U2 PI: PLATE-1 N: 27", "U:U2 PI: PLATE-1 N: 28", "U:U2 PI: PLATE-1 N: 29", "U:U2 PI: PLATE-1 N: 30", "U:U2 PI: PLATE-1 N: 31", "U:U2 PI: PLATE-1 N: 32", "U:U2 PI: PLATE-1 N: 33", "U:U2 PI: PLATE-1 N: 34", "U:U2 PI: PLATE-1 N: 35", "U:U2 PI: PLATE-1 N: 36", "U:U2 PI: PLATE-1 N: 37", "U:U2 PI: PLATE-1 N: 288", "U:U2 PI: PLATE-1 N: 289", "U:U2 PI: PLATE-1 N: 290", "U:U2 PI: PLATE-1 N: 291", "U:U2 PI: PLATE-1 N: 292", "U:U2 PI: PLATE-1 N: 293", "U:U2 PI: PLATE-1 N: 294", "U:U2 PI: PLATE-1 N: 5550", "U:U2 PI: PLATE-1 N: 5552", "U:U2 PI: PLATE-1 N: 5763", "U:U2 PI: PLATE-1 N: 5797", "U:U2 PI: PLATE-1 N: 5822", "U:U2 PI: PLATE-1 N: 5839", "U:U2 PI: PLATE-1 N: 5855", "U:U2 PI: PLATE-1 N: 5859", "U:U2 PI: PLATE-1 N: 5883", "U:U2 PI: PLATE-1 N: 5887", "U:U2 PI: PLATE-1 N: 5898", "U:U2 PI: PLATE-1 N: 5899", "U:U2 PI: PLATE-1 N: 5902", "U:U2 PI: PLATE-1 N: 5904", "U:U2 PI: PLATE-1 N: 5939", "U:U2 PI: PLATE-1 N: 5978", "U:U2 PI: PLATE-1 N: 5979", "U:U2 PI: PLATE-1 N: 6014", "U:U2 PI: PLATE-1 N: 6049", "U:U2 PI: PLATE-1 N: 6084", "U:U2 PI: PLATE-1 N: 6700", "U:U2 PI: PLATE-1 N: 6701", "U:U2 PI: PLATE-1 N: 6702", "U:U2 PI: PLATE-1 N: 10975", "U:U2 PI: PLATE-1 N: 11021" ) )')
tmpName = xy52.name
session.xyDataObjects.changeKey(tmpName, 'displacement')
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
x0 = session.xyDataObjects['displacement']
x1 = session.xyDataObjects['force']
session.xyReportOptions.setValues(numDigits=9)
session.writeXYReport(fileName='FD_curve.txt_0p1_wtppm.txt', appendMode=OFF, 
    xyData=(x0, x1))
