# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#input
workingFolder='D:\\Kei_Nam\\Documents\\DGSolver\\DG2DSolver_newSendRecv\\CASES'
caseName='testDurst_2'
techplotCase='cylindercellCentered'
readByTechplotCaseName = False

# Stagnation line
Pt1 = [-0.11999999731779099, 2.081669990969516e-17, 0.0]
Pt2 = [0.10000000149011612, 2.081669990969516e-17, 0.0]

#Variables
dataArray=[]


#Doc time
def getTime():
	with open(workingFolder + '\\' + caseName + '\\' + 'time.txt') as f:
		array = []
		for line in f: # read rest of lines
			for x in line.split():
				array.append(x)
	#print(array[0])
	time=int(array[0])
	return time

#Doc processor
def getProcessor():
	totalProc=0
	with open(workingFolder + '\\' + caseName + '\\' + 'System\\DGOptions.txt') as f:
		for line in f: # read rest of lines
			array = []
			for x in line.split():
				array.append(x)
			if array[0] == 'totalProcess':
				totalProc=int(array[1])
	print(totalProc)
	return totalProc

var='T'

def readCase(inode):
	# create a new 'Tecplot Reader'
	if readByTechplotCaseName == True:
		techplotCellData = TecplotReader(FileNames=[workingFolder + '\\' + caseName + '\\Processor' + str(inode) + '\\TecplotFile\\' + str(time) + '\\' + techplotCase + '.dat'])
	else:
		techplotCellData = TecplotReader(FileNames=[workingFolder + '\\' + caseName + '\\Processor' + str(inode) + '\\TecplotFile\\' + str(time) + '\\' + caseName + 'cellCentered.dat'])
	
	dataArray.append(techplotCellData)

time=getTime()
#time=42000
totalProcesses=getProcessor()
for i in range(totalProcesses):
	readCase(i)



# Group data
groupDatasets1 = GroupDatasets(registrationName='GroupDatasets1', Input=dataArray)
# create a new 'Merge Blocks'
mergeBlocks1 = MergeBlocks(registrationName='MergeBlocks1', Input=groupDatasets1)
# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=mergeBlocks1)
cellDatatoPointData1.CellDataArraytoprocess = ['P', 'Rho', 'T', 'Theta1', 'Theta2', 'Velocity_mag', 'dRhoX', 'dRhoY', 'u', 'v']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1627, 821]

# show data in view
flatPlate2_12CorescellCentereddatDisplay = Show(cellDatatoPointData1, renderView1)

# trace defaults for the display properties.
flatPlate2_12CorescellCentereddatDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.03597879968583584, 0.012537300121039152, 10000.0]
renderView1.CameraFocalPoint = [0.03597879968583584, 0.012537300121039152, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(flatPlate2_12CorescellCentereddatDisplay, ('POINTS', var))

# rescale color and/or opacity maps used to include current data range
flatPlate2_12CorescellCentereddatDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
flatPlate2_12CorescellCentereddatDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'P'
pLUT = GetColorTransferFunction(var)

# get opacity transfer function/opacity map for 'P'
pPWF = GetOpacityTransferFunction(var)



# Plot stagnation line
# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=cellDatatoPointData1,
    Source='Line')

# init the 'Line' selected for 'Source'
plotOverLine1.Source.Point1 = Pt1
plotOverLine1.Source.Point2 = Pt2


# Plot on cylinder surface
# create a new 'Plot On Intersection Curves'
plotOnIntersectionCurves1 = PlotOnIntersectionCurves(registrationName='PlotOnIntersectionCurves1', Input=cellDatatoPointData1)
plotOnIntersectionCurves1.SliceType = 'Cylinder'

# Properties modified on plotOnIntersectionCurves1.SliceType
plotOnIntersectionCurves1.SliceType.Center = [0.2524, 0.0, 0.0]
plotOnIntersectionCurves1.SliceType.Axis = [0.0, 0.0, 1.0]
plotOnIntersectionCurves1.SliceType.Radius = 0.153

UpdatePipeline(time=0.0050001, proxy=plotOverLine1)
UpdatePipeline(time=0.0050001, proxy=plotOnIntersectionCurves1)