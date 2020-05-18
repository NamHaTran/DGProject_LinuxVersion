# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#input
workingFolder='D:\\Kei_Nam\\Documents\\Private\\DG2DSolver_checkSlipBC\\CASES\\'
#caseName='ramp_Kn005'
caseName='w_Kn001_2'
#caseName='c14'
#caseName='c5'

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
	flatPlate2_12CorescellCentereddat = TecplotReader(FileNames=[workingFolder + '\\' + caseName + '\\Processor' + str(inode) + '\\TecplotFile\\' + str(time) + '\\' + caseName + 'cellCentered.dat'])

	# get active view
	renderView1 = GetActiveViewOrCreate('RenderView')
	# uncomment following to set a specific view size
	# renderView1.ViewSize = [1627, 821]

	# show data in view
	flatPlate2_12CorescellCentereddatDisplay = Show(flatPlate2_12CorescellCentereddat, renderView1)

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
	ColorBy(flatPlate2_12CorescellCentereddatDisplay, ('CELLS', var))

	# rescale color and/or opacity maps used to include current data range
	flatPlate2_12CorescellCentereddatDisplay.RescaleTransferFunctionToDataRange(True, False)

	# show color bar/color legend
	flatPlate2_12CorescellCentereddatDisplay.SetScalarBarVisibility(renderView1, True)

	# get color transfer function/color map for 'P'
	pLUT = GetColorTransferFunction(var)

	# get opacity transfer function/opacity map for 'P'
	pPWF = GetOpacityTransferFunction(var)

time=getTime()
totalProcesses=getProcessor()
for i in range(totalProcesses):
	readCase(i)