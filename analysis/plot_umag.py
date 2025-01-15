import sys
sys.path.append("/scratch/sfw/visit/3.1.4/3.1.4/linux-x86_64/lib/site-packages")
from visit import *
base_path = "/scratch/bistable_helix/data/giesekus"
alphas = ["0.01", "0.05", "0.1", "0.2", "0.3"]
Des = ["0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"]
#alphas = ["0.3"]
#Des = ["1.5"]
Launch()
image_path = "/scratch/bistable_helix/data/umag/"
DefineScalarExpression("umag", "magnitude(<U>)")
# Newtonian fluids
full_path = "/scratch/bistable_helix/data/De=0.0"
OpenDatabase("localhost:" + full_path + "/viz/dumps.visit", 0)
DeleteAllPlots()

AddPlot("Pseudocolor", "umag", 1, 1)
psAtts = PseudocolorAttributes()
psAtts.minFlag = 1
psAtts.min = 0
psAtts.maxFlag = 1
psAtts.max = 100
SetPlotOptions(psAtts)

AddOperator("Box", 0)
boxAtts = BoxAttributes()
boxAtts.minx = -2.5
boxAtts.maxx = 2.5
boxAtts.miny = -2.5
boxAtts.maxy = 2.5
boxAtts.minz = 0
boxAtts.maxz = 6
SetOperatorOptions(boxAtts, 0, 0)

AddOperator("Slice", 0)
slAtts = SliceAttributes()
slAtts.originType = slAtts.Intercept
slAtts.originIntercept = 3
slAtts.axisType = slAtts.ZAxis
slAtts.project2d = 0
SetOperatorOptions(slAtts, 2, 0)

OpenDatabase("localhost:" + full_path + "/viz/lag_data.visit", 0)
AddPlot("Mesh", "Higdon_helix_101_mesh", 1, 0)
meshAtts = MeshAttributes()
meshAtts.legendFlag = 0
meshAtts.lineWidth = 9
SetPlotOptions(meshAtts)

DrawPlots()

View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0, 0, 1)
View3DAtts.focus = (0, 0, 2)
View3DAtts.viewUp = (0, 1, 0)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 17.3205
View3DAtts.nearPlane = -34.641
View3DAtts.farPlane = 34.641
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 4.59497
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0, 0, 2)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)

annAtts = AnnotationAttributes()
annAtts.axes3D.visible = 0
annAtts.axes3D.triadFlag = 0
annAtts.userInfoFlag = 0
annAtts.databaseInfoFlag = 0
annAtts.legendInfoFlag = 0
SetAnnotationAttributes(annAtts)

SetQueryOutputToValue()
nts = TimeSliderGetNStates()
TimeSliderSetState(nts - 1)
saveAtts = SaveWindowAttributes()
saveAtts.outputToCurrentDirectory = 0
saveAtts.outputDirectory = image_path
saveAtts.fileName = "Nhigh"
saveAtts.family = 0
saveAtts.width = 2048
SetSaveWindowAttributes(saveAtts)
SaveWindow()

DeleteAllPlots()
CloseDatabase("localhost:" + full_path + "/viz/dumps.visit")

full_path = "/scratch/bistable_helix/data/De=0.0_scaled"
OpenDatabase("localhost:" + full_path + "/viz/dumps.visit", 0)
DeleteAllPlots()

AddPlot("Pseudocolor", "umag", 1, 1)
SetPlotOptions(psAtts)
AddOperator("Box", 0)
SetOperatorOptions(boxAtts, 0, 0)

AddOperator("Slice", 0)
SetOperatorOptions(slAtts, 2, 0)

OpenDatabase("localhost:" + full_path + "/viz/lag_data.visit", 0)
AddPlot("Mesh", "Higdon_helix_101_mesh", 1, 0)
SetPlotOptions(meshAtts)

DrawPlots()

SetView3D(View3DAtts)
SetAnnotationAttributes(annAtts)

SetQueryOutputToValue()
nts = TimeSliderGetNStates()
TimeSliderSetState(nts - 1)
saveAtts.fileName = "Nlow"
SetSaveWindowAttributes(saveAtts)
SaveWindow()

DeleteAllPlots()
CloseDatabase("localhost:" + full_path + "/viz/dumps.visit")

# Non newtonian fluids
for alpha in alphas:
    for De in Des:
        full_path = base_path + "/alpha=" + alpha + "/De=" + De
        OpenDatabase("localhost:" + full_path + "/viz/dumps.visit", 0)
        DeleteAllPlots()

        AddPlot("Pseudocolor", "umag", 1, 1)
        psAtts = PseudocolorAttributes()
        psAtts.minFlag = 1
        psAtts.min = 0
        psAtts.maxFlag = 1
        psAtts.max = 100
        SetPlotOptions(psAtts)

        AddOperator("Box", 0)
        boxAtts = BoxAttributes()
        boxAtts.minx = -2.5
        boxAtts.maxx = 2.5
        boxAtts.miny = -2.5
        boxAtts.maxy = 2.5
        boxAtts.minz = 0
        boxAtts.maxz = 6
        SetOperatorOptions(boxAtts, 0, 0)

        AddOperator("Slice", 0)
        slAtts = SliceAttributes()
        slAtts.originType = slAtts.Intercept
        slAtts.originIntercept = 3
        slAtts.axisType = slAtts.ZAxis
        slAtts.project2d = 0
        SetOperatorOptions(slAtts, 2, 0)

        OpenDatabase("localhost:" + full_path + "/viz/lag_data.visit", 0)
        AddPlot("Mesh", "Higdon_helix_101_mesh", 1, 0)
        meshAtts = MeshAttributes()
        meshAtts.legendFlag = 0
        meshAtts.lineWidth = 9
        SetPlotOptions(meshAtts)

        DrawPlots()

        View3DAtts = View3DAttributes()
        View3DAtts.viewNormal = (0, 0, 1)
        View3DAtts.focus = (0, 0, 2)
        View3DAtts.viewUp = (0, 1, 0)
        View3DAtts.viewAngle = 30
        View3DAtts.parallelScale = 17.3205
        View3DAtts.nearPlane = -34.641
        View3DAtts.farPlane = 34.641
        View3DAtts.imagePan = (0, 0)
        View3DAtts.imageZoom = 4.59497
        View3DAtts.perspective = 1
        View3DAtts.eyeAngle = 2
        View3DAtts.centerOfRotationSet = 0
        View3DAtts.centerOfRotation = (0, 0, 2)
        View3DAtts.axis3DScaleFlag = 0
        View3DAtts.axis3DScales = (1, 1, 1)
        View3DAtts.shear = (0, 0, 1)
        View3DAtts.windowValid = 1
        SetView3D(View3DAtts)

        annAtts = AnnotationAttributes()
        annAtts.axes3D.visible = 0
        annAtts.axes3D.triadFlag = 0
        annAtts.userInfoFlag = 0
        annAtts.databaseInfoFlag = 0
        annAtts.legendInfoFlag = 0
        SetAnnotationAttributes(annAtts)

        SetQueryOutputToValue()
        nts = TimeSliderGetNStates()
        TimeSliderSetState(nts - 1)
        saveAtts = SaveWindowAttributes()
        saveAtts.outputToCurrentDirectory = 0
        saveAtts.outputDirectory = image_path
        saveAtts.fileName = "alpha_" + str(alpha) + "_De_" + str(De)
        saveAtts.family = 0
        saveAtts.width = 2048
        SetSaveWindowAttributes(saveAtts)
        SaveWindow()

        DeleteAllPlots()
        CloseDatabase("localhost:" + full_path + "/viz/dumps.visit")

Close()
