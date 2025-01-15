import sys
sys.path.append("/scratch/sfw/visit/3.1.4/3.1.4/linux-x86_64/lib/site-packages")
from visit import *
base_path = "/scratch/bistable_helix/data/giesekus"
alphas = ["0.05", "0.1", "0.2", "0.3"]
Des = ["0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"]
#alphas = ["0.3"]
#Des = ["1.5"]
Launch()
image_path = "/scratch/bistable_helix/data/vel/"

# Open Newtonian
OpenDatabase("localhost:/scratch/bistable_helix/data/De=0.0/viz/dumps.visit", 0)
for alpha in alphas:
    for De in Des:
        full_path = base_path + "/alpha=" + alpha + "/De=" + De
        DeleteAllPlots()

        OpenDatabase("localhost:/scratch/bistable_helix/data/De=0.0/viz/lag_data.visit", 0)
        AddPlot("Mesh", "Higdon_helix_101_mesh", 1, 0)
        meshAtts = MeshAttributes()
        meshAtts.legendFlag = 0
        meshAtts.lineWidth = 9
        SetPlotOptions(meshAtts)

        ActivateDatabase("localhost:/scratch/bistable_helix/data/De=0.0/viz/dumps.visit")

        DefineScalarExpression("vel_diff", "pos_cmfe(</scratch/bistable_helix/data/giesekus/alpha=" + alpha + "/De=" + De + "/viz/dumps.visit[0]id:U_magnitude>, <amr_mesh>, 0) - U_magnitude")

        AddPlot("Volume", "vel_diff", 1, 1)
        volAtts = VolumeAttributes()
        volAtts.legendFlag = 0
        volAtts.colorControlPoints.GetControlPoints(0).colors = (0, 0, 255, 255)
        volAtts.colorControlPoints.GetControlPoints(0).position = 0
        volAtts.colorControlPoints.GetControlPoints(1).colors = (255, 255, 255, 255)
        volAtts.colorControlPoints.GetControlPoints(1).position = 0.5
        volAtts.colorControlPoints.GetControlPoints(2).colors = (255, 0, 0, 255)
        volAtts.colorControlPoints.GetControlPoints(2).position = 1
        volAtts.colorControlPoints.GetControlPoints(3).colors = (255, 0, 0, 255)
        volAtts.colorControlPoints.GetControlPoints(3).position = 1
        volAtts.colorControlPoints.GetControlPoints(4).colors = (255, 255, 0, 255)
        volAtts.colorControlPoints.GetControlPoints(4).position = 1
        volAtts.opacityAttenuation = 1
        volAtts.opacityMode = volAtts.FreeformMode
        volAtts.freeformOpacity = (255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255)
        volAtts.useColorVarMin = 1
        volAtts.colorVarMin = -10
        volAtts.useColorVarMax = 1
        volAtts.colorVarMax = 10
        SetPlotOptions(volAtts)

        AddOperator("Box", 0)
        boxAtts = BoxAttributes()
        boxAtts.minx = -1
        boxAtts.maxx = 1
        boxAtts.miny = -1
        boxAtts.maxy = 1
        boxAtts.minz = 0
        boxAtts.maxz = 6
        SetOperatorOptions(boxAtts, 0, 0)

        DrawPlots()

        View3DAtts = View3DAttributes()
        View3DAtts.viewNormal = (0, 1, 0)
        View3DAtts.focus = (0, 0, 2)
        View3DAtts.viewUp = (0, 0, 1)
        View3DAtts.viewAngle = 30
        View3DAtts.parallelScale = 17.3205
        View3DAtts.nearPlane = -34.641
        View3DAtts.farPlane = 34.641
        View3DAtts.imagePan = (0, 0)
        View3DAtts.imageZoom = 4.59498
        View3DAtts.perspective = 1
        View3DAtts.eyeAngle = 2
        View3DAtts.centerOfRotationSet = 0
        View3DAtts.centerOfRotation = (0, 0, 2)
        View3DAtts.axis3DScaleFlag = 0
        View3DAtts.axis3DScales = (1, 1, 1)
        View3DAtts.shear = (0, 0, 1)
        View3DAtts.windowValid = 1
        SetView3D(View3DAtts)

        SetQueryOutputToValue()
        nts = TimeSliderGetNStates()
        TimeSliderSetState(159)
        saveAtts = SaveWindowAttributes()
        saveAtts.outputToCurrentDirectory = 0
        saveAtts.outputDirectory = image_path
        saveAtts.fileName = "alpha_" + str(alpha) + "_De_" + str(De)
        saveAtts.family = 0
        saveAtts.width = 4096
        SetSaveWindowAttributes(saveAtts)
        SaveWindow()

        DeleteAllPlots()
        CloseDatabase("localhost:" + full_path + "/viz/dumps.visit")

Close()
