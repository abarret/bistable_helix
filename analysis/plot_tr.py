import sys
sys.path.append("/scratch/sfw/visit/3.1.4/3.1.4/linux-x86_64/lib/site-packages")
from visit import *
base_path = "/scratch/bistable_helix/data/giesekus"
alphas = ["0.05", "0.1", "0.2", "0.3"]
Des = ["0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"]
#alphas = ["0.3"]
#Des = ["1.5"]
Launch()
image_path = "/scratch/bistable_helix/data/tr/"
DefineScalarExpression("tr", "trace(Conformation_Tensor)")
for alpha in alphas:
    for De in Des:
        full_path = base_path + "/alpha=" + alpha + "/De=" + De
        OpenDatabase("localhost:" + full_path + "/viz/dumps.visit", 0)
        DeleteAllPlots()

        AddPlot("Pseudocolor", "tr", 1, 1)
        psAtts = PseudocolorAttributes()
        psAtts.scaling = psAtts.Log
        psAtts.minFlag = 1
        psAtts.min = 3
        psAtts.maxFlag = 1
        psAtts.max = 30
        SetPlotOptions(psAtts)

        AddOperator("Box", 0)
        boxAtts = BoxAttributes()
        boxAtts.minx = -1
        boxAtts.maxx = 1
        boxAtts.miny = -1
        boxAtts.maxy = 1
        boxAtts.minz = 0
        boxAtts.maxz = 6
        SetOperatorOptions(boxAtts, 0, 0)

        AddOperator("Resample", 0)
        reAtts = ResampleAttributes()
        reAtts.samplesX = 25
        reAtts.samplesY = 25
        reAtts.samplesZ = 25
        SetOperatorOptions(reAtts, 1, 0)

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
        meshAtts.lineWidth = 5
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
        View3DAtts.imageZoom = 11.9182
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
        TimeSliderSetState(nts - 1)
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
