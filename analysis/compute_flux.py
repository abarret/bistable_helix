import sys
sys.path.append("/scratch/sfw/visit/3.1.4/3.1.4/linux-x86_64/lib/site-packages")
from visit import *
base_path = "/scratch/bistable_helix/data/giesekus"
#alphas = ["0.3"]
#Des = ["0.01", "0.1", "0.3", "2.0"]
alphas = ["0.01", "0.05", "0.1", "0.2", "0.3"]
Des = ["0.01", "0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"]
Launch()
#base_path = "/scratch/ibamr/objs-opt/examples/cut_cell/ex2/10_cycles/both_leaflets/1000x"
DefineScalarExpression("udotn", "U[2]")
for alpha in alphas:
    for De in Des:
        full_path = base_path + "/alpha=" + alpha + "/De=" + De
        OpenDatabase("localhost:" + full_path + "/viz/dumps.visit", 0)
        DeleteAllPlots()

        AddPlot("Pseudocolor", "udotn", 1, 1)

        AddOperator("Clip", 0)
        clipAtts = ClipAttributes()
        clipAtts.funcType = clipAtts.Sphere
        clipAtts.center = (0.0, 0.0, 5.5)
        clipAtts.radius = 1.5
        clipAtts.sphereInverse = 1
        SetOperatorOptions(clipAtts, 0, 0)

        AddOperator("Slice", 0)
        sliceAtts = SliceAttributes()
        sliceAtts.originType = sliceAtts.Intercept
        sliceAtts.originIntercept = 5.5
        sliceAtts.axisType = sliceAtts.ZAxis
        sliceAtts.project2d = 0
        SetOperatorOptions(sliceAtts, 1, 0)

        DrawPlots()

        SetQueryOutputToValue()
        nts = TimeSliderGetNStates()
        filename = full_path + "/advective_flux_r_1.5.txt"
        file = open(filename, "w")
        for ts in range(nts):
            TimeSliderSetState(ts)
            time = Query("Time")
            print("Getting data at time " + str(time))
            flux = Query("Weighted Variable Sum")
            file.write(str(time) + " " + str(flux) + "\n")
        file.close()
        DeleteAllPlots()
        CloseDatabase("localhost:" + full_path + "/viz/dumps.visit")

Close()
